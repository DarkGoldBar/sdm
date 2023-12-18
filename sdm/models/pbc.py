# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/13
from __future__ import annotations
from typing import Optional, Sequence, Tuple, Union, Iterator, List, Dict
import math
import itertools
import collections
import spglib
import numpy as np
import networkx as nx
from numpy import ndarray, dot, transpose
from ..utils import latt62Matrix, matrix2Latt6
from ..analysis.similarity import kabsch
from ..analysis.find_bonds import ddosap_py, find_empiric_bonds


class PBC():
    """周期性边界条件的基础功能类

    Properties:
        lattice: Lattice
    """
    def __init__(self, lattice: Union[ndarray, Lattice], symprec: float = 1e-4) -> None:
        if isinstance(lattice, Lattice):
            self.lattice = lattice
        else:
            self.lattice = Lattice(lattice)
        self._frac_coords = None
        self.symprec = symprec

    @property
    def frac_coords(self) -> ndarray:
        """切片赋值将会丢失！"""
        self._frac_coords = self.lattice.getFractionalCoords(self.coords)
        return self._frac_coords

    @frac_coords.setter
    def frac_coords(self, value: ndarray):
        self._frac_coords = value
        self.coords = self.lattice.getCartesianCoords(self._frac_coords)

    @property
    def volume(self) -> float:
        return self.lattice.volume

    def moveAtomsInCell(self) -> None:
        """把原子移动到胞内"""
        self.frac_coords -= np.floor(self.frac_coords)

    def moveMolCenterInCell(self, groups=None) -> None:
        """把分子的几何中心移动到胞内"""
        groups = (self.groups if groups is None else groups)
        fcoords = self.frac_coords
        for g in groups:
            offset = np.floor(np.average(fcoords[list(g)], axis=0))
            fcoords[list(g)] -= offset
        self.frac_coords = fcoords

    def moveAtomsByMolecule(self, detect_cut=1.25, symprec=1e-5) -> None:
        radius = np.array([at.covelent_radius for at in self.elements], float)
        affine_matrix = np.eye(4).reshape((1, 4, 4))
        col, vec = ddosap_py(self.frac_coords, radius * detect_cut,
                             self.lattice.parameters,
                             affine_matrix,
                             symprec=symprec)

        bonds = [b[:2] for b in find_empiric_bonds(self.coords, radius, detect_cut)]
        mask = (col[:, 0] < col[:, 1])
        pbe_bonds = [tuple(b) for b in col[mask]]
        LTlattice = Lattice.fromParameters(*self.lattice.parameters)
        pbe_vec = vec[mask] @ (LTlattice.inv_matrix @ self.lattice.matrix)

        G = nx.Graph()
        G.add_nodes_from(range(len(self)))
        G.add_edges_from(bonds)
        group = None
        for idx, bd in enumerate(pbe_bonds):
            if bd not in G.edges:
                for group in nx.connected_components(G):
                    if bd[1] in group:
                        break
                vector = self.coords[bd[0]] - self.coords[bd[1]] + pbe_vec[idx]
                self.coords[list(group)] += vector
                bonds = [b[:2] for b in find_empiric_bonds(
                    self.coords, radius, detect_cut)]
                G.add_edges_from(bonds)

    def setLatticeLowerTriangular(self) -> None:
        fcoords = self.frac_coords
        self.lattice = Lattice.fromParameters(*self.lattice.parameters)
        self.frac_coords = fcoords

    def setLattice(self, matrix: ndarray, fix_frac=False,
                   translate=False, rotate=False) -> None:
        """重建晶胞

        默认保持笛卡尔坐标

        fix_frac: 保持分数坐标
        translate: 平移分子
        rotate: 旋转分子
        """
        groups = self.groups
        coords = None if self.coords is None else self.coords
        fcoords = self.frac_coords
        geom_cent1 = [self.getGeometryCenter(g) for g in groups]
        # Update lattice
        self.lattice = Lattice(matrix)

        if fix_frac or translate or rotate:
            self.frac_coords = fcoords
        if translate or rotate:
            geom_cent2 = [self.getGeometryCenter(g) for g in groups]
            U_mat = None
            if rotate:
                U_mat = [kabsch(coords[gp] - gc1, self.coords[gp] - gc2)
                         for gp, gc1, gc2 in zip(groups, geom_cent1, geom_cent2)]

            self.coords = coords

            if rotate:
                for gp, gc1, U in zip(groups, geom_cent1, U_mat):
                    self.rotateWithMatrix(U, center=gc1, group=gp)
            if translate:
                for gp, gc1, gc2 in zip(groups. geom_cent1, geom_cent2):
                    self.translate(gc2 - gc1, group=gp)

    def genKpoints(self, ktype: str, sym_kps=True, h_max=15., kspacing=0.5,
                   kmesh: Optional[Sequence[int]] = None,
                   kshift: Optional[Sequence[float]] = None) -> list:
        """生成K点

        Args:
            ktype : Kpoints type, enum:["vasp", "dftb"]
            h_max : max height, (ktype=='dftb')
            kspacing : kpoints spacing, (ktype=='vasp')
            kmesh & kshift : Specify k-mesh and shift
            sym_kps : If symmetrize kpoints

        Return:
            if sym_kps == False:
                List[int]
            if sym_kps == True:
                List[((kx, ky, kz), k_weight)]
        """
        # ---> mesh + shift
        if kmesh is not None:
            mesh, shift = kmesh, kshift
            if shift is None:
                shift = [(0.0 if x % 2 == 1 else 0.5) for x in mesh]
        elif ktype == 'dftb':
            height = self.lattice.height
            mesh = [int(x) for x in np.ceil(h_max / height)]
            shift = [(0.0 if x % 2 == 1 else 0.5) for x in mesh]
        elif ktype == 'vasp':
            rmat = self.lattice.reciprocal_lattice.matrix
            mesh = [max(1, int(x)) for x in np.ceil(np.linalg.norm(rmat, axis=1) / kspacing)]
            shift = [(0.0 if x % 2 == 1 else 0.5) for x in mesh]
        else:
            raise ValueError('Kpoints type must be vasp or dftb')

        # mesh + shift ---> rmesh
        if sym_kps:
            is_shift = np.array([1 if i else 0 for i in shift])
            mapping, grid = spglib.get_ir_reciprocal_mesh(
                np.array(mesh), self.toSpglibCell(), is_shift, symprec=self.symprec)
            rmesh = []
            tmp_map = list(mapping)
            for i in np.unique(mapping):
                rmesh.append(((grid[i] + is_shift * 0.5) / mesh, tmp_map.count(i)))
        else:
            rmesh = mesh
        return rmesh


class Lattice():
    def __init__(self, matrix: ndarray):
        matrix = getattr(matrix, 'matrix', matrix)
        self.matrix = np.array(matrix, dtype=float).reshape((3, 3))
        self.matrix.setflags(write=False)
        self._inv_matrix = None
        self._params = None
        self._height = None

    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix.copy())

    @classmethod
    def fromParameters(cls, a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
        return cls(latt62Matrix([a, b, c, alpha, beta, gamma]))

    @property
    def parameters(self) -> Tuple[float, float, float, float, float, float]:
        """
        Returns: (a, b, c, alpha, beta, gamma).
        """
        if self._params is None:
            self._params = matrix2Latt6(self.matrix)
        return self._params

    @property
    def abc(self) -> Tuple[float, float, float]:
        """
        :return: The lengths (a, b, c) of the lattice.
        """
        return self.parameters[:3]

    @property
    def length(self) -> Tuple[float, float, float]:
        """
        :return: The lengths (a, b, c) of the lattice.
        """
        return self.parameters[:3]

    @property
    def angles(self) -> Tuple[float, float, float]:
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        return self.parameters[3:]

    @property
    def inv_matrix(self) -> np.ndarray:
        """
        Inverse of lattice matrix.
        """
        if self._inv_matrix is None:
            self._inv_matrix = np.linalg.inv(self.matrix)
            self._inv_matrix.setflags(write=False)
        return self._inv_matrix

    @property
    def reciprocal_lattice(self) -> "Lattice":
        """
        Return the reciprocal lattice. Note that this is the standard
        reciprocal lattice used for solid state physics with a factor of 2 *
        pi. If you are looking for the crystallographic reciprocal lattice,
        use the reciprocal_lattice_crystallographic property.
        The property is lazily generated for efficiency.
        """
        v = np.linalg.inv(self.matrix).T
        return Lattice(v * 2 * np.pi)

    @property
    def reciprocal_lattice_crystallographic(self) -> "Lattice":
        """
        Returns the *crystallographic* reciprocal lattice, i.e., no factor of
        2 * pi.
        """
        return Lattice(self.reciprocal_lattice.matrix / (2 * np.pi))

    @property
    def volume(self) -> float:
        """
        返回体积
        """
        return np.dot(np.cross(self.matrix[0], self.matrix[1]), self.matrix[2])

    @property
    def height(self) -> float:
        """
        返回高
        """
        if self._height is None:
            a, b, c = self.abc
            bottom = (a * b * c) / np.array(self.abc) * np.sin(np.radians(self.angles))
            self._height = self.volume / bottom
        return self._height

    @property
    def metric_tensor(self) -> np.ndarray:
        """
        The metric tensor of the lattice.
        """
        return dot(self.matrix, self.matrix.T)

    def getCartesianCoords(self, fractional_coords: ndarray) -> ndarray:
        """从分数坐标计算直角坐标
        """
        return np.dot(fractional_coords, self.matrix)

    def getFractionalCoords(self, cart_coords: ndarray) -> ndarray:
        """从直角坐标计算分数坐标
        """
        return np.dot(cart_coords, self.inv_matrix)

    def getSupercellPacking(self,
                            min_size: Optional[float] = None,
                            max_volume: Optional[float] = None) -> List[int, int, int]:
        """计算合适的超胞堆积

        Args:
            min_size: 宽度的最小值
            max_volume: 体积的最大值, 选择接近正方的胞

        Return:
            packing: List[int, int, int]
        """
        packing = [1, 1, 1]

        if min_size is not None:
            packing = np.ceil(min_size / self.height).astype(int).tolist()

        if max_volume is not None:
            ncell = max_volume / self.volume
            packlist = []
            for i in range(packing[0], int(ncell) + 1):
                for j in range(packing[1], int(ncell / i) + 1):
                    k = int(ncell / i / j)
                    packlist.append([i, j, k])
            if packlist:
                stdh = np.std(packlist * self.height, axis=1)
                packing = packlist[np.argmin(stdh)]
            else:
                msg = "PBC: cannot find packing with volume < {}, current {}"
                msg = msg.format(max_volume, self.volume * packing[0] * packing[1] * packing[2])
                print(msg)
        return packing

    def getLatticeSystem(self, int_number=None, atol=0.1) -> str:
        """判断当前结构的晶系
        """
        if not int_number:
            if hasattr(self, "space_group"):
                int_number = self.space_group.int_number
            else:
                raise ValueError("无法读取 int_number")
        if not 1 <= int_number <= 230:
            raise ValueError("错误的 int_number: {}".format(int_number))

        # 大部分晶系 == 晶格
        system_list = [
            (2, "Triclinic"),
            (15, "Monoclinic"),
            (74, "Orthorhombic"),
            (142, "Tetragonal"),
            (167, "Trigonal"),
            (194, "Hexagonal"),
            (230, "Cubic"),
        ]
        for k, v in system_list:
            lattice_system = v
            if int_number <= k:
                break

        # 但是, 三方晶系 => 六方晶格 / 棱方晶格
        if lattice_system == 'Trigonal':
            angles = np.array(self.angles)
            if np.allclose(angles, angles[[1, 2, 0]], atol=atol):
                lattice_system = 'Rhombohedral'
            else:
                lattice_system = 'Hexagonal'
        return lattice_system

    def getLatticeSystemByParameters(self, atol=0.01) -> str:
        """从晶格参数判断晶系
        """
        from numpy import allclose, isclose
        param = np.array(self.parameters)
        lattice_system = "Triclinic"
        if allclose(param[3:], 90, atol=atol):
            # 立方
            if allclose(param[:3], param[0], atol=atol):
                lattice_system = "Cubic"
            # 正方: 3
            elif any(isclose(param[[0, 1, 2]], param[[1, 2, 0]], atol=atol)):
                lattice_system = "Tetragonal"
            # 正交
            else:
                lattice_system = "Orthorhombic"
        # 棱方
        elif allclose(param, param[[0, 0, 0, 3, 3, 3]], atol=atol):
            lattice_system = "Rhombohedral"
        elif sum(isclose(param[3:6], 90, atol=atol)) >= 2:
            # 六方
            if allclose(param[[1, 3]], [param[2], 120], atol=atol):
                lattice_system = "Hexagonal"
            elif allclose(param[[0, 4]], [param[2], 120], atol=atol):
                lattice_system = "Hexagonal"
            elif allclose(param[[0, 5]], [param[1], 120], atol=atol):
                lattice_system = "Hexagonal"
            # 单斜
            else:
                lattice_system = "Monoclinic"
        return lattice_system

    def setParametersAuto(self, int_number=None, atol=0.01) -> "Lattice":
        """根据晶系微调晶胞参数
        """
        if int_number is None:
            lattice_system = self.getLatticeSystemByParameters(atol=atol)
        else:
            lattice_system = self.getLatticeSystem(int_number, atol=atol)
        param = np.array(self.parameters)
        if lattice_system == "Monoclinic":
            k = np.argmax(np.abs(param[3:] - 90))
            idx = (np.arange(3) + k) % 3
            param[idx + 3] = [param[k+3], 90, 90]
        elif lattice_system == "Hexagonal":
            k = np.argmax(np.abs(param[3:] - 90))
            idx = (np.arange(3) + k) % 3
            param[idx + 3] = [120, 90, 90]
            param[idx[1:]] = np.average(param[idx[1:]])
        elif lattice_system == "Tetragonal":
            k = np.argmax(np.abs(param[:3] - np.average(param[:3])))
            idx = (np.arange(3) + k) % 3
            param[idx[1:]] = np.average(param[idx[1:]])
            param[3:] = 90
        elif lattice_system == "Cubic":
            param[:3] = np.average(param[:3])
            param[3:] = 90
        elif lattice_system == "Orthorhombic":
            param[3:] = 90
        elif lattice_system == "Rhombohedral":
            param[:3] = np.average(param[:3])
            param[3:] = np.average(param[3:])
        else:
            return self
        disparity = sum(abs(param - self.parameters))
        if disparity > atol:
            msg = "disparity too far:\n  {}\n  {}\n  {}".format(lattice_system, self.parameters, param)
            raise RuntimeError(msg)
        return self.setParameters(param)

    def setParameters(self, param: ndarray) -> "Lattice":
        """从param生成matrix, 并旋转至当前矩阵最近的角度
        """
        matrix = latt62Matrix(param)
        rot = kabsch(matrix, self.matrix)
        mat_final = np.dot(matrix, rot)
        return self.__class__(mat_final)

    def transformCentering(self, old='P', new='P') -> "Lattice":
        """转换晶胞中心类型, 更新晶胞

        old: 原晶胞中心类型
        new: 目标晶胞中心类型
        """
        transform_matrix = {
            "A": [1,  0,  0,  0,  1,  1,  0, -1,  1],
            "B": [1,  0, -1,  0,  1,  0,  1,  0,  1],
            "C": [1,  1,  0, -1,  1,  0,  0,  0,  1],
            "F": [-1, 1,  1,  1, -1,  1,  1,  1, -1],
            "I": [0,  1,  1,  1,  0,  1,  1,  1,  0],
            "H": [1,  0,  1, -1,  1,  1,  0, -1,  1],
        }
        old = old.upper()
        new = new.upper()
        if old in "PR" and new in "ABCFIH":
            tmat = np.array(transform_matrix.get(new)).reshape((3, 3)).T
        elif old in "ABCFIH" and new in "PR":
            tmat = np.array(transform_matrix.get(old)).reshape((3, 3)).T
            tmat = np.linalg.inv(tmat)
        else:
            raise ValueError("Not implemented!")
        return self.__class__(np.dot(tmat, self.matrix))

    def d_hkl(self, miller_index: ndarray) -> float:
        """
        Returns the distance between the hkl plane and the origin
        Args:
            miller_index ([h,k,l]): Miller index of plane
        Returns:
            d_hkl (float)
        """
        gstar = self.reciprocal_lattice_crystallographic.metric_tensor
        hkl = np.array(miller_index)
        return 1 / ((np.dot(np.dot(hkl, gstar), hkl.T)) ** (1 / 2))

    def __repr__(self):
        outs = [
            "<SDM Lattice>",
            "    abc : {:10.6f} {:10.6f} {:10.6f}".format(*self.parameters[:3]),
            " angles : {:10.6f} {:10.6f} {:10.6f}".format(*self.parameters[3:]),
            " volume : {:10.6f}".format(self.volume),
            "      A : " + " ".join(map(repr, self.matrix[0])),
            "      B : " + " ".join(map(repr, self.matrix[1])),
            "      C : " + " ".join(map(repr, self.matrix[2])),
        ]
        return "\n".join(outs)

    def __eq__(self, other):
        """
        A lattice is considered to be equal to another if the internal matrix
        representation satisfies np.allclose(matrix1, matrix2) to be True.
        """
        if other is None:
            return False
        # shortcut the np.allclose if the memory addresses are the same
        # (very common in Structure.from_sites)
        return self is other or np.allclose(self.matrix, other.matrix)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return "\n".join([" ".join(["%.6f" % i for i in row]) for row in self.matrix])

    def getFastReducedCell(self) -> "Lattice":
        matrix = latt62Matrix(self.parameters)
        rot = kabsch(matrix, self.matrix)
        mul = np.floor(matrix[1, 0] / matrix[0, 0] + 0.5)
        matrix[1] -= mul * matrix[0]
        mul = np.floor(matrix[2, 1] / matrix[1, 1] + 0.5)
        matrix[2] -= mul * matrix[1]
        mul = np.floor(matrix[2, 0] / matrix[0, 0] + 0.5)
        matrix[2] -= mul * matrix[0]
        mat_final = np.dot(matrix, rot)
        return self.__class__(mat_final)

    def getNiggliReducedLattice(self, tol: float = 1e-5) -> "Lattice":
        """
        Get the Niggli reduced lattice using the numerically stable algo
        proposed by R. W. Grosse-Kunstleve, N. K. Sauter, & P. D. Adams,
        Acta Crystallographica Section A Foundations of Crystallography, 2003,
        60(1), 1-6. doi:10.1107/S010876730302186X
        Args:
            tol (float): The numerical tolerance. The default of 1e-5 should
                result in stable behavior for most cases.
        Returns:
            Niggli-reduced lattice.
        """
        matrix = self.getFastReducedCell().matrix
        e = tol * self.volume ** (1 / 3)

        # Define metric tensor
        G = np.dot(matrix, matrix.T)

        # This sets an upper limit on the number of iterations.
        for count in range(100):
            # The steps are labelled as Ax as per the labelling scheme in the
            # paper.
            (A, B, C, E, N, Y) = (
                G[0, 0],
                G[1, 1],
                G[2, 2],
                2 * G[1, 2],
                2 * G[0, 2],
                2 * G[0, 1],
            )

            if A > B + e or (abs(A - B) < e and abs(E) > abs(N) + e):
                # A1
                M = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
                G = dot(transpose(M), dot(G, M))
            if (B > C + e) or (abs(B - C) < e and abs(N) > abs(Y) + e):
                # A2
                M = [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
                G = dot(transpose(M), dot(G, M))
                continue

            l = 0 if abs(E) < e else E / abs(E)
            m = 0 if abs(N) < e else N / abs(N)
            n = 0 if abs(Y) < e else Y / abs(Y)
            if l * m * n == 1:
                # A3
                i = -1 if l == -1 else 1
                j = -1 if m == -1 else 1
                k = -1 if n == -1 else 1
                M = [[i, 0, 0], [0, j, 0], [0, 0, k]]
                G = dot(transpose(M), dot(G, M))
            elif l * m * n == 0 or l * m * n == -1:
                # A4
                i = -1 if l == 1 else 1
                j = -1 if m == 1 else 1
                k = -1 if n == 1 else 1

                if i * j * k == -1:
                    if n == 0:
                        k = -1
                    elif m == 0:
                        j = -1
                    elif l == 0:
                        i = -1
                M = [[i, 0, 0], [0, j, 0], [0, 0, k]]
                G = dot(transpose(M), dot(G, M))

            (A, B, C, E, N, Y) = (
                G[0, 0],
                G[1, 1],
                G[2, 2],
                2 * G[1, 2],
                2 * G[0, 2],
                2 * G[0, 1],
            )

            # A5
            if abs(E) > B + e or (abs(E - B) < e and 2 * N < Y - e) or (abs(E + B) < e and Y < -e):
                M = [[1, 0, 0], [0, 1, -E / abs(E)], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A6
            if abs(N) > A + e or (abs(A - N) < e and 2 * E < Y - e) or (abs(A + N) < e and Y < -e):
                M = [[1, 0, -N / abs(N)], [0, 1, 0], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A7
            if abs(Y) > A + e or (abs(A - Y) < e and 2 * E < N - e) or (abs(A + Y) < e and N < -e):
                M = [[1, -Y / abs(Y), 0], [0, 1, 0], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A8
            if E + N + Y + A + B < -e or (abs(E + N + Y + A + B) < e < Y + (A + N) * 2):
                M = [[1, 0, 1], [0, 1, 1], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            break

        A = G[0, 0]
        B = G[1, 1]
        C = G[2, 2]
        E = 2 * G[1, 2]
        N = 2 * G[0, 2]
        Y = 2 * G[0, 1]
        a = math.sqrt(A)
        b = math.sqrt(B)
        c = math.sqrt(C)
        alpha = math.acos(E / 2 / b / c) / math.pi * 180
        beta = math.acos(N / 2 / a / c) / math.pi * 180
        gamma = math.acos(Y / 2 / a / b) / math.pi * 180

        latt = Lattice.fromParameters(a, b, c, alpha, beta, gamma)

        mapped = self.findMapping(latt, e, skip_rotation_matrix=True)
        if mapped is not None:
            if np.linalg.det(mapped[0].matrix) > 0:
                return mapped[0]
            return Lattice(-mapped[0].matrix)

        raise ValueError("can't find niggli")

    def findAllMappings(
        self,
        other_lattice: "Lattice",
        ltol: float = 1e-5,
        atol: float = 1,
        skip_rotation_matrix: bool = False,
    ) -> Iterator[Tuple["Lattice", Optional[np.ndarray], np.ndarray]]:
        """
        Finds all mappings between current lattice and another lattice.
        Args:
            other_lattice (Lattice): Another lattice that is equivalent to
                this one.
            ltol (float): Tolerance for matching lengths. Defaults to 1e-5.
            atol (float): Tolerance for matching angles. Defaults to 1.
            skip_rotation_matrix (bool): Whether to skip calculation of the
                rotation matrix
        Yields:
            (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
            found. aligned_lattice is a rotated version of other_lattice that
            has the same lattice parameters, but which is aligned in the
            coordinate system of this lattice so that translational points
            match up in 3D. rotation_matrix is the rotation that has to be
            applied to other_lattice to obtain aligned_lattice, i.e.,
            aligned_matrix = np.inner(other_lattice, rotation_matrix) and
            op = SymmOp.from_rotation_and_translation(rotation_matrix)
            aligned_matrix = op.operate_multi(latt.matrix)
            Finally, scale_matrix is the integer matrix that expresses
            aligned_matrix as a linear combination of this
            lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)
            None is returned if no matches are found.
        """
        lengths = other_lattice.abc
        (alpha, beta, gamma) = other_lattice.angles

        frac, dist, _, _ = self.getPointsInSphere_py(
            [[0, 0, 0]], [0, 0, 0], max(lengths) * (1 + ltol), zip_results=False
        )
        cart = self.getCartesianCoords(frac)  # type: ignore
        # this can't be broadcast because they're different lengths
        inds = [np.logical_and(dist / l < 1 + ltol, dist / l > 1 / (1 + ltol)) for l in lengths]  # type: ignore
        c_a, c_b, c_c = (cart[i] for i in inds)
        f_a, f_b, f_c = (frac[i] for i in inds)
        l_a, l_b, l_c = (np.sum(c ** 2, axis=-1) ** 0.5 for c in (c_a, c_b, c_c))

        def get_angles(v1, v2, l1, l2):
            x = np.inner(v1, v2) / l1[:, None] / l2
            x[x > 1] = 1
            x[x < -1] = -1
            angles = np.arccos(x) * 180.0 / np.pi
            return angles

        alphab = np.abs(get_angles(c_b, c_c, l_b, l_c) - alpha) < atol
        betab = np.abs(get_angles(c_a, c_c, l_a, l_c) - beta) < atol
        gammab = np.abs(get_angles(c_a, c_b, l_a, l_b) - gamma) < atol

        for i, all_j in enumerate(gammab):
            inds = np.logical_and(all_j[:, None], np.logical_and(alphab, betab[i][None, :]))
            for j, k in np.argwhere(inds):
                scale_m = np.array((f_a[i], f_b[j], f_c[k]), dtype=np.int_)  # type: ignore
                if abs(np.linalg.det(scale_m)) < 1e-8:
                    continue

                aligned_m = np.array((c_a[i], c_b[j], c_c[k]))

                if skip_rotation_matrix:
                    rotation_m = None
                else:
                    rotation_m = np.linalg.solve(aligned_m, other_lattice.matrix)

                yield Lattice(aligned_m), rotation_m, scale_m

    def findMapping(
        self,
        other_lattice: "Lattice",
        ltol: float = 1e-5,
        atol: float = 1,
        skip_rotation_matrix: bool = False,
    ) -> Optional[Tuple["Lattice", Optional[np.ndarray], np.ndarray]]:
        """
        Finds a mapping between current lattice and another lattice. There
        are an infinite number of choices of basis vectors for two entirely
        equivalent lattices. This method returns a mapping that maps
        other_lattice to this lattice.
        Args:
            other_lattice (Lattice): Another lattice that is equivalent to
                this one.
            ltol (float): Tolerance for matching lengths. Defaults to 1e-5.
            atol (float): Tolerance for matching angles. Defaults to 1.
        Returns:
            (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
            found. aligned_lattice is a rotated version of other_lattice that
            has the same lattice parameters, but which is aligned in the
            coordinate system of this lattice so that translational points
            match up in 3D. rotation_matrix is the rotation that has to be
            applied to other_lattice to obtain aligned_lattice, i.e.,
            aligned_matrix = np.inner(other_lattice, rotation_matrix) and
            op = SymmOp.from_rotation_and_translation(rotation_matrix)
            aligned_matrix = op.operate_multi(latt.matrix)
            Finally, scale_matrix is the integer matrix that expresses
            aligned_matrix as a linear combination of this
            lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)
            None is returned if no matches are found.
        """
        for x in self.findAllMappings(other_lattice, ltol, atol, skip_rotation_matrix=skip_rotation_matrix):
            return x
        return None

    def getPointsInSphere_py(
        self,
        frac_points: ndarray,
        center: ndarray,
        r: float,
        zip_results=True,
    ) -> Union[List[Tuple[np.ndarray, float, int, np.ndarray]], List[np.ndarray]]:
        """
        Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic
        images.
        Algorithm:
        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelpiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.
           Nxmax = r * length_of_b_1 / (2 Pi)
        2. keep points falling within r.
        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                 point, or return the raw fcoord, dist, index arrays
        Returns:
            if zip_results:
                [(fcoord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                fcoords, dists, inds, image
        """
        cart_coords = self.getCartesianCoords(frac_points)
        neighbors = get_points_in_spheres(
            all_coords=cart_coords,
            center_coords=np.array([center]),
            r=r,
            pbc=True,
            numerical_tol=1e-8,
            lattice=self,
            return_fcoords=True,
        )[0]
        if len(neighbors) < 1:
            return [] if zip_results else [()] * 4  # type: ignore
        if zip_results:
            return neighbors
        return [np.array(i) for i in list(zip(*neighbors))]


def get_points_in_spheres(
    all_coords: np.ndarray,
    center_coords: np.ndarray,
    r: float,
    pbc: Union[bool, List[bool]] = True,
    numerical_tol: float = 1e-8,
    lattice: Lattice = None,
    return_fcoords: bool = False,
) -> List[List[Tuple[np.ndarray, float, int, np.ndarray]]]:
    """
    For each point in `center_coords`, get all the neighboring points in `all_coords` that are within the
    cutoff radius `r`.
    Args:
        all_coords: (list of cartesian coordinates) all available points
        center_coords: (list of cartesian coordinates) all centering points
        r: (float) cutoff radius
        pbc: (bool or a list of bool) whether to set periodic boundaries
        numerical_tol: (float) numerical tolerance
        lattice: (Lattice) lattice to consider when PBC is enabled
        return_fcoords: (bool) whether to return fractional coords when pbc is set.
    Returns:
        List[List[Tuple[coords, distance, index, image]]]
    """
    if isinstance(pbc, bool):
        pbc = [pbc] * 3
    pbc = np.array(pbc, dtype=np.bool_)  # type: ignore
    if return_fcoords and lattice is None:
        raise ValueError("Lattice needs to be supplied to compute fractional coordinates")
    center_coords_min = np.min(center_coords, axis=0)
    center_coords_max = np.max(center_coords, axis=0)
    # The lower bound of all considered atom coords
    global_min = center_coords_min - r - numerical_tol
    global_max = center_coords_max + r + numerical_tol
    if np.any(pbc):
        if lattice is None:
            raise ValueError("Lattice needs to be supplied when considering periodic boundary")
        recp_len = np.array(lattice.reciprocal_lattice.abc)
        maxr = np.ceil((r + 0.15) * recp_len / (2 * math.pi))
        frac_coords = lattice.getFractionalCoords(center_coords)
        nmin_temp = np.floor(np.min(frac_coords, axis=0)) - maxr
        nmax_temp = np.ceil(np.max(frac_coords, axis=0)) + maxr
        nmin = np.zeros_like(nmin_temp)
        nmin[pbc] = nmin_temp[pbc]
        nmax = np.ones_like(nmax_temp)
        nmax[pbc] = nmax_temp[pbc]
        all_ranges = [np.arange(x, y, dtype="int64") for x, y in zip(nmin, nmax)]
        matrix = lattice.matrix
        # temporarily hold the fractional coordinates
        image_offsets = lattice.getFractionalCoords(all_coords)
        all_fcoords = []
        # only wrap periodic boundary
        for k in range(3):
            if pbc[k]:  # type: ignore
                all_fcoords.append(np.mod(image_offsets[:, k : k + 1], 1))
            else:
                all_fcoords.append(image_offsets[:, k : k + 1])
        all_fcoords = np.concatenate(all_fcoords, axis=1)
        image_offsets = image_offsets - all_fcoords
        coords_in_cell = np.dot(all_fcoords, matrix)
        # Filter out those beyond max range
        valid_coords = []
        valid_images = []
        valid_indices = []
        for image in itertools.product(*all_ranges):
            coords = np.dot(image, matrix) + coords_in_cell
            valid_index_bool = np.all(
                np.bitwise_and(coords > global_min[None, :], coords < global_max[None, :]),
                axis=1,
            )
            ind = np.arange(len(all_coords))
            if np.any(valid_index_bool):
                valid_coords.append(coords[valid_index_bool])
                valid_images.append(np.tile(image, [np.sum(valid_index_bool), 1]) - image_offsets[valid_index_bool])
                valid_indices.extend([k for k in ind if valid_index_bool[k]])
        if len(valid_coords) < 1:
            return [[]] * len(center_coords)
        valid_coords = np.concatenate(valid_coords, axis=0)
        valid_images = np.concatenate(valid_images, axis=0)

    else:
        valid_coords = all_coords  # type: ignore
        valid_images = [[0, 0, 0]] * len(valid_coords)
        valid_indices = np.arange(len(valid_coords))

    # Divide the valid 3D space into cubes and compute the cube ids
    all_cube_index = _compute_cube_index(valid_coords, global_min, r)  # type: ignore
    nx, ny, nz = _compute_cube_index(global_max, global_min, r) + 1
    all_cube_index = _three_to_one(all_cube_index, ny, nz)
    site_cube_index = _three_to_one(_compute_cube_index(center_coords, global_min, r), ny, nz)
    # create cube index to coordinates, images, and indices map
    cube_to_coords = collections.defaultdict(list)  # type: Dict[int, List]
    cube_to_images = collections.defaultdict(list)  # type: Dict[int, List]
    cube_to_indices = collections.defaultdict(list)  # type: Dict[int, List]
    for i, j, k, l in zip(all_cube_index.ravel(), valid_coords, valid_images, valid_indices):
        cube_to_coords[i].append(j)
        cube_to_images[i].append(k)
        cube_to_indices[i].append(l)

    # find all neighboring cubes for each atom in the lattice cell
    site_neighbors = find_neighbors(site_cube_index, nx, ny, nz)
    neighbors = []  # type: List[List[Tuple[np.ndarray, float, int, np.ndarray]]]

    for i, j in zip(center_coords, site_neighbors):
        l1 = np.array(_three_to_one(j, ny, nz), dtype=int).ravel()
        # use the cube index map to find the all the neighboring
        # coords, images, and indices
        ks = [k for k in l1 if k in cube_to_coords]
        if not ks:
            neighbors.append([])
            continue
        nn_coords = np.concatenate([cube_to_coords[k] for k in ks], axis=0)
        nn_images = itertools.chain(*[cube_to_images[k] for k in ks])
        nn_indices = itertools.chain(*[cube_to_indices[k] for k in ks])
        dist = np.linalg.norm(nn_coords - i[None, :], axis=1)
        nns: List[Tuple[np.ndarray, float, int, np.ndarray]] = []
        for coord, index, image, d in zip(nn_coords, nn_indices, nn_images, dist):
            # filtering out all sites that are beyond the cutoff
            # Here there is no filtering of overlapping sites
            if d < r + numerical_tol:
                if return_fcoords and (lattice is not None):
                    coord = np.round(lattice.getFractionalCoords(coord), 10)
                nn = (coord, float(d), int(index), image)
                nns.append(nn)
        neighbors.append(nns)
    return neighbors


def _compute_cube_index(coords: np.ndarray, global_min: float, radius: float) -> np.ndarray:
    """
    Compute the cube index from coordinates
    Args:
        coords: (nx3 array) atom coordinates
        global_min: (float) lower boundary of coordinates
        radius: (float) cutoff radius
    Returns: (nx3 array) int indices
    """
    return np.array(np.floor((coords - global_min) / radius), dtype=int)


def _three_to_one(label3d: np.ndarray, ny: int, nz: int) -> np.ndarray:
    """
    The reverse of _one_to_three
    """
    return np.array(label3d[:, 0] * ny * nz + label3d[:, 1] * nz + label3d[:, 2]).reshape((-1, 1))


def _one_to_three(label1d: np.ndarray, ny: int, nz: int) -> np.ndarray:
    """
    Convert a 1D index array to 3D index array
    Args:
        label1d: (array) 1D index array
        ny: (int) number of cells in y direction
        nz: (int) number of cells in z direction
    Returns: (nx3) int array of index
    """
    last = np.mod(label1d, nz)
    second = np.mod((label1d - last) / nz, ny)
    first = (label1d - last - second * nz) / (ny * nz)
    return np.concatenate([first, second, last], axis=1)


def find_neighbors(label: np.ndarray, nx: int, ny: int, nz: int) -> List[np.ndarray]:
    """
    Given a cube index, find the neighbor cube indices
    Args:
        label: (array) (n,) or (n x 3) indice array
        nx: (int) number of cells in y direction
        ny: (int) number of cells in y direction
        nz: (int) number of cells in z direction
    Returns: neighbor cell indices
    """

    array = [[-1, 0, 1]] * 3
    neighbor_vectors = np.array(list(itertools.product(*array)), dtype=int)
    if np.shape(label)[1] == 1:
        label3d = _one_to_three(label, ny, nz)
    else:
        label3d = label
    all_labels = label3d[:, None, :] - neighbor_vectors[None, :, :]
    filtered_labels = []
    # filter out out-of-bound labels i.e., label < 0
    for labels in all_labels:
        ind = (labels[:, 0] < nx) * (labels[:, 1] < ny) * (labels[:, 2] < nz) * np.all(labels > -1e-5, axis=1)
        filtered_labels.append(labels[ind])
    return filtered_labels
