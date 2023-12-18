# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/13
import numpy as np
import networkx as nx
from .molecule import Molecule
from .pbc import PBC, Lattice
from .symmetry import SpaceGroup, SymmOp
from ..analysis.find_bonds import find_empiric_bonds_pbc
from ..utils import frac_distance2


class Structure(Molecule, PBC):
    """晶体结构

    构建对象：
    Structure.new(
        atom_symbols=['H', 'C', 'C', 'H'],
        coords=[[0, 0, 1], [0, 0, 0.6], [0, 0, -0.6], [0, 0, -1]],
        matrix=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
        int_number=1)
    """
    def __init__(
            self, 
            graph: nx.Graph,
            coords: np.ndarray,
            lattice: Lattice,
            space_group: SpaceGroup,
            title='Untitled',
            symprec=0.001
        ):
        super(Structure, self).__init__(graph, coords, title=title)
        PBC.__init__(self, lattice, symprec=symprec)
        self.space_group = space_group
        self.symprec = symprec

    def __repr__(self):
        s_mol = Molecule.__repr__(self).replace('Molecule', 'Structure', 1).splitlines()
        s_lat = repr(self.lattice).splitlines()
        s_spg = repr(self.space_group).splitlines()
        s_mol[-1] += ' × ' + s_spg[1]
        s = s_mol + s_lat[1:3] + s_spg[2:4]
        return '\n'.join(s)

    def copy(self):
        return super(Molecule, self).copy(
            title=self.title,
            lattice=self.lattice.copy(),
            space_group=self.space_group.copy(),
            symprec=self.symprec)

    def findBonds(self, method='empiric', detect_cut=1.25):
        """找键"""
        alias = {
            'empiric-pbc': 'ep',
        }
        if method is None:
            self.delBonds(self.graph.edges)
            self.delAtomAttribute('halfbonds')
        elif not method:
            return
        else:
            method = alias.get(method, method)
        if method == 'ep':
            radius = np.array([at.covelent_radius for at in self.elements])
            affine_matrices = [op.affine_matrix for op in self.space_group.symmetry_ops]
            halfbonds = find_empiric_bonds_pbc(self.frac_coords, radius, self.lattice.parameters, affine_matrices)
            self.setAtomAttribute("halfbonds", halfbonds)
        else:
            super(Structure, self).findBonds(method=method, detect_cut=detect_cut)

    def getP1Structure(self, nodup=True):
        """转换为 P1 对称性"""
        struct = self.copy()
        struct.title = self.title
        struct.space_group = SpaceGroup.fromHallNumber(1)
        new_frac_coords = [self.frac_coords]
        op_xyz = SymmOp(np.eye(4))
        for op in self.space_group.symmetry_ops:
            if op == op_xyz:
                continue
            struct._addGraph(self.graph)
            new_frac_coords.append(op.operate_multi(self.frac_coords))
        struct.frac_coords = np.vstack(new_frac_coords)
        if nodup:
            remove_index = struct.delDuplicatedAtoms()
            if len(remove_index) > 0:
                print('[getP1Structure] atoms removed:', len(remove_index))
        return struct

    def getAsymmetricStructure(self, angle_tol=-1.0, p1=True):
        """
        调用Spglib, 找出当前结构的非对称单元

        Args:
            angle_tol: (float) 晶胞角度误差
            p1: (bool) 是否转为p1, 在不想应用现有对称操作时传入Flase
        """
        struct = self.getP1Structure() if p1 else self
        space_group = SpaceGroup.fromSpglib(struct.toSpglibCell(p1=False), symprec=self.symprec, angle_tol=angle_tol)
        eqatom = space_group.data_with_struct['equivalent_atoms']
        remove_index = [eqatom[x] in eqatom[:x] for x in range(eqatom.shape[0])]
        remove_index = np.where(remove_index)[0]
        struct.delAtoms(remove_index)
        struct.space_group = space_group
        # struct.title += '_asymmetric'
        return struct

    def getPrimitiveStructure(self, angle_tol=-1.0):
        """调用Spglib, 找出当前结构的原胞"""
        struct = self.getP1Structure()
        space_group = SpaceGroup.fromSpglib(struct.toSpglibCell(p1=False), symprec=self.symprec, angle_tol=angle_tol)
        prime_map = space_group.data_with_struct['mapping_to_primitive'].tolist()
        lattice = Lattice(space_group.data_with_struct['primitive_lattice'])
        is_primitive = [prime_map.index(i) for i in sorted(set(prime_map))]
        remove_index = [i for i in struct.atom_index if i not in is_primitive]
        struct.delAtoms(remove_index)
        struct.space_group = space_group
        struct.lattice = lattice
        # struct.title += '_primitive'
        return struct

    def getSupercell(self, index, index_low=[0, 0, 0], group=None):
        """扩超胞, 不含展开对称性

        Args:
            index: (int[3]) packing上标
            index_low: (int[3]) packing下标
            group: 仅展开选中的原子
        """
        struct = self.copy()
        st_add = self.copy()
        if group is not None:
            st_add.delAtoms(sorted(set(range(len(st_add))) - set(group)))
        for i in range(index_low[0], index[0]):
            for j in range(index_low[1], index[1]):
                for k in range(index_low[2], index[2]):
                    if (i, j, k) != (0, 0, 0):
                        struct._addGraph(st_add.graph)
                        struct._addCoords(st_add.coords + np.dot([i, j, k], st_add.lattice.matrix))
        return struct

    def findEquivalentAtomAndSymmetry(self):
        """找对称等价原子，和对应的对称操作"""
        stp1 = self.getP1Structure()
        spg = SpaceGroup.fromSpglib(stp1.toSpglibCell())
        equal_atom = spg.data_with_struct['equivalent_atoms']
        equal_atom_sym = - np.ones(equal_atom.shape[0], int)
        fcoord = stp1.frac_coords
        fcoord -= np.floor(fcoord)
        for idx, op in enumerate(self.space_group.symmetry_ops):
            oped = op.operate_multi(fcoord[equal_atom])
            oped -= np.floor(oped)
            frac_vec = np.abs(oped - fcoord)
            frac_vec = 0.5 - np.abs(frac_vec - 0.5)
            cart_vec = frac_vec @ self.lattice.matrix
            is_close = (cart_vec**2).sum(axis=1) < self.symprec**2
            equal_atom_sym[is_close] = idx

        assert np.all(equal_atom_sym >= 0), 'some atoms can`t find symmetry'
        return equal_atom, equal_atom_sym

    def fixSymmetryByReference(self, reference):
        """根据参照结构修复对称性，要求当前结构为P1，且原子顺序与参照结构展开P1相同"""
        equal_atom, equal_atom_sym = reference.findEquivalentAtomAndSymmetry()
        invops = [op.inverse for op in reference.space_group.symmetry_ops]
        fcoords_as = {}
        for i, (ipos, isym) in enumerate(zip(equal_atom, equal_atom_sym)):
            crd = invops[isym].operate(self.frac_coords[i])
            if ipos in fcoords_as:
                fcoords_as[ipos].append(crd)
            else:
                fcoords_as[ipos] = [crd]

        fcoords = []
        for i in range(len(reference)):
            if i in fcoords_as:
                crd = np.array(fcoords_as[i])
                crd -= np.floor(crd)
                fdiff = crd - crd[0]
                fdiff[fdiff > 0.5] -= 1
                fdiff[fdiff < -0.5] += 1
                crd = crd[0] + np.average(fdiff, axis=0)
                crd -= np.floor(crd)
            else:
                crd = reference.space_group.symmetry_ops[equal_atom_sym[i]].operate(fcoords[equal_atom[i]])
            fcoords.append(crd)

        st = reference.copy()
        st.lattice = self.lattice
        st.frac_coords = fcoords
        st.moveAtomsByMolecule()
        return st

    def delDuplicatedAtoms(self):
        """删除同一个坐标重复出现的原子（保留原子序靠前的）"""
        frac = self.frac_coords
        mat = self.lattice.matrix
        is_close = frac_distance2(frac, frac, mat) < self.symprec**2
        remove_bool = np.any(np.tril(is_close, -1), axis=1)
        remove_index = np.where(remove_bool)[0]
        self.delAtoms(remove_index)
        return remove_index

    def transformSpaceGroup(self, new_space_group: SpaceGroup):
        """仅支持相同int_number的空间群互相转换"""
        if new_space_group.int_number != self.space_group.int_number:
            raise ValueError('Unsupported transformation')
        mat, vec = self.space_group.getTransformationMatrixAndVector(new_space_group.universal_symbol)
        if mat is None:
            raise ValueError('Unsupported transformation')
        self.lattice = Lattice(np.dot(mat, self.lattice.matrix))
        self.frac_coords -= vec
        self.space_group = new_space_group

    def transformCentering(self, old='P', new='P'):
        """转换晶胞中心类型, 更新晶胞、对称操作

        old: 原晶胞中心类型
        new: 目标晶胞中心类型
        """
        self.lattice = self.lattice.transformCentering(old, new)
        self.space_group = self.space_group.transformCentering(old, new)

    def _getVaspStyleMap(self):
        symbol_type = sorted(set([at.symbol for at in self.elements]))
        mask = []
        idx = 0
        for symbol in symbol_type:
            for i, ele in self.atoms('element'):
                if ele.symbol == symbol:
                    mask.append(i)
                    idx += 1
        return mask

    @classmethod
    def new(cls, frac_coords=None, coords=None, lattice=None, latt6=(1, 1, 1, 90, 90, 90),
            atom_symbols=None, atom_numbers=None, atom_titles=None, atom_charges=None,
            affine_matrices=None, symmops_xyz=None, spg_symbol=None, hall_number=None, int_number=1,
            bonds=None, bond_method=None, title='Untitled', symprec=1e-4, **kwargs):
        if lattice is not None:
            lattice = Lattice(lattice)
        else:
            lattice = Lattice.fromParameters(*latt6)

        if frac_coords is not None:
            coords = lattice.getCartesianCoords(frac_coords)
        elif coords is None:
            raise ValueError('No coords info')

        if affine_matrices is not None:
            space_group = SpaceGroup([SymmOp(aff) for aff in affine_matrices])
        elif symmops_xyz is not None:
            space_group = SpaceGroup.fromXyzOps(symmops_xyz)
        elif spg_symbol is not None:
            space_group = SpaceGroup.fromSymbol(spg_symbol)
        elif hall_number is not None:
            space_group = SpaceGroup.fromHallNumber(hall_number)
        else:
            space_group = SpaceGroup.fromIntNumber(int_number)

        if atom_symbols is not None:
            atoms = atom_symbols
        elif atom_numbers is not None:
            atoms = [int(x) for x in atom_numbers]

        obj = cls(None, None, lattice=lattice, space_group=space_group, title=title, symprec=symprec)
        obj.addAtoms(atoms, coords=coords, title=atom_titles, charge=atom_charges)
        if bonds is not None:
            obj.setBonds(bonds)
        elif bond_method:
            obj.findBonds(bond_method)
        return obj

    @classmethod
    def fromMolecule(cls, molecule: Molecule, lattice=None, latt6=(25, 25, 25, 90, 90, 90),
                     affine_matrices=None, symmops_xyz=None, spg_symbol=None, hall_number=None, int_number=1,
                     symprec=1e-4):
        if lattice is not None:
            lattice = Lattice(lattice)
        else:
            lattice = Lattice.fromParameters(*latt6)

        if affine_matrices is not None:
            space_group = SpaceGroup([SymmOp(aff) for aff in affine_matrices])
        elif symmops_xyz is not None:
            space_group = SpaceGroup.fromXyzOps(symmops_xyz)
        elif spg_symbol is not None:
            space_group = SpaceGroup.fromSymbol(spg_symbol)
        elif hall_number is not None:
            space_group = SpaceGroup.fromHallNumber(hall_number)
        else:
            space_group = SpaceGroup.fromIntNumber(int_number)
        obj = cls(molecule.graph, molecule.coords, lattice=lattice, space_group=space_group, title=molecule.title, symprec=symprec)
        return obj

    @classmethod
    def fromDict(cls, info_dict, bond_method='empiric', symprec=0.001):
        ext_dict = {k: info_dict.pop(k) for k in info_dict if k.startswith('_')}

        obj = cls.new(**info_dict, bond_method=bond_method, symprec=symprec)

        # 特殊处理
        if ext_dict.get('_lattice_centering'):
            center = ext_dict['_lattice_centering']
            if center != 'P':
                obj.transformCentering('P', center)
        elif ext_dict.get('_symmops_centring'):
            center = ext_dict['_symmops_centring']
            if center != 'P':
                obj.space_group = obj.space_group.transformCentering('P', center)
        return obj

    def toDict(self, less=True, vasp_style=False, crystal_style=False):
        vasp_map = None
        crystal_cartops = None
        struct = self

        if vasp_style:
            struct = self.getP1Structure()
            vasp_map = struct._getVaspStyleMap()
            struct = struct.getAtoms(vasp_map)
        if crystal_style:
            crystal_cartops = self.space_group.getCartesianOps(struct.lattice.matrix)

        info_dict = {
            "title": struct.title,
            "lattice": struct.lattice.matrix,
            "symmops_xyz_string": struct.space_group.xyz_string,
            "coords": struct.coords,
            "atom_symbols": [at.symbol for at in struct.elements],
            "atom_titles": struct.atom_titles,
            "bonds": [(i, j, t.value) for i, j, t in struct.bonds(data='type')],
        }
        if less:
            return info_dict

        info_dict_ext = {
            "latt6": struct.lattice.parameters,
            "spg_int_number": struct.space_group.int_number,
            "spg_int_symbol": struct.space_group.int_symbol_full,
            "affine_matrices": [op.affine_matrix for op in struct.space_group.symmetry_ops],
            "lattice_system": struct.lattice.getLatticeSystem(struct.space_group.int_number),
            "frac_coords": struct.frac_coords,
            "atom_numbers": [at.Z for at in struct.elements],
            "VASP_map": vasp_map,
            "CRYSTAL_cartops": crystal_cartops,
        }
        info_dict.update(info_dict_ext)
        return info_dict

    def toString(self, fmt, vasp_style=False, crystal_style=False, **kwargs):
        from sdm.io import getWriter
        writter = getWriter(fmt)
        info_dict = self.toDict(less=False, vasp_style=vasp_style, crystal_style=crystal_style)
        return writter(info_dict, **kwargs)

    @classmethod
    def fromString(cls, string, fmt, **kwargs):
        from sdm.io import getParser
        parser = getParser(fmt)
        info_dict = parser(string, **kwargs)
        return cls.fromDict(info_dict, **kwargs)

    @classmethod
    def fromFile(cls, filename, **kwargs):
        with open(filename) as f:
            return cls.fromString(f.read(), filename, **kwargs)

    def toSpglibCell(self, p1=True):
        """转换为 Spglib Cell 数据结构"""
        st_p1 = self.getP1Structure() if p1 else self
        return (
            st_p1.lattice.matrix.tolist(),
            st_p1.frac_coords.tolist(),
            [a.Z for a in st_p1.elements]
        )

    def toMgStruct(self, sorted_primitive_reduced=False):
        """转换为 Pymatgen Structure 对象"""
        from pymatgen.core import Structure as PmgStructure
        st_p1 = self.getP1Structure() if self.space_group.int_number != 1 else self
        st_mg = PmgStructure(
            self.lattice.matrix,
            [atom.symbol for atom in st_p1.elements],
            st_p1.frac_coords,
        )
        if sorted_primitive_reduced:
            st_mg = st_mg.get_sorted_structure()
            st_mg = st_mg.get_primitive_structure()
            st_mg = st_mg.get_reduced_structure()
        return st_mg
