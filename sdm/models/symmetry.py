# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/13
from __future__ import annotations
from typing import Sequence
import re
import copy
import spglib
import hashlib
import itertools
import numpy as np
import warnings
from fractions import Fraction
from numpy import ndarray
from ..data import HM_INFO, NUM_INFO, SIGN_TO_HM, HALL_TO_HM, HM_SWITCH


class SpaceGroup():
    def __init__(self, symmetry_ops: Sequence[SymmOp], symprec=1e-5) -> None:
        self.symprec = symprec
        self.symmetry_ops = symmetry_ops
        self._int_symbol = None
        self._data = None
        self.data_with_struct = {}

    def _get_sign(self):
        """空间群签名"""
        signs = sorted(op.sign for op in self.symmetry_ops)
        md5str = hashlib.md5(b"".join(signs)).hexdigest()
        return md5str

    def __eq__(self, other):
        if isinstance(other, SpaceGroup):
            return self._get_sign() == other._get_sign()
        return NotImplemented

    def __repr__(self):
        info = {
            "id": hex(id(self)),
            "ops": len(self.symmetry_ops),
            "int_number": self.int_number,
            "int_symbol": self.int_symbol_full,
            "hall_number": self.hall_number,
            "hall_symbol": self.hall_symbol,
        }
        s = [
            "<SDM Spacegroup @ {id}>",
            "{ops} symmetry operation",
            "    ITA : [{int_number:3d}] {int_symbol}",
            "   hall : [{hall_number:3d}] {hall_symbol}"
        ]
        return "\n".join(s).format(**info)

    @property
    def int_symbol(self):
        """Hermann-Mauguin 表示"""
        if self._int_symbol is None:
            sign = self._get_sign()
            symbol = SIGN_TO_HM.get(sign, None)

            if symbol is None:
                rotations = []
                translations = []
                for op in self.symmetry_ops:
                    rotations.append(op.rotation_matrix)
                    trans = (op.affine_matrix[:3, 3] * 12 % 12).round() / 12
                    translations.append(trans)
                number = spglib.get_hall_number_from_symmetry(
                    rotations,
                    translations,
                    symprec=self.symprec,
                )
                if number in HALL_TO_HM:
                    warnings.warn(f'数据库中没有匹配的操作: {sign}')
                    symbol = HALL_TO_HM[number]
                else:
                    warnings.warn(f'找不到对应的空间群信息: {sign}')
                    symbol = ''
            self._int_symbol = symbol
        return self._int_symbol

    @property
    def data(self):
        """根据 hall number 从 spglib 数据库获取数据"""
        if self._data is None:
            self._data = copy.deepcopy(HM_INFO.get(self.int_symbol, {}))
        return self._data

    @property
    def int_number(self):
        return self.data.get('sg_num', 0)

    @property
    def hall_number(self):
        return self.data.get("hall_num", 0)

    @property
    def int_symbol_full(self):
        return self.data.get('full_h_m', 'none')

    @property
    def universal_symbol(self):
        return self.data.get('universal_h_m', 'none')

    @property
    def hall_symbol(self):
        return self.data.get('hall', 'none')

    @property
    def lattice_system(self) -> str:
        return self.data.get('lattice_system', 'none')

    @property
    def xyz_string(self) -> list:
        return [symmop.as_xyz_string() for symmop in self.symmetry_ops]

    def copy(self):
        return SpaceGroup([op.copy() for op in self.symmetry_ops], symprec=self.symprec)

    def transformCentering(self, old='P', new='P') -> SpaceGroup:
        """转换晶胞中心类型, 更新对称操作

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
        symmop_translate_vector = {
            "A": [(0, 0, 0), (0.0, 0.5, 0.5)],
            "B": [(0, 0, 0), (0.5, 0.0, 0.5)],
            "C": [(0, 0, 0), (0.5, 0.5, 0.0)],
            "F": [(0, 0, 0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
            "I": [(0, 0, 0), (0.5, 0.5, 0.5)],
            "H": [(0, 0, 0), (2 / 3., 1 / 3., 1 / 3.), (1 / 3., 2 / 3., 2 / 3.)],
        }
        old = old.upper()
        new = new.upper()
        new_symmops = []
        if old in "PR" and new in "ABCFIH":
            # M-1 * aff * M + translate
            M = np.eye(4)
            Mrot = np.array(transform_matrix.get(new)).reshape((3, 3))
            M[:3, :3] = Mrot
            M_1 = np.linalg.inv(M)
            tvecs = symmop_translate_vector[new]
            for op in self.symmetry_ops:
                for vec in tvecs:
                    affine = op.affine_matrix.copy()
                    affine = M_1.dot(affine).dot(M)
                    affine[:3, 3] += vec
                    affine[:3, 3] -= np.floor(affine[:3, 3])
                    new_symmops.append(SymmOp(affine))

        elif old in "ABCFIH" and new in "PR":
            #  M * aff * M-1 & 去重
            M = np.eye(4)
            Mrot = np.array(transform_matrix.get(old)).reshape((3, 3))
            M[:3, :3] = Mrot
            M_1 = np.linalg.inv(M)
            for op in self.symmetry_ops:
                affine = op.affine_matrix.copy()
                affine = M.dot(affine).dot(M_1)
                affine[:3, 3] -= np.floor(affine[:3, 3])
                new_op = SymmOp(affine)
                if any(new_op == _op for _op in new_symmops):
                    continue
                new_symmops.append(new_op)
        else:
            raise ValueError("Not implemented!")
        return SpaceGroup(new_symmops)

    def getCartesianOps(self, matrix: ndarray) -> list[SymmOp]:
        matinv = np.linalg.inv(matrix)
        cart_ops = []
        for op in self.symmetry_ops:
            affine = np.eye(4)
            affine[:3, :3] = np.dot(np.dot(matinv, op.affine_matrix[:3, :3]), matrix).round(10)
            affine[:3, 3] = np.dot(op.affine_matrix[:3, 3], matrix).round(10)
            cart_ops.append(SymmOp(affine))
        return cart_ops

    @classmethod
    def fromCartesianOps(cls, matrix: ndarray, affine_cart: ndarray, **kwargs) -> SpaceGroup:
        matinv = np.linalg.inv(matrix)
        symmetry_ops = []
        for mat in affine_cart:
            affine = np.eye(4)
            affine[:3, :3] = np.dot(np.dot(matrix, mat[:3, :3]), matinv).round(10)
            affine[:3, 3] = np.dot(mat[:3, 3], matinv).round(10)
            symmetry_ops.append(SymmOp(affine))
        return cls(symmetry_ops, **kwargs)

    @classmethod
    def fromXyzOps(cls, ops: Sequence[str], **kwargs) -> SpaceGroup:
        """从操作字符串构造"""
        return cls([SymmOp.from_xyz_string(op) for op in ops], **kwargs)

    @classmethod
    def fromIntSymbol(cls, symbol, **kwargs) -> SpaceGroup:
        """从准确的空间群名构造"""
        data = HM_INFO[symbol]
        obj = cls.fromXyzOps(data['symops'])
        obj._int_symbol = symbol
        return obj

    @classmethod
    def fromHallNumber(cls, hall_number: int, **kwargs) -> SpaceGroup:
        """从Hall群号构造空间群"""
        symbol = HALL_TO_HM[hall_number]
        return cls.fromIntSymbol(symbol)

    @classmethod
    def fromIntNumber(cls, int_number: int, **kwargs) -> SpaceGroup:
        """从H-M群号构造空间群"""
        symbol = NUM_INFO[int_number]['default_hm']
        return cls.fromIntSymbol(symbol)

    @classmethod
    def fromSymbol(cls, symbol: str) -> SpaceGroup:
        """从空间群名称构造

        检索键: ['hermann_mauguin', 'full_h_m', 'universal_h_m', 'hall', 'schoenflies']
        """
        aka = {
            "C2": "C121",
            "C2/c": "C12/c1",
            "C2/m": "C12/m1",
            "Cc": "C1c1",
            "Cm": "C1m1",
            "P2": "P121",
            "P2/c": "P12/c1",
            "P2/m": "P12/m1",
            "P2_1": "P12_11",
            "P2_1/c": "P12_1/c1",
            "P2_1/m": "P12_1/m1",
            "Pc": "P1c1",
            "Pm": "P1m1"
        }
        symbol = symbol[0].upper() + symbol[1:].lower()
        short = symbol.strip().replace(" ", "")
        if short in aka:
            hm_symbol = aka[short]
        else:
            for k, v in HM_INFO.items():
                if (symbol in [v['full_h_m'], v['hall']]) or (short in [v['hermann_mauguin'], v['universal_h_m'], v['schoenflies']]):
                    hm_symbol = k
                    break
            else:
                raise ValueError("没有找到空间群符号 {}".fomrat(symbol))
        return cls.fromIntSymbol(hm_symbol)

    @classmethod
    def fromSpglib(cls, spglib_cell: list, angle_tol=-1.0, symprec=1e-5) -> SpaceGroup:
        """
        Interface to spglib

        Find symmetry operations from a crystal structure

        spglib_cell : Structure.toSpglibCell()

        symprec : float
            Symmetry search tolerance in the unit of length.
        angle_tolerance : float
            Symmetry search tolerance in the unit of angle deg. If the value is
            negative, an internally optimized routine is used to judge symmetry.
        """
        symmetry_dataset = spglib.get_symmetry_dataset(spglib_cell, symprec=symprec, angle_tolerance=angle_tol)
        obj = cls([SymmOp.from_rotation_and_translation(r, t) for r, t in zip(symmetry_dataset['rotations'], symmetry_dataset['translations'])], symprec=symprec)
        obj.data_with_struct = symmetry_dataset
        return obj

    def getTransformationMatrixAndVector(self, universal_symbol: str) -> tuple:
        '''找出同一空间群下不同对称操作组合的变换矩阵和平移向量'''
        key = (self.universal_symbol, universal_symbol)
        op = SymmOp.from_xyz_string(HM_SWITCH[key])
        return op.rotation_matrix, op.translation_vector

    def sort(self) -> None:
        '''把原始操作放在第一个'''
        op0 = SymmOp(np.eye(4))
        self.symmetry_ops.sort(key=lambda op: op == op0, reverse=True)

class SymmOp():
    def __init__(self, affine_transformation_matrix: ndarray):
        """
        对称操作对象, 接受 4x4 对称操作矩阵.

        使用 from_rotation_and_translation 构造函数, 从旋转矩阵和平移向量构造

        Args:
            affine_transformation_matrix (4x4 array)
        """
        affine_transformation_matrix = np.array(affine_transformation_matrix)
        if affine_transformation_matrix.shape != (4, 4):
            raise ValueError("Affine Matrix must be a 4x4 numpy array!")
        self.affine_matrix = affine_transformation_matrix
        self.affine_matrix.setflags(write=False)
        self._sign = None

    def copy(self):
        return SymmOp(self.affine_matrix.copy())

    @classmethod
    def from_rotation_and_translation(cls, rot, trans):
        """
        从旋转矩阵和平移向量构造对称操作对象

        Args:
            rotation_matrix (3x3 array): Rotation matrix.
            translation_vec (3x1 array): Translation vector.
        Returns:
            SymmOp object
        """
        rotation_matrix = np.array(rot)
        translation_vec = np.array(trans)
        if rotation_matrix.shape != (3, 3):
            raise ValueError("Rotation Matrix must be a 3x3 numpy array.")
        if translation_vec.shape != (3,):
            raise ValueError("Translation vector must be a rank 1 numpy array " "with 3 elements.")
        affine_matrix = np.eye(4)
        affine_matrix[0:3][:, 0:3] = rotation_matrix
        affine_matrix[0:3][:, 3] = translation_vec
        return cls(affine_matrix)

    def __eq__(self, other):
        return self.sign == other.sign

    def __mul__(self, other):
        """乘法 合成两个对称操作"""
        new_matrix = np.dot(self.affine_matrix, other.affine_matrix)
        return self.__class__(new_matrix)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        output = [
            "Rot:",
            str(self.affine_matrix[0:3][:, 0:3]),
            "tau",
            str(self.affine_matrix[0:3][:, 3]),
        ]
        return "\n".join(output)

    @property
    def rotation_matrix(self) -> np.ndarray:
        """3x3 旋转矩阵"""
        return self.affine_matrix[0:3][:, 0:3]

    @property
    def translation_vector(self) -> np.ndarray:
        """1x3 平移向量"""
        return self.affine_matrix[0:3][:, 3]

    @property
    def inverse(self) -> "SymmOp":
        """逆操作"""
        invr = np.linalg.inv(self.affine_matrix)
        return SymmOp(invr)

    @property
    def sign(self) -> tuple:
        """操作签名, 精度为 1/12, 用于比较和哈希"""
        if self._sign is None:
            self._sign = self.rotation_matrix.round().astype('int8').tobytes()
            self._sign += ((self.translation_vector * 12).round() % 12).astype('int8').tobytes()
        return self._sign

    def operate(self, point):
        """
        Apply the operation on a point.
        Args:
            point: Cartesian coordinate.
        Returns:
            Coordinates of point after operation.
        """
        affine_point = np.array([point[0], point[1], point[2], 1])
        return np.dot(self.affine_matrix, affine_point)[0:3]

    def operate_multi(self, points):
        """
        Apply the operation on a list of points.
        Args:
            points: List of Cartesian coordinates
        Returns:
            Numpy array of coordinates after operation
        """
        points = np.array(points)
        affine_points = np.concatenate([points, np.ones(points.shape[:-1] + (1,))], axis=-1)
        return np.inner(affine_points, self.affine_matrix)[..., :-1]

    def as_xyz_string(self, components=("x", "y", "z"), delim=", ") -> str:
        """
        Args:
            components: either ('x', 'y', 'z') or ('a', 'b', 'c')
            delim: delimiter

        Returns a string of the form 'x, y, z', '-x, -y, z', '-y+1/2, x+1/2, z+1/2'
        """
        parts = []
        for rot, t in zip(self.rotation_matrix, self.translation_vector):
            s = []
            for r, dim in zip(rot, components):
                if r != 0:
                    f = self._fraction(r)
                    f = f.replace("/", dim + "/") if "/" in f else f + dim  # 加符号
                    f = re.sub(r"^(-?)1" + dim, r"\1" + dim, f)  # 去1
                    s.append(f)
            s.append(self._fraction(t))
            s = "+".join(s) if s else "0"
            s = s.replace("+-", "-").replace("+0", "")
            parts.append(s)
        return delim.join(parts)

    @classmethod
    def from_xyz_string(cls, xyz_string: str):
        """
        Args:
            xyz_string: string of the form 'x, y, z', '-x, -y, z', '-y+1/2, x+1/2, z+1/2'

        Returns:
            SymmOp
        """
        rot_matrix = np.zeros((3, 3))
        trans = np.zeros(3)
        toks = xyz_string.strip().replace(" ", "").lower().split(",")
        re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
        re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
        for i, tok in enumerate(toks):
            # build the rotation matrix
            for m in re_rot.finditer(tok):
                factor = -1.0 if m.group(1) == "-" else 1.0
                if m.group(2) != "":
                    factor *= float(m.group(2)) / float(m.group(3)) if m.group(3) != "" else float(m.group(2))
                j = ord(m.group(4)) - 120
                rot_matrix[i, j] = factor
            # build the translation vector
            for m in re_trans.finditer(tok):
                factor = -1 if m.group(1) == "-" else 1
                num = float(m.group(2)) / float(m.group(3)) if m.group(3) != "" else float(m.group(2))
                trans[i] = num * factor
        return cls.from_rotation_and_translation(rot_matrix, trans)

    @staticmethod
    def _fraction(number):
        frac = Fraction(number)
        f12 = str(frac.limit_denominator(12))
        f1e4 = str(frac.limit_denominator(10000))
        if f12 != f1e4:
            msg = '对称操作小数精度不足, 建议修复对称性 {} -> {}'.format(f12, f1e4)
            warnings.warn(msg)
        return f12
