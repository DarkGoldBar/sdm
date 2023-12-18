# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/13
from __future__ import annotations
from typing import Optional, Sequence
from numpy import ndarray
import numpy as np

from sdm.analysis.find_bonds import detectDistance
from ..utils import rotate
_type_SeqInt = Optional[Sequence[int]]


class Space():
    """空间坐标的基础功能类

    Properties:
        coords: np.ndarray (shape=n*3)
    """
    def __init__(self, coords: ndarray = None) -> None:
        if coords is None:
            coords = []
        self.coords = np.array(coords).reshape((-1, 3))

    def _addCoords(self, coords: ndarray) -> None:
        """增加坐标"""
        self.coords = np.vstack([self.coords, coords])

    def _delCoords(self, index_list: Sequence[int]) -> None:
        """删除坐标"""
        mask = np.ones(self.coords.shape[0], bool)
        mask[index_list] = False
        self.coords = self.coords[mask]

    def _reorderCoords(self, old2new: dict) -> None:
        """坐标重排序"""
        mask = [k for k, v in sorted(old2new.items(), key=lambda kv: kv[1])]
        self.coords = self.coords[mask]

    def getDistance(self, atom1: int, atom2: int) -> float:
        """返回两个原子之间的距离"""
        return np.linalg.norm(self.coords[atom1] - self.coords[atom2])

    def getAngle(self, atom1: int, atom2: int, atom3: int) -> float:
        """返回三个原子之间的键角"""
        v1 = self.coords[atom1] - self.coords[atom2]
        v2 = self.coords[atom3] - self.coords[atom2]
        d = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
        d = max(min(d, 1), -1)
        angle = np.arccos(d)
        return angle / np.pi * 180

    def getDihedral(self, atom1: int, atom2: int, atom3: int, atom4: int) -> float:
        """返回四个原子之间的二面角"""
        v1 = self.coords[atom2] - self.coords[atom1]
        v2 = self.coords[atom3] - self.coords[atom2]
        v3 = self.coords[atom4] - self.coords[atom3]
        v12 = np.cross(v1, v2)
        v23 = np.cross(v2, v3)
        dihe = np.arctan2(np.dot(v12, v3) * np.linalg.norm(v2), np.dot(v12, v23))
        return dihe / np.pi * 180

    def getCentroid(self, group: _type_SeqInt = None) -> ndarray:
        """返回一组原子的质心"""
        group = slice(None) if group is None else list(group)
        coords = self.coords[group]
        weight = np.array([e.atomic_mass for e in self.elements])[group]
        return np.average(coords, axis=0, weights=weight)

    def getGeometryCenter(self, group: _type_SeqInt = None) -> ndarray:
        """返回一组原子的几何中心"""
        group = slice(None) if group is None else list(group)
        return np.average(self.coords[group], axis=0)

    def getBorderBox(self, group: _type_SeqInt = None, expand=0) -> ndarray:
        """返回一组原子的边界框"""
        group = slice(None) if group is None else list(group)
        pos_min = np.min(self.coords[group], axis=0) - expand
        pos_max = np.max(self.coords[group], axis=0) + expand
        mask = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ])
        bbox = pos_max * mask + pos_min * (1-mask)
        return bbox

    def translate(self, vector: ndarray, group: _type_SeqInt = None) -> None:
        """平移"""
        group = slice(None) if group is None else list(group)
        self.coords[group] += vector

    def rotateWithAxisAngle(self, axis_p0: ndarray, axis_p1: ndarray,
                            theta: float, group: _type_SeqInt = None) -> None:
        """轴角旋转"""
        group = slice(None) if group is None else list(group)
        self.coords[group] = rotate(self.coords[group], axis_p0, axis_p1, theta)

    def rotateWithMatrix(self, matrix, center: Optional[ndarray] = None, group: _type_SeqInt = None) -> None:
        """矩阵旋转"""
        group = slice(None) if group is None else list(group)
        if center is None:
            center = np.zeros(3)
        self.coords[group] = np.dot(self.coords[group] - center, matrix) + center

    def GUItranslate(
            self,
            atoms: list[int],
            target: float,
            groups: list[list[int]],
            config: list[str]) -> None:
        """GUI 的平移接口

        atoms: 两个原子的index,
        target: 目标长度,
        groups: 随动原子组,
        config: 每个原子的配置 [[TG|TA|FIX], [TG|TA|FIX]]

        config 选项:
            TG: Translate Group
            TA: Translate Atom
            FIX: Fix Atom
        """
        assert not all(x == "FIX" for x in config), "bad config"
        dist = self.getDistance(*atoms)
        vector = self.coords[atoms[1]] - self.coords[atoms[0]]
        vector *= (target - dist) / dist
        if config[0] == "FIX":
            vecs = [None, vector]
        elif config[1] == "FIX":
            vecs = [-vector, None]
        else:
            vecs = [-vector / 2, vector / 2]

        for i in range(2):
            if vecs[i] is not None:
                group = [atoms[i]] + groups[i] if config[i] == "TG" else [atoms[i]]
                self.translate(vecs[i], group)

    def GUIrotate(
            self,
            atoms: list[int],
            target: float,
            groups: list[list[int]],
            config: list[str]) -> None:
        """GUI 的平移接口

        atoms: 三个原子的index,
        target: 目标长度,
        groups: 随动原子组,
        config: 每个原子的配置 [[RG|TG|RA|FIX], [TG|TA|FIX], [RG|TG|RA|FIX]]

        config 选项:
            TG: Translate Group
            TA: Translate Atom
            RG: Rotate Group
            RA: Rotate Atom
            FIX: Fix Atom
        """
        assert not (config[0] == config[2] == 'FIX'), "bad config"
        theta = (target - self.getAngle(*atoms)) / 180 * np.pi
        vec10 = self.coords[atoms[1]] - self.coords[atoms[0]]
        vec21 = self.coords[atoms[2]] - self.coords[atoms[1]]
        v_axis = np.cross(vec10, vec21)
        axis_p0 = self.coords[atoms[1]]
        axis_p1 = axis_p0 + v_axis
        old_coords = self.coords[atoms]
        groups_r1 = [atoms[0]] + groups[0] if config[0] == "RG" else [atoms[0]]
        groups_r3 = [atoms[2]] + groups[2] if config[2] == "RG" else [atoms[2]]

        if config[1] == "FIX" and config[0] == "FIX":
            self.rotateWithAxisAngle(axis_p0, axis_p1, theta, groups_r3)
        elif config[1] == "FIX" and config[2] == "FIX":
            self.rotateWithAxisAngle(axis_p0, axis_p1, -theta, groups_r1)
        elif config[0] == "FIX" and config[2] == "FIX":
            raise NotImplementedError("FIX, xx, FIX not implemented")
        else:
            self.rotateWithAxisAngle(axis_p0, axis_p1, theta / 2, groups_r1)
            self.rotateWithAxisAngle(axis_p0, axis_p1, -theta / 2, groups_r3)
            if config[0] == "FIX":
                vector = old_coords[0] - self.coords[atoms[0]]
            elif config[1] == "FIX":
                vector = np.zeros(3)
            elif config[2] == "FIX":
                vector = old_coords[2] - self.coords[atoms[2]]
            else:
                vector = np.average(old_coords - self.coords[atoms], axis=0)
            group = groups_r1 + [atoms[1]] + groups_r3
            self.translate(vector, group)
        for i in range(3):
            if config[i] == "TG":
                self.translate(self.coords[atoms[i]] - old_coords[i], groups[i])

    def GUItwist(
            self,
            atoms: list[int],
            target: float,
            groups: list[list[int]],
            config: list[str]) -> None:
        """GUI 的扭转接口

        value: 目标角度
        groups: 随动原子
        config: [[RS|RG|RA|FIX], [RS|RG|RA|FIX]]

        config 选项:
            RS: Rotate Groups (including groups at atom2,atom3)
            RG: Rotate Group
            RA: Rotate Atom
            FIX: Fix Atom
        """
        assert not all(x == "fix" for x in config), "bad config"
        dihedral = self.getDihedral(*atoms)
        theta = ((target - dihedral) / 180 * np.pi)
        # print(f"{theta=}, {target=}, {dihedral=}")
        if config[0] == "RA":
            group1 = [atoms[0]]
        elif config[0] == "RG":
            group1 = groups[0] + [atoms[0]]
        elif config[0] == "RS":
            group1 = groups[0] + groups[1] + [atoms[0]]

        if config[1] == "RA":
            group2 = [atoms[3]]
        elif config[1] == "RG":
            group2 = groups[3] + [atoms[3]]
        elif config[1] == "RS":
            group2 = groups[3] + groups[2] + [atoms[3]]

        axis_p0, axis_p1 = self.coords[atoms[1:3]]
        if config[0] == "FIX":
            self.rotateWithAxisAngle(axis_p0, axis_p1, theta, group2)
        elif config[1] == "FIX":
            self.rotateWithAxisAngle(axis_p0, axis_p1, -theta, group1)
        else:
            self.rotateWithAxisAngle(axis_p0, axis_p1, -theta / 2, group1)
            self.rotateWithAxisAngle(axis_p0, axis_p1, theta / 2, group2)
