# -*- coding: utf-8 -*-
from __future__ import annotations
from collections import deque
import abc
import copy
import json
import numpy as np
from numpy import ndarray

import sdm
from sdm.utils import rotation_matrix
from sdm.models import Structure, Molecule
from .components import Model3D
from . import logger as root_logger
logger = root_logger.getChild('VSession')


class ViewSession:
    _class = "ViewSession"

    _default_style = {
        'Rf_atoms': 1.0,
        'R_bonds': 0.1,
        'R_select': 0.1,
        'R_axes': 0.0,
        'R_cell': 0.0,
        'R_box': None,
        'C_axis_x': (1.0, 0.0, 0.0),
        'C_axis_y': (0.0, 1.0, 0.0),
        'C_axis_z': (0.0, 0.0, 1.0),
        'C_box': (1.0, 1.0, 1.0),
        'C_cell': (1.0, 1.0, 1.0),
        'C_background': (0.0, 0.0, 0.0, 1.0)
    }

    _style_data = {
        "Wireframe": {
            "Rf_atoms": None,
            "R_bonds": 0.0
        },
        "Stick": {
            "Rf_atoms": 0.0,
            "R_bonds": 0.08,
        },
        "Ball-Stick": {
            "Rf_atoms": 1.0,
            "R_bonds": 0.08
        },
        "Spacefill": {
            "Rf_atoms": 3.0,
            "R_bonds": None
        },
        "Dark": {
            "C_background": (0.0, 0.0, 0.0, 1.0),
            "C_cell": (1.0, 1.0, 1.0),
            "C_box": (1.0, 1.0, 1.0)
        },
        "Light": {
            "C_background": (1.0, 1.0, 1.0, 1.0),
            "C_cell": (0.0, 0.0, 0.0),
            "C_box": (0.0, 0.0, 0.0)
        },
    }

    _default_scene = {
        "scale": 1.0,
        "view": (-10, 10, -10, 10, -500, 500),
        "camera": ((0, 0, 1), (0, 0, 0), (0, 1, 0))
    }

    _default_array = {
        "triangles": None,
        "circles": None,
        "lines": None,
        "text": None
    }

    _default_info = {
        "source": "",
    }

    def __init__(self, data=None, info=None, style=None, scene=None):
        # 配置
        self.style = self._default_style.copy()
        self.info = self._default_info.copy()
        self.scene = self._default_scene.copy()
        self.arrays = self._default_array.copy()
        if info:
            self.info.update(info)
        if style:
            self.style.update(style)
        if scene:
            self.scene.update(scene)
        self.scene['view'] = np.array(self.scene['view'], 'f')
        self.scene['camera'] = np.array(self.scene['camera'], 'f')
        # 数据
        self.dataqueue = deque([data], maxlen=20)
        self.dataindex = 0
        self.datareset = copy.copy(data)
        self.selection = []

    # 抽象方法 & 抽象属性
    @abc.abstractproperty
    def coords(self) -> ndarray:
        """select 需要用到坐标"""
        pass

    @abc.abstractproperty
    def radius(self) -> ndarray:
        """select 需要用到半径"""
        pass

    @abc.abstractproperty
    def groups(self) -> list[list]:
        """select 需要用到选择组"""
        pass

    @abc.abstractmethod
    def makeArray(self):
        """制造传递给 VBO 用的 array"""
        pass

    @abc.abstractmethod
    def getSelectInfo(self) -> str:
        pass

    @abc.abstractmethod
    def getDataInfo(self) -> str:
        pass

    # 数据存取 Start
    @property
    def data(self):
        """dataqueue:
        正常: 旧旧旧旧旧当
        撤销: 旧旧当新新新
        新加: 旧旧当加 -> 旧旧旧当
        复位: 当(复)新新新
        """
        return self.dataqueue[self.dataindex]

    @data.setter
    def data(self, val):
        self.dataqueue[self.dataindex] = val

    def dataPush(self, data, preview=False):
        if preview:
            self.data = data
        else:
            for _ in range(self.dataindex+1, len(self.dataqueue)):
                self.dataqueue.pop()
            self.dataqueue.append(data)
            self.dataindex = len(self.dataqueue) - 1
        self.makeArray()

    def dataUndo(self):
        if self.dataindex > 0:
            self.dataindex -= 1
            self.makeArray()

    def dataRedo(self):
        if self.dataindex+1 < len(self.dataqueue):
            self.dataindex += 1
            self.makeArray()

    def dataReset(self):
        self.dataqueue.popleft()
        self.dataqueue.appendleft(copy.copy(self.datareset))
        self.dataindex = 0
        self.makeArray()

    # 原子选择
    def select(self, cur: int, shift=False, dclick=False):
        if cur is None:
            if not self.selection or shift:
                return
            else:
                self.selection = []
        elif dclick:
            curl = [sorted(g) for g in self.groups if cur in g][0]
            if shift:
                if cur in self.selection:
                    self.selection = self.selection + [x for x in curl if x not in self.selection]
                else:
                    self.selection = [x for x in self.selection if x not in curl]
            else:
                self.selection = curl
        elif shift:
            if cur in self.selection:
                self.selection.remove(cur)
            else:
                self.selection.append(cur)
        else:
            self.selection = [cur]
        self.makeSelectArray()

    def makeSelectArray(self):
        if self.style['Rf_atoms']:
            radius = self.radius * self.style['Rf_atoms']
        else:
            radius = self.radius * 0
        self.arrays["circles"] = []
        for idx in self.selection:
            self.addFixedCircle(self.coords[idx], radius[idx] + 0.1)
        if len(self.arrays["circles"]) > 0:
            self.arrays["circles"] = np.vstack(self.arrays["circles"])

    def calcSelectByRay(self, l0: ndarray, l1: ndarray) -> int:
        """鼠标发出一条射线，穿过原子球，选中距离va方向交点最近的原子
                P
            b / |vd
             /  |
            l0------->l1 va
               --> vc
        """
        va = l1 - l0
        vb = self.coords - l0
        vc = (np.dot(va, vb.T) / np.dot(va, va)).reshape((-1, 1)) * np.tile(va, (vb.shape[0], 1))
        vd = vb - vc
        r2d2 = self.radius ** 2 - np.sum(vd * vd, axis=1)
        inrange = (r2d2 >= 0)
        if np.sum(inrange) > 0:
            r2d2i = r2d2[inrange]
            vec = vc[inrange] + np.sqrt(r2d2i).reshape((-1, 1)) / np.linalg.norm(va) * np.tile(va, (r2d2i.shape[0], 1))
            amin = np.argmin(np.dot(va, vec.T))
            return int(np.where(inrange)[0][amin])
        return None

    # 修改样式
    def setStyle(self, name):
        self.style.update(self._style_data.get(name, {}))

    # 模型写入
    def addBall(self, radius: float, center: ndarray, color=[1, 1, 1]):
        vertnorm = Model3D.ball.copy()
        vertnorm[:, :3] *= radius
        vertnorm[:, :3] += np.array(center, 'f')
        vertnorm[:, 6:] = np.array(color, 'f')
        self.arrays['triangles'].append(vertnorm)

    def addStick(self, radius: float, src: ndarray, dst: ndarray, color=[1, 1, 1]):
        vector = dst - src
        length = np.linalg.norm(vector)
        rot_matrix = rotation_matrix(np.cross([0, 0, 1], vector), np.arccos(vector[2] / length))

        vertnorm = Model3D.stick.copy()
        vertnorm[:, :3] *= np.array([[radius, radius, length]], 'f')
        vertnorm[:, :3] = np.dot(vertnorm[:, :3], rot_matrix)
        vertnorm[:, 3:6] = np.dot(vertnorm[:, 3:6], rot_matrix)
        vertnorm[:, :3] += np.array(src, 'f')
        vertnorm[:, 6:] = np.array(color, 'f')
        self.arrays['triangles'].append(vertnorm)

    def addLine(self, src: list, dst: list, color=[1, 1, 1]):
        self.arrays['lines'] += [[src, color], [dst, color]]

    def addFixedCircle(self, center: ndarray, radius: float):
        vertex = np.array([
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, -1, 1],
            [1, 1, 1, -1, -1, 1],
            [1, 1, 1, -1, -1, 1],
            [1, 1, 1, -1, 1, 1],
            [1, 1, 1, 1, 1, 1]
        ], 'f')
        vertex[:, :3] = center
        vertex[:, 3:] *= radius
        self.arrays['circles'].append(vertex)

    def addLattice(self, matrix: ndarray, p_ori=(0, 0, 0)):
        cellorder = [1, 5, 1, 6, 2, 4, 2, 6, 3, 4, 3, 5, 7, 4, 7, 5, 7, 6]
        cellorder2 = [[1, 5], [1, 6], [2, 4], [2, 6], [3, 4], [3, 5], [7, 4], [7, 5], [7, 6]]
        peek = np.sum(matrix, axis=0)
        p = np.vstack([[0, 0, 0], matrix, peek-matrix, peek]).astype('f') + p_ori
        vertex = []
        if self.style['R_axes'] > 0:
            radius = self.style['R_axes']
            self.addStick(radius, p[0], p[1], self.style['C_axis_x'])
            self.addStick(radius, p[0], p[2], self.style['C_axis_y'])
            self.addStick(radius, p[0], p[3], self.style['C_axis_z'])
        else:
            vertex += [
                [p[0], self.style['C_axis_x']], [p[1], self.style['C_axis_x']],
                [p[0], self.style['C_axis_y']], [p[2], self.style['C_axis_y']],
                [p[0], self.style['C_axis_z']], [p[3], self.style['C_axis_z']],
            ]

        if self.style['R_cell'] > 0:
            radius = self.style['R_axes']
            color = self.style['C_cell']
            for i, j in cellorder2:
                self.addStick(radius, p[i], p[j], color)
        else:
            color = self.style['C_cell']
            vertex += [[p[i], color] for i in cellorder]
        self.arrays['lines'] += vertex
    # 模型写入 End


class ViewStruct(ViewSession):
    _class = "ViewStruct"

    @classmethod
    def load_string(cls, filename, string=None, **kwargs):
        """
        kwargs = {info=None, style=None, scene=None}
        """
        if filename == 'json':
            info_dict = json.loads(string)
            if "lattice" in info_dict:
                data = Structure.fromDict(info_dict)
            else:
                data = Molecule.fromDict(info_dict)
        else:
            data = sdm.open(filename, string=string)
        obj = cls(data=data, **kwargs)
        obj.datareset = data.copy()
        return obj

    @property
    def coords(self):
        """select 需要用到坐标"""
        return self.data.coords

    @property
    def radius(self):
        """select 需要用到半径"""
        return np.array([at.vradius for at in self.data.elements])

    @property
    def groups(self):
        """select 需要用到选择组"""
        return self.data.groups

    def getSelectInfo(self) -> list:
        selection = self.selection
        atom_titles = self.data.atom_titles
        info = "%d"%len(selection)
        if self.data is not None:
            if len(selection) == 1:
                info = "[%10.6f, %10.6f, %10.6f]"%tuple(self.coords[selection[0]])
            elif len(selection) == 2:
                info = "%10.6f"%self.data.getDistance(*selection)
                bond = self.data.getBonds([selection])[0]
                bond = "No Bond" if bond is None else bond["type"].name
                info += ", " + bond
            elif len(selection) == 3:
                info = "%10.6f"%self.data.getAngle(*selection)
            elif len(selection) == 4:
                info = "%10.6f"%self.data.getDihedral(*selection)
        infos = [', '.join([str(x) for x in selection]),
                 ', '.join([atom_titles[i] for i in selection]),
                 info]
        return infos

    def getDataInfo(self) -> str:
        return repr(self.data).split('\n', 1)[-1]

    def setCamera(self, along=None, old_camera=None, target=None):
        if along in ("a", "b", "c", "a*", "b*", "c*"):
            if "*" in along:
                mat = self.data.lattice.reciprocal_lattice.matrix
            else:
                mat = self.data.lattice.matrix
            i, j = {"a": (0, 1), "b": (1, 2), "c": (2, 0)}.get(along[0], (0, 1))
            camera = - np.array([mat[i], [0, 0, 0], np.cross(mat[i], mat[j])])
        elif old_camera is not None:
            camera = old_camera
        else:
            camera = np.array(self._default_scene['camera'], dtype='f')

        if target is None:
            target = np.average(self.coords, axis=0)
        camera[0] += target - camera[1]
        camera[1] = target
        self.scene['camera'] = camera
        self.scene['view'][:] = self._default_scene['view']

    def makeArray(self):
        struct = self.data
        elements = struct.elements
        coords = struct.coords
        radius = self.radius
        hexcolor = [at.vcolor for at in elements]
        hexcolor = [[int(c[:2], 16), int(c[2:4], 16), int(c[4:], 16)] for c in hexcolor]
        color = np.array(hexcolor, 'f') / 255
        lattice = getattr(struct, 'lattice', None)

        # 键信息
        bonds = []
        if struct.bonds:
            for i, j, rank in struct.bonds(data='type'):
                ci, cj = coords[[i, j]]
                cm = (ci + cj) / 2
                bonds.append((ci, cm, rank.value, i))
                bonds.append((cj, cm, rank.value, j))
        for i, halfbonds in struct.atoms('halfbonds'):
            if halfbonds:
                src = coords[i]
                for j, vec in halfbonds:
                    bonds.append((i, src, src + vec / 2, 1))

        # 制作 array
        camera = self.scene['camera'][0]
        for k in self.arrays:
            self.arrays[k] = []

        # atoms
        if self.style['Rf_atoms'] is not None:
            if self.style['Rf_atoms'] > 0:
                ra = radius * self.style['Rf_atoms']
            elif self.style['Rf_atoms'] == 0:
                ra = np.ones(len(radius), 'f') * self.style['R_bonds']
            for r, cd, cl in zip(ra, coords, color):
                self.addBall(r, cd, cl)

        # bonds
        if self.style['R_bonds'] is not None and len(bonds) > 0:
            rb = self.style['R_bonds']
            for ci, cm, rank, cl in bonds:
                if rank == 3:
                    drift = np.cross(camera, cm - ci)
                    if rb == 0:
                        drift *= 0.4 / np.linalg.norm(drift)
                        self.addLine(ci + drift, cm + drift, color[cl])
                        self.addLine(ci - drift, cm - drift, color[cl])
                        self.addLine(ci, cm, color[cl])
                    else:
                        drift *= 1.8 * rb / np.linalg.norm(drift)
                        self.addStick(rb, ci + drift, cm + drift, color[cl])
                        self.addStick(rb, ci - drift, cm - drift, color[cl])
                        self.addStick(rb, ci, cm, color[cl])
                elif rank == 2:
                    drift = np.cross(camera, cm - ci)
                    if rb == 0:
                        drift *= 0.2 / np.linalg.norm(drift)
                        self.addLine(ci + drift, cm + drift, color[cl])
                        self.addLine(ci - drift, cm - drift, color[cl])
                    else:
                        drift *= 1.2 * rb / np.linalg.norm(drift)
                        self.addStick(rb, ci + drift, cm + drift, color[cl])
                        self.addStick(rb, ci - drift, cm - drift, color[cl])
                elif rank == 1:
                    if rb == 0:
                        self.addLine(ci, cm, color[cl])
                    else:
                        self.addStick(rb, ci, cm, color[cl])
                else:  # rank == other
                    if rb == 0:
                        pass
                    else:
                        self.addStick(rb * 0.5, ci, cm, color[cl])

        if lattice is not None:
            self.addLattice(lattice.matrix)

        # concentrate
        for k, v in self.arrays.items():
            if len(v) >= 1:
                self.arrays[k] = np.vstack(v).astype('f')
            else:
                self.arrays[k] = None
        self.makeSelectArray()


class ViewCompare(ViewSession):
    _class = "ViewCompare"

    def updateCompare(self, com, source=None) -> None:
        self.data.compare = com
        atoms = com.entity1.atoms + com.entity2.atoms
        coords = np.vstack([com.entity1.coords, com.entity2.coords])
        len1 = len(com.entity1)
        self.data.atom_titles = com.entity1.atom_titles + com.entity2.atom_titles
        self.data.coords = coords
        self.data.radius = np.array([at.size for at in atoms])
        self.data.color = np.array([at.color for at in atoms])
        self.data.color[:len1] = [0, 1, 0]

        self.data.bonds = []
        for i, j, rank in com.entity1.bonds:
            ci, cj = coords[[i, j]]
            cm = (ci + cj) / 2
            self.data.bonds.append((i, j, ci, cm, rank.value))
            self.data.bonds.append((j, i, cj, cm, rank.value))

        for i, j, rank in com.entity2.bonds:
            i += len1
            j += len1
            ci, cj = coords[[i, j]]
            cm = (ci + cj) / 2
            self.data.bonds.append((i, j, ci, cm, rank.value))
            self.data.bonds.append((j, i, cj, cm, rank.value))

        if source:
            self.data.source = source
            self.origin_data = copy.deepcopy(self.data)

    def getSelectStatus(self):
        if self.data.compare is not None:
            if len(self.selection) == 1:
                data = "[%10.6f, %10.6f, %10.6f]"%tuple(self.data.coords[self.selection[0]])
            elif len(self.selection) == 2:
                data = "%10.6f"%self.data.compare.getDistance(*self.selection)
            elif len(self.selection) == 3:
                data = "%10.6f"%self.data.compare.getAngle(*self.selection)
            elif len(self.selection) == 4:
                data = "%10.6f"%self.data.compare.getDihedral(*self.selection)
        info = [', '.join([str(x) for x in self.selection]),
                ', '.join([self.data.atom_titles[i] for i in self.selection]),
                data]
        return info
