# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/13
import json
import numpy as np
from collections.abc import Iterable
from .atom_graph import AtomGraph
from .space import Space
from ..analysis.find_bonds import find_babel_bonds, find_empiric_bonds
from ..analysis.internal_coord import find_internal_coord_path, find_internal_coord_order, find_initial_4_chain, calc_internal_coords, calc_cart_from_internal


class Molecule(AtomGraph, Space):
    """分子对象, 包含 拓扑 + 空间 的功能, 以及拓扑空间结合的功能

    mol[0]       -> infomation of atom
    mol[[1,2,3]] -> part of molecule
    mol1 + mol2  -> combine of molecule
    """
    VECTOR_PROPS = ['coords']

    def __init__(self, graph=None, coords=None, title='Untitled'):
        super(Molecule, self).__init__(graph=graph, coords=coords)
        self.title = title

    def __repr__(self):
        info = {'id': hex(id(self)), 'title': self.title, 'natom': len(self)}
        s = ['<SDM Molecule at {id}>', '{title}', '{natom} atoms']
        return '\n'.join(s).format(**info)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.graph.nodes[key], self.graph[key], self.coords[key]
        if isinstance(key, slice):
            mask = self.atom_index[key]
        elif isinstance(key, Iterable):
            mask = sorted(key)
        else:
            raise TypeError('key must be Iterable')
        obj = self.getAtoms(mask)
        return Molecule(obj.graph, obj.coords)

    def __add__(self, other) -> "Molecule":
        mol = self.copy()
        mol.addAtomGraph(other)
        mol.title += ' + ' + other.title
        return mol

    def copy(self):
        return super(Molecule, self).copy(title=self.title)

    def findBonds(self, method='empiric', *args, **kwargs):
        """找键
            'empiric': 'e',  # 完全按照距离判断
            'idatm': 'i',  # 距离判断 + 杂化判断
            'babel': 'b',  # 调用 openbabel
            'obabel': 'b3',  # 调用 openbabel3
            'babel3': 'b3',
        """
        alias = {
            'empiric': 'e',
            'idatm': 'i',
            'babel': 'b',
            'obabel': 'b3',
            'babel3': 'b3',
        }
        if method is None:
            self.graph.remove_edges_from(self.graph.edges)
        if not method:
            return
        else:
            method = alias.get(method, method)

        if method == 'e':
            radius = np.array([at.covelent_radius for at in self.elements])
            bonds = find_empiric_bonds(self.coords, radius, *args, **kwargs)
            self.setBonds(bonds)
        elif method == 'i':
            from sdm.analysis.find_bonds_idatm import find_idatm_bonds
            d = self.toDict()
            bonds = find_idatm_bonds(d['atom_symbols'], d['coords'])
            self.setBonds(bonds)
        elif method == 'b':
            bonds = find_babel_bonds(self.toString("xyz", remove_titles=True))
            self.setBonds(bonds)
        elif method == 'b3':
            raise NotImplementedError()
        else:
            raise ValueError("No such bond method {}".format(method))

    def rotateAloneRotatableBond(self, atom1, atom2, theta=0., use_rad=False):
        """
        idx1: 原子编号，固定的
        idx2: 原子编号，旋转的
        theta: 旋转角度
        use_rad: True: 弧度; False: 角度
        """
        assert self.isRotatableBond(atom1, atom2)
        if not use_rad:
            theta = np.deg2rad(theta)
        groups = self.getConnectedGroups(exclude_bonds=[(atom1, atom2)])
        group = [g for g in groups if atom2 in g]
        self.rotateWithAxisAngle(self.coords[atom1], self.coords[atom2], theta, group=group)

    def getInternalCoords(self) -> tuple:
        # 起始点避开柔性键的端点
        _dihedrals = self.findRotatableBonds(return4=False)
        exclude = list(set(i for ii in _dihedrals for i in ii))
        init4 = find_initial_4_chain(self.graph, exclude=exclude)
        # 找到内坐标的原子排序
        order = find_internal_coord_order(self.graph, init4=init4)
        mol = self.getAtoms(order)
        # 找出内坐标路径, 有柔性角的按照柔性角顺序
        paths = find_internal_coord_path(mol.graph)
        dihedrals = mol.findRotatableBonds(return4=True)
        for d in dihedrals:
            if d[3] > d[0]:
                d = d[::-1]
            paths[d[0]] = list(d)
        # 计算内坐标数值
        internal_paths = [paths[i] for i in mol.atom_index]
        internal_coords = calc_internal_coords(mol, internal_paths)
        return mol, internal_paths, internal_coords

    @classmethod
    def new(cls, coords=None, zmatrix=None, zpath=None,
            atom_symbols=None, atom_numbers=None, atom_titles=None, atom_charges=None,
            bonds=None, bond_method=None, title='Untitled', **kwargs):
        if zmatrix is not None:
            coords = calc_cart_from_internal(zpath, zmatrix)
        if atom_symbols is not None:
            atoms = atom_symbols
        elif atom_numbers is not None:
            atoms = [int(x) for x in atom_numbers]

        obj = cls(title=title)
        obj.addAtoms(atoms, coords=coords, titles=atom_titles, charges=atom_charges)
        if bonds is not None:
            obj.setBonds(bonds)
        elif bond_method:
            obj.findBonds(bond_method)
        return obj

    @classmethod
    def fromDict(cls, info_dict, bond_method=None):
        return cls.new(**info_dict, bond_method=bond_method)

    def toDict(self, less=False):
        ele = self.elements
        return {
            'title': self.title,
            'coords': self.coords,
            'atom_titles': self.atom_titles,
            'atom_symbols': [at.symbol for at in ele],
            'atom_numbers': [at.Z for at in ele],
            'atom_charges': self.atom_charges,
            'bonds': [(i, j, t.value) for i, j, t in self.bonds(data='type')],
        }

    def toJson(self):
        info_dict = self.toDict(less=True)
        return json.dumps({k: v.tolist() if hasattr(v, "tolist") else v for k, v in info_dict.items()})

    def toString(self, fmt, **kwargs):
        from sdm.io import getWriter
        writter = getWriter(fmt)
        info_dict = self.toDict()
        return writter(info_dict, **kwargs)

    def save(self, filename, **kwargs):
        string = self.toString(filename, **kwargs)
        with open(filename, 'w') as f:
            f.write(string)

    def show(self, host=None, port=None):
        """把结构发送到窗口"""
        import os
        from sdm_gui.server import MyClient
        data = json.dumps({
            'string': self.toJson(),
            'origin': f"PID:{os.getpid()}",
            'type': 'json'
        })
        MyClient(host, port).sendall(data)
