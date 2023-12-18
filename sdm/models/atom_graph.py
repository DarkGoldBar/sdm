# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/13
from __future__ import annotations
from typing import Optional, Sequence
from enum import IntEnum
from networkx import Graph
from numpy import concatenate
import numpy as np
import networkx as nx
from networkx.classes.reportviews import NodeView, EdgeView
from .periodic_table import Element, PERIODIC_INDEX
from ..analysis.graph import isRotatableBond as _isRotatableBond
from ..utils import atom_idx_to_title


class BondType(IntEnum):
    UNDEFINED = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMA = 4


class BaseGraph():
    """拓扑的基础功能类

    Args:
        graph: nx.Graph

    Properties:
        atoms -> NodeView
        bonds -> EdgeView
        atom_index -> list
        atom_titles -> list[str]
        elements -> list
        natom -> int
        mass -> float
        cycles -> list[list[int]]
        groups -> list[list]
    """
    def __init__(self, graph: Optional[Graph] = None) -> None:
        self.graph = Graph(graph)

    def __len__(self) -> int:
        return len(self.graph)

    @property
    def atoms(self) -> NodeView:
        """
        >>> st.atoms
        NodeView((0, 1, 2, ...))

        >>> st.atoms[0]
        {'element': < N @ 139810084190064 >, 'title': 'N_00'}

        >>> st.atoms('title')
        NodeDataView({0: 'N_00', 1: 'C_01', 2: 'N_02', ...})
        """
        return self.graph.nodes

    @property
    def bonds(self) -> EdgeView:
        """
        >>> self.bonds
        EdgeView([(0, 1), (0, 2), (0, 3), ...])

        >>> self.bonds[0, 1]
        {'type': <BondType.SINGLE: 1>}

        >>> self.bonds(1)
        EdgeDataView([(1, 0), (1, 4), (1, 7)])

        >>> self.bonds([1, 3])
        EdgeDataView([(1, 0), (1, 4), (1, 7), (3, 0)])

        >>> self.bonds(1, data='type')
        EdgeDataView([(1, 0, <BondType.SINGLE: 1>), (1, 4, <BondType.SINGLE: 1>), (1, 7, <BondType.SINGLE: 1>)])
        """
        return self.graph.edges

    @property
    def atom_index(self) -> list:
        """Return 从小到大的原子序"""
        return sorted(self.graph.nodes)

    @property
    def atom_titles(self) -> list[str]:
        """Return 原子名称列表"""
        return self.getAtomAttribute('title', '')

    @property
    def atom_charges(self) -> list[int]:
        """Return 原子名称列表"""
        return self.getAtomAttribute('charge', 0)

    @property
    def elements(self) -> list:
        """Return 原子的元素列表"""
        return self.getAtomAttribute('element')

    @property
    def natom(self) -> int:
        """Return 原子数"""
        return len(self.graph)

    @property
    def mass(self) -> float:
        """Return 原子总质量"""
        return sum(e.atomic_mass for e in self.elements)

    @property
    def cycles(self) -> list[list[int]]:
        """Return 包含的环"""
        return list(nx.algorithms.cycles.cycle_basis(self.graph))

    @property
    def groups(self) -> list[list]:
        """Return 按键连分组"""
        return [sorted(g) for g in self.getConnectedGroups()]

    def _addNodes(self, input_atoms, **kwargs) -> list:
        """增加原子
        input_atoms: int | str | Element
        """
        start = 0 if not self.natom else (max(self.graph.nodes) + 1)
        new_nodes = []

        for idx, iat in enumerate(input_atoms):
            attrib = {"element": None, "title": None}
            if isinstance(iat, str):
                attrib['element'] = Element(iat)
            elif isinstance(iat, int):
                attrib['element'] = Element(PERIODIC_INDEX[iat])
            else:
                attrib['element'] = iat
            new_nodes.append((idx + start, attrib))
        self.graph.add_nodes_from(new_nodes)

        for k, v in kwargs.items():
            self.setAtomAttribute(k, v)
        self._setNodeTitles()
        return [i[0] for i in new_nodes]

    def _delNodes(self, index_list: Sequence[int]) -> None:
        """删除原子与对应连接"""
        self.graph.remove_nodes_from(index_list)
        mapping = {v: k for k, v in enumerate(sorted(self.graph.nodes))}
        nx.relabel_nodes(self.graph, mapping, copy=False)

    def _addGraph(self, graph: Graph) -> None:
        """合并另一个图"""
        start = 0 if not self.natom else (max(self.graph.nodes) + 1)
        new_map = {i: i + start for i in graph.nodes}
        new_graph = nx.relabel_nodes(graph, new_map)
        self.graph = nx.algorithms.operators.union(self.graph, new_graph)

    def _getNodes(self, mask):
        graph = self.graph.subgraph(list(mask))
        mapping = {v: k for k, v in enumerate(sorted(graph.nodes))}
        nx.relabel_nodes(graph, mapping, copy=False)
        return self.__class__(graph)

    def _setNodeTitles(self):
        values = {k: atom_idx_to_title(self.atoms[k]['element'], k) for k, v in self.atoms('title') if v is None}
        self.setAtomAttribute('title', values)

    def getConnectedGroups(self, exclude_atoms=[], exclude_bonds=[]) -> list:
        graph = nx.Graph(self.graph)
        graph.remove_nodes_from(exclude_atoms)
        graph.remove_edges_from(exclude_bonds)
        return list(nx.connected_components(graph))

    def getNeighbor(self, atoms: Sequence[int]) -> set:
        return set(j for i in atoms for j in self.graph.adj[i]) - set(atoms)

    def setAtomAttribute(self, attrib, values) -> None:
        if isinstance(values, dict):
            for idx, val in values.items():
                self.graph.nodes[idx][attrib] = val
        elif values:
            for idx, val in zip(self.atom_index, values):
                self.graph.nodes[idx][attrib] = val

    def getAtomAttribute(self, attrib, default=None) -> list:
        return [self.graph.nodes[idx].get(attrib, default) for idx in self.atom_index]

    def delAtomAttribute(self, attrib) -> None:
        for i in self.graph.nodes:
            self.graph.nodes[i].pop(attrib, None)

    def getBonds(self, bonds_i_j: Sequence[list[int]]) -> list[dict]:
        return [self.graph.get_edge_data(i, j) for i, j in bonds_i_j]

    def setBonds(self, bonds_i_j_t: Sequence[list[int]]) -> None:
        self.graph.add_edges_from([(i, j, {"type": BondType(t)}) for i, j, t in bonds_i_j_t])

    def delBonds(self, bonds_i_j: Sequence[list[int]]) -> None:
        self.graph.remove_edges_from(bonds_i_j)

    def findAtomByTitle(self, atom_title: str) -> int:
        r = [k for k, v in self.atoms(data='title') if v == atom_title]
        return r[0] if r else None

    def findRotatableBonds(self, return4=False) -> list[list[int]]:
        """柔性键的判断标准：
        1. 键: 单键 & 非环
        2. 两端原子: 非顶点 & 非全连H & 非全连F
        3. 不是特殊的 CN 或 CC 三键形式 (线性结构)
        4. 非特殊结构-COOH, -NH2

        Args:
            return4: True -> 返回4原子表示的二面角; False -> 返回2原子表示的键
        """
        cycles = self.cycles
        atom_symbols = [atom.symbol for atom in self.elements]
        rotatable = []
        for i, j, t in self.bonds(data='type'):
            if t.value != 1:
                continue
            if len([c for c in cycles if i in c and j in c]) > 0:
                continue
            if not _isRotatableBond(self.graph.adj, i, j, atom_symbols):
                continue
            if not _isRotatableBond(self.graph.adj, j, i, atom_symbols):
                continue
            rotatable.append((i, j))

        if return4:
            dihedrals = []
            for i, j in rotatable:
                ni = min(set(self.graph.adj[i]) - {j})
                nj = min(set(self.graph.adj[j]) - {i})
                dihedrals.append((ni, i, j, nj))
            return dihedrals
        else:
            return rotatable

    def isRotatableBond(self, idx1: int, idx2: int) -> bool:
        atom_symbols = [atom.symbol for atom in self.elements]
        if (idx1, idx2) not in self.bonds:
            return False
        if self.graph.edges[(idx1, idx2)]['type'].value != 1:
            return False
        if any(c for c in self.cycles if idx1 in c and idx2 in c):
            return False
        if not _isRotatableBond(self.graph.adj, idx1, idx2, atom_symbols):
            return False
        if not _isRotatableBond(self.graph.adj, idx2, idx1, atom_symbols):
            return False
        return True

    def isChiralAtom(self, idx: int) -> bool:
        from ..analysis.graph import BFSTraversal
        adj = self.graph.adj
        ele = [e.name for e in self.elements]
        if self.atoms[idx]['element'].symbol in ('C', 'Si'):
            if adj[idx] != 4:
                return False
        elif self.atoms[idx]['element'].symbol in ('N', 'P'):
            if adj[idx] != 3:
                return False
        else:
            return False
        bfst = BFSTraversal(adj, idx, data=ele)
        for i in range(1, len(bfst)):
            for j in range(i):
                if bfst[i] == bfst[j]:
                    return True
        return False

    def GUIgetFollowerGroups(self, atoms: list[int]) -> list[list[int]]:
        """GUI 变换前获取每个原子的随动原子
        """
        conn_gs = self.getConnectedGroups(exclude_atoms=atoms)
        follow_gs = {k: [] for k in atoms}
        for g in conn_gs:
            g_nei = self.getNeighbor(g) & set(atoms)
            if len(g_nei) == 1:
                k = g_nei.pop()
                follow_gs[k] += list(g)
        return [follow_gs[k] for k in atoms]


class AtomGraph(BaseGraph):
    """原子图对象
    包含向量属性管理, 每个向量属性应该和图的长度相等
    """
    VECTOR_PROPS = []

    def __init__(self, graph=None, **kwargs):
        self.graph = Graph(graph)
        for vkey in self.VECTOR_PROPS:
            vval = kwargs.get(vkey)
            vval = np.empty((0, 3)) if vval is None else vval
            if len(vval) != len(self.graph):
                msg = '向量属性 {} 与原子长度不符: {} != {}'.format(vkey, len(vval), len(self.graph))
                raise ValueError(msg)
            setattr(self, vkey, vval)

    def addAtoms(self, input_atoms, **kwargs):
        """添加新原子, 必须包含所有向量参数"""
        natom = len(input_atoms)
        vp_new = {}
        for vkey in self.VECTOR_PROPS:
            vval = kwargs.pop(vkey)
            if vval is None:
                msg = '缺少向量属性 {}'.format(vkey)
                raise ValueError(msg)
            if len(vval) != natom:
                msg = '向量属性 {} 与原子长度不符: {} != {}'.format(vkey, len(vval), natom)
                raise ValueError(msg)
            vp_new[vkey] = concatenate([getattr(self, vkey), vval])
        self._addNodes(input_atoms, **kwargs)
        for k, v in vp_new.items():
            setattr(self, k, v)

    def addAtomGraph(self, other):
        """合并两个原子图, 把传入参数放在后面"""
        vp_new = {}
        for vkey in self.VECTOR_PROPS:
            vp_new[vkey] = concatenate([getattr(self, vkey), getattr(other, vkey)])
        self._addGraph(other.graph)
        for k, v in vp_new.items():
            setattr(self, k, v)

    def delAtoms(self, indexes):
        """从原子图中删除节点, 并删除节点对应的向量值"""
        mask = np.ones(len(self.graph), dtype=bool)
        mask[indexes] = False
        vp_new = {vkey: getattr(self, vkey)[mask] for vkey in self.VECTOR_PROPS}
        self._delNodes(indexes)
        for k, v in vp_new.items():
            setattr(self, k, v)

    def getAtoms(self, mask, reorder=True):
        """从原子图中截取一部分, 可用于原子重排序"""
        obj = self.copy()
        G = self.graph.subgraph(list(mask))
        # G = nx.Graph()
        # node_view = self.graph.nodes(data=True)
        # nodes = [(i, node_view[i]) for i in mask]
        # G.add_nodes_from(nodes)
        # G.add_edges_from(self.graph.edges(data=True))

        if reorder:
            mapping = {m: i for i, m in enumerate(mask)}
            G = nx.relabel_nodes(G, mapping)
            for vkey in self.VECTOR_PROPS:
                vval = getattr(self, vkey)[mask]
                setattr(obj, vkey, vval)
        obj.graph = G
        return obj

    def copy(self, **kwargs):
        vp_new = {vkey: getattr(self, vkey).copy() for vkey in self.VECTOR_PROPS}
        return self.__class__(self.graph.copy(), **vp_new, **kwargs)
