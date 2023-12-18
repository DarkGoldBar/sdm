# -*- coding: utf-8 -*-
# Author : bochen.li
# Date : 2022/06/15
from itertools import combinations, permutations, product
import networkx as nx
from networkx import Graph


"""针对RMSD计算前, 寻找map路径的优化

当前问题: 使用DFS暴搜会在特定情况下效率非常低

1.1 给原子分组, 分组进行排列组合
1.2 大小组嵌套时, 只能乘法遍历
# TODO: 1.3 组内选择时, 按照连接性生成组合 (苯, 立方烷的情况)

# TODO: 2. 鞭炮状分子, 虽然每一对对称原子都有两种选项, 但是这些选择之前完全没有关联,
#    应该按顺序用贪心法解决, 而不是搜索 (需要个别处理, 需要坐标)
"""


class PathFinder:
    """Molecule compare path finder

    All atom_attr and bond_attr must be hashable

    Node Attributes:
        b: bfs-label
        d: hashable data
        e: equivalent-atoms

    Example:
    >>> p = PathFinder(mol, True, {'element': None}, {'type': None})
    >>> p = PathFinder(mol, True, {'element': None})
    """
    def __init__(self, graph,
                 without_H=True,
                 atom_attr: dict = None,
                 bond_attr: dict = None) -> None:
        self.atom_attr = atom_attr if atom_attr else {}
        self.bond_attr = bond_attr if bond_attr else {}
        self.without_H = without_H
        self._graph = None
        self._equal_atoms = []
        self._equal_groups = []
        self.groups = []

        self.graph = graph

    @property
    def graph(self) -> Graph:
        return self._graph

    @graph.setter
    def graph(self, graph: Graph):
        self._graph = Graph()
        nodes = []
        for i, d in graph.nodes.data():
            if self.without_H and d['element'].number == 1:
                continue
            data = [d.get(k, self.atom_attr[k]) for k in d if k in self.atom_attr]
            data.append(len(graph[i]))
            nodes.append((i, {'d': tuple(data)}))
        self._graph.add_nodes_from(nodes)
        edges = []
        for i, j, d in graph.edges.data():
            if i not in self._graph or j not in self._graph:
                continue
            data = [d.get(k, self.bond_attr[k]) for k in d if k in self.bond_attr]
            edges.append((i, j, {'d': tuple(data)}))
        self._graph.add_edges_from(edges)

    def _set_bfs_label(self, node):
        G = self._graph
        assert node in G
        old = set([node])
        cur = set([node])
        nxt = set(G.adj[node])
        label = [[(G.nodes[node]['d'], )], ]
        while nxt:
            nxt_label = []
            for i in nxt:
                nxt_label_i = sorted([G.edges[(i, j)]['d'] for j in cur if G.has_edge(i, j)])
                nxt_label_i.insert(0, G.nodes[i]['d'])
                nxt_label.append(tuple(nxt_label_i))
            nxt_label = sorted(nxt_label, key=str)
            label.append(nxt_label)
            old = old | nxt
            cur = nxt
            nxt = set(i for j in cur for i in G.adj[j]) - old
        G.nodes[node]['bfs'] = label
        G.nodes[node]['b'] = str(label)

    def set_bfs_label(self):
        for i in self._graph.nodes:
            self._set_bfs_label(i)

    def set_groups(self):
        G = self._graph
        H = Graph()
        for i, j in combinations(G.nodes, 2):
            if G.nodes[i]['b'] == G.nodes[j]['b']:
                H.add_edge(i, j)
        I = Graph(G.subgraph(H.nodes))
        I.remove_edges_from(H.edges)
        self._equal_atoms = list(nx.connected_components(H))
        self._equal_groups = list(nx.connected_components(I))

        self.groups = []
        for cc in nx.connected_components(H):
            f = any(j in nx.node_connected_component(I, i) for i, j in combinations(cc, 2))
            if not f:
                group = [nx.node_connected_component(I, i) for i in cc]
                for i, j in combinations(cc, 2):
                    I.add_edge(i, j)
                self.groups.append(group)

        for g in self._equal_groups:
            if len(g) > 2:
                sub_G = G.subgraph(list(g))
                p = PathFinder(sub_G, False)
                p._graph = sub_G
                p.set_groups()
                self.groups += p.groups

    @property
    def path_count(self):
        count = 1
        for g in self.groups:
            count *= len(g)
        return count

    def path_iter(self):
        groups = [[sorted(x) for x in g] for g in self.groups]
        chained = [[i for g in gg for i in g] for gg in groups]
        _mp = {i: i for i in self.graph.nodes}
        pm_list = [permutations(range(len(x)), len(x)) for x in groups]
        prod = product(*pm_list)
        for choice in prod:
            mp = _mp.copy()
            for cho, gp, src in zip(choice, groups, chained):
                tgt = [i for ci in cho for i in gp[ci]]
                upd = {i: mp[j] for i, j in zip(src, tgt)}
                mp.update(upd)
            yield mp
