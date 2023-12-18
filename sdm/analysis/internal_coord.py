# -*- coding:utf-8 -*-
from numpy import cross, array, linalg, ndarray, deg2rad, sin, cos, dot, eye
from networkx import Graph
from networkx.algorithms import single_source_shortest_path


def find_internal_coord_path(G: Graph):
    """ 找到每个原子的 m-b-a-d 路径
    要求保证图G中, 后面的节点一定与前面的节点相连
    """
    paths = {}
    nodes = sorted(G.nodes)
    for i, node in enumerate(nodes):
        path = single_source_shortest_path(G, node, cutoff=3)
        length = min(i+1, 4)
        for j in nodes[:i+1]:
            if j in path and len(path[j]) == length:
                paths[node] = path[j]
                break
        else:
            raise ValueError("无法找到合适的路径 {}".format(node))
    return paths


def find_internal_coord_order(G: Graph, init4=None):
    """找到合适的内坐标分子的原子顺序, 使每一个原子都连接在它前面的原子上

    要求: init4 给出一串的4个原子
    """
    if init4 is None:
        init4 = find_initial_4_chain(G)
    qout = list(init4)
    qin = sorted((i for i in G.nodes if i not in qout), reverse=True)
    wait_set = set()
    while qin or wait_set:
        connected = set(j for i in qout for j in G.adj[i]) - set(qout)
        work = wait_set & connected
        i = None
        if work:
            i = min(work)
            wait_set.remove(i)
        elif qin:
            i = qin.pop()
            if i not in connected:
                wait_set.add(i)
                continue
        else:
            raise RuntimeError("未找到链接的原子 {}".format(wait_set))
        qout.append(i)
    return qout


def find_initial_4_chain(G: Graph, exclude=[]):
    """找到序号最低的一串的4个原子

    Args:
        exclude: 需要排除的原子, 比如在柔性键两端的原子
    """
    for i in sorted(G.nodes):
        if i in exclude:
            continue
        path = single_source_shortest_path(G, i, cutoff=3)
        path4 = sorted([p for p in path.values() if len(p) == 4])
        if path4:
            return path4[0]


def calc_internal_coords(mol, paths: list) -> list:
    """从分子对象和内坐标路径计算内坐标"""
    icoords = []
    for path in paths:
        lp = len(path)
        ic = []
        if lp > 1:
            ic.append(mol.getDistance(*path[:2]))
        if lp > 2:
            ic.append(mol.getAngle(*path[:3]))
        if lp > 3:
            ic.append(mol.getDihedral(*path))
        icoords.append(ic)
    return icoords


def calc_cart_from_internal(zpath: list, icoords: list) -> ndarray:
    """从内坐标路径和内坐标计算笛卡尔坐标"""
    if zpath is None:
        zpath = [[0], [0, 1]] + [[i-2, i-1, i] for i in range(2, len(icoords))]
    coords = []
    for path, ic in zip(zpath, icoords):
        if len(path) == 1:
            pm = [0, 0, 0]
            coords.append(pm)
            continue
        elif len(path) == 2:
            mat = eye(3)
            pb = coords[path[0]]
            nb = ic[0]
            na = 180
            nd = 0
            vc = [1, 0, 0]
            va = [0, 1, 0]
            vb = [0, 0, 1]
        elif len(path) == 3:
            pb = coords[path[0]]
            pa = coords[path[1]]
            nb, na = ic
            nd = 0
            vc = pb - pa
            va = [0, 0, 1]
            vb = cross(vc, va)
        else:
            pb = coords[path[0]]
            pa = coords[path[1]]
            pd = coords[path[2]]
            nb, na, nd = ic
            vc = pb - pa
            va = pd - pa
            vb = cross(vc, va)
            va = cross(vb, vc)
        mat = array([va, vb, vc])
        mat = mat / linalg.norm(mat, axis=1).reshape((-1, 1))
        phi = deg2rad(180 - na)
        theta = deg2rad(nd)
        rm = array([sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)]) * nb
        pm = dot(rm, mat) + pb
        coords.append(pm)
    return array(coords)
