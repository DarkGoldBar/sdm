# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from .find_bonds import detectDistance
from ..models import Molecule

# IDATM
# https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/idatm.html
# https://doi.org/10.1002/jcc.540120716


def find_idatm_bonds(atom_symbols, coords):
    mol = Molecule.fromDict({'atom_symbols': atom_symbols, 'coords': coords}, None)
    # 1. 经验找键
    bonds = empiric(mol)
    mol.setBonds(bonds)

    # 2. 设置属性 di, 根据 di 校正连接性
    remove = atom_info_di(mol)
    mol.delBonds(remove)

    # 3.1 杂化判断, 设置属性 hybr
    # 3.2 根据杂化类型, 设置单键 (有一个原子有杂化的键, 设置为单键)
    bonds = assign_initial_hybridization(mol)
    mol.setBonds(bonds)

    # 3.3 根据柔性角设置单键
    bonds = assign_single_bond_use_dihe(mol)
    mol.setBonds(bonds)

    # 3.4 根据键长判断单键
    bonds = assign_single_bond_use_bond_length(mol)
    mol.setBonds(bonds)

    # 3.5 把所有未标识原子且它的所有键连都为单键, 标识为 sp3
    bonds = percieve_hybr_and_bond_order(mol)
    mol.setBonds(bonds)

    # 3.6 把没有匹配的键统一设置为 single bond
    bonds = [[i, j, 1] for i, j in mol.bonds if mol.bonds[(i, j)]['type'] == 0]
    mol.setBonds(bonds)

    # 4. 特殊情况校正 (S=O(=O))
    bonds = correct_bonds_by_groups(mol)
    mol.setBonds(bonds)
    return mol.bonds(data='type')


def empiric(mol) -> list[int, int, int]:
    """1. 经验找键
    0.1 < rij < Ri + Rj + 0.4,  rij = |xi - xj|, Ri = Covalent radius
    """
    radius = np.array([ele.covelent_radius for ele in mol.elements]) + 0.4
    col, _ = detectDistance(mol.coords, mol.coords, radius, radius, 0.1)
    unique_col = np.unique(col[col[:, 0] < col[:, 1]], axis=1)
    bonds = [(i, j, 0) for i, j in unique_col]
    return bonds


def atom_info_di(mol) -> list[int, int]:
    """2. 设置属性 di, 根据 di 校正连接性
    """
    adj = mol.graph.adj
    for i in mol.atoms:
        if len(adj[i]) < 2:
            di = len(adj[i])
        else:
            # 计算原子与邻接原子的整体协方差矩阵 (covariance matrix)
            c = mol.coords[list(adj[i]) + [i]]
            c -= np.average(c, axis=0)
            cov_mat = sum([np.outer(p, p) for p in c])
            eig_values = np.linalg.eigvals(cov_mat)
            # square root >= 0.2
            di = sum(eig_values > 0.04)
        mol.atoms[i]['di'] = di

    remove = []
    for i in mol.atoms:
        max_conn = determine_atom_max_bond_number(mol.atoms[i]['di'], mol.elements[i].Z)
        k = len(adj[i])
        if k > max_conn:
            neighbor_dist = [(mol.getDistance(i, inb), inb) for inb in adj[i]]
            remove_list = list(zip(*sorted(neighbor_dist)))[1][max_conn:]
            remove += [[i, inb] for inb in remove_list]
    return remove


def assign_initial_hybridization(mol: Molecule) -> list[int, int, int]:
    """3.1 杂化判断, 设置属性 hybr
    探测明显的杂化类型
    这一步完了以后, 就只剩下 d < 3, z = {C, N, O, Si, P, S, Se}, Q < 4, 没有处理. (在有机里面, 相当于根本没处理)
    """
    adj = mol.graph.adj
    # 杂化初始猜测
    hybr_list = [determine_atom_hybr(
        mol.elements[i].Z,  # 原子序号
        len(adj[i]),  # 临接原子数
        mol.atoms[i]['di']      # di
    ) for i in mol.atoms]

    # 所有的邻接原子杂化已经标定, 则标定本原子为 sp3
    find_new_atom = True
    while find_new_atom:
        find_new_atom = False
        for i in mol.atoms:
            if hybr_list[i] == '?' and all(hybr_list[idx] != '?' for idx in adj[i]):
                hybr_list[i] = 'sp3'
                find_new_atom = True
    mol.setAtomAttribute('hybr', hybr_list)
    # 3.2 根据杂化类型, 设置单键 (有一个原子有杂化的键, 设置为单键)
    bonds = []
    for i, j in mol.bonds:
        if hybr_list[i] != '?' or hybr_list[j] != '?':
            bonds.append([i, j, 1])
    return bonds


def assign_single_bond_use_dihe(mol: Molecule) -> list[int, int, int]:
    """3.3 根据柔性角标识单键
    """
    adj = mol.graph.adj
    bonds = []
    for i, j in mol.bonds:
        if mol.atoms[i]['di'] == 1 or mol.atoms[j]['di'] == 1:
            continue
        dihe_list = []
        for ni in set(adj[i]) - {j}:
            for nj in set(adj[j]) - {i}:
                dihe = abs(mol.getDihedral(ni, i, j, nj))
                dihe_list.append(min(dihe, 180 - dihe))
        if dihe_list and max(dihe_list) >= 15.0:
            bonds.append([i, j, 1])
    return bonds


def assign_single_bond_use_bond_length(mol: Molecule) -> list[int, int, int]:
    """3.4 根据键长判断单键
    """
    bonds = []
    for i, j in mol.bonds:
        if mol.bonds[i, j]['type'] != 0:
            continue
        rij = mol.getDistance(i, j)
        label = "-".join([e.symbol for e in sorted([mol.elements[i], mol.elements[j]], key=lambda e:e.Z)])
        if rij > BONDLENGTH[label] - 0.05:
            bonds.append([i, j, 1])
    return bonds


def percieve_hybr_and_bond_order(mol: Molecule) -> list[int, int, int]:
    for i in mol.atoms:
        if mol.atoms[i]['hybr'] == '?' and all(mol.bonds[b]['type'] == 1 for b in mol.bonds(i)):
            mol.atoms[i]['hybr'] = 'sp3'
    # 双键,三键分析
    atoms_weight = [calc_atom_weight(mol, i) for i in mol.atoms]
    bonds_weight = []
    for i, j in mol.bonds:
        if mol.bonds[(i, j)]['type'] == 1:
            continue
        wij = atoms_weight[i] + atoms_weight[j]
        rij = mol.getDistance(i, j)
        label = "-".join([e.symbol for e in sorted([mol.elements[i], mol.elements[j]], key=lambda e:e.Z)])
        if rij < BONDLENGTH[label] - 0.25:
            wij += 3.0*(BONDLENGTH[label] - rij)  # 键长越短, 权重越大
            bonds_weight.append((i, j, wij))
        elif rij < BONDLENGTH[label] - 0.11:
            wij += 2.0*(BONDLENGTH[label] - rij)  # 键长越短, 权重越大
            bonds_weight.append((i, j, wij))

    match_bonds = mwma(bonds_weight)
    total_w = 0.
    for i, j in match_bonds:
        for idx1, idx2, w in bonds_weight:
            if (i == idx1 and j == idx2):
                total_w += w
                break

    adj = mol.graph.adj
    bonds = []
    for i, j in match_bonds:
        di = mol.atoms[i]['di']
        dj = mol.atoms[j]['di']
        # 线性原子处理
        if (di == 1 and len(adj[i]) == 2) or (dj == 1 and len(adj[j]) == 2):
            if di == 1 and len(adj[i]) == 2:
                linear_atom = i
                nbs = list(adj[i])
            elif dj == 1 and len(adj[j]) == 2:
                linear_atom = j
                nbs = list(adj[j])
            i, j, k = nbs[0], linear_atom, nbs[1]
            rij = mol.getDistance(i, j)
            rjk = mol.getDistance(j, k)
            label_ij = "-".join([e.symbol for e in sorted([mol.elements[i], mol.elements[j]], key=lambda e:e.Z)])
            label_jk = "-".join([e.symbol for e in sorted([mol.elements[j], mol.elements[k]], key=lambda e:e.Z)])
            if label_ij in BONDLENGTH.keys() and rij < BONDLENGTH[label_ij] - 0.25:
                bonds.append([i, j, 3])
                bonds.append([j, k, 1])
            elif label_jk in BONDLENGTH.keys() and rjk < BONDLENGTH[label_jk] - 0.25:
                bonds.append([i, j, 1])
                bonds.append([j, k, 3])
            elif label_ij in BONDLENGTH.keys() and rij < BONDLENGTH[label_ij] - 0.11:
                bonds.append([i, j, 2])
                bonds.append([j, k, 2])
            elif label_jk in BONDLENGTH.keys() and rjk < BONDLENGTH[label_jk] - 0.11:
                bonds.append([i, j, 2])
                bonds.append([j, k, 2])
            else:
                print("!!!!!!!!!!!!!!!!!!!!!!! bond length two long!!!  %s: %s - %s" % (label_ij, rij, rjk))
                bonds.append([i, j, 1])
                bonds.append([j, k, 1])
        else:
            rij = mol.getDistance(i, j)
            label = "-".join([e.symbol for e in sorted([mol.elements[i], mol.elements[j]], key=lambda e:e.Z)])
            if rij < BONDLENGTH[label] - 0.25:
                # TODO: 处理线性原子的两个双键问题
                order = 3
            elif rij < BONDLENGTH[label] - 0.11:
                order = 2
            else:
                # raise Exception("bond length two long!!!")
                print("!!!!!!!!!!!!!!!!!!!!!!! %s - %s: %s - %s bond length two long!!!" % (i+1, j+1, rij, BONDLENGTH[label]))
                order = 1
            bonds.append([i, j, order])
    return bonds


def correct_bonds_by_groups(mol: Molecule):
    """4. 特殊情况校正

    情况说明:
      1. 四度 S 接 2 个 O, S-O 校正为双键.
        TODO:
          是否要对这一族都进行相同校正?
          是否要对键长进行限制?
      2. 四度 P 接 3 个 O, 其中一个 P-O 校正为双键
      2. H3PO4 (未实现)
    """
    adj = mol.graph.adj
    bonds = []
    s_oo_list = []
    p_ooo_list = []
    for i, e in mol.atoms('element'):
        if e.symbol == 'S':
            nbs = list(adj[i])
            o_nbs = []
            if len(nbs) == 4:
                for nb in nbs:
                    if mol.elements[nb].symbol == 'O':
                        o_nbs.append(nb)
                if len(o_nbs) == 2:
                    s_oo_list.append((i, o_nbs))

        if e.symbol == 'P':
            nbs = list(adj[i])
            if len(nbs) == 4:
                o_nbs = [nb for nb in nbs if mol.elements[nb].symbol == 'O']
                if len(o_nbs) == 3:
                    p_ooo_list.append((i, o_nbs))

    for s_oo in s_oo_list:
        s_idx, (o_idx1, o_idx2) = s_oo
        bonds.append([s_idx, o_idx1, 2])
        bonds.append([s_idx, o_idx2, 2])

    for p_ooo in p_ooo_list:
        p_idx, (o_idx1, o_idx2, o_idx3) = p_ooo
        r1 = mol.getDistance(p_idx, o_idx1),
        r2 = mol.getDistance(p_idx, o_idx2),
        r3 = mol.getDistance(p_idx, o_idx3)
        o_idx = sorted([(r1, o_idx1), (r2, o_idx2), (r3, o_idx3)])[0][1]
        bonds.append([p_idx, o_idx, 2])
    return bonds


def determine_atom_max_bond_number(di, zi):
    """判断最大成键数
    Args:
      di(int):
      zi(int): 原子数(即原子在周期表的排位)
    Return:
      int : 最大成键数
    """
    if di == 0:  # isolation atoms
        bi = 0
    elif zi < 3:  # H and He
        bi = 1
    elif zi > 2 and di == 1:  # sp or linear atoms
        bi = 2
    elif zi < 11 and di == 2:  # sp2, second line of periodic table
        bi = 3
    elif (zi > 10 and di == 2) or (zi < 11 and di == 3):  # sp3 or square-planar
        bi = 4
    else:
        bi = 7
    return bi


def determine_atom_hybr(zi, qi, di):
    """原子杂化类型判断
    Args:
      zi(int): 原子数(元素周期表序号)
      qi(int): 连接数(图结点的边数)
      di(int):
    """
    if zi in [1, 2]:
        hybr = 'sp3'
    elif (qi > 4 and (zi in GROUPTABLE[5])) or (qi == 5 and (zi in GROUPTABLE[4] + GROUPTABLE[5] + GROUPTABLE[6] + GROUPTABLE[7] + GROUPTABLE[8])):
        hybr = 'dsp3'
    elif (qi > 4 and (zi in GROUPTABLE[6])) or (qi == 6 and (zi in GROUPTABLE[4] + GROUPTABLE[5] + GROUPTABLE[6] + GROUPTABLE[7] + GROUPTABLE[8])):
        hybr = 'd2sp3'
    elif (qi > 4 and (zi in GROUPTABLE[7])) or (qi == 7 and (zi in GROUPTABLE[4] + GROUPTABLE[5] + GROUPTABLE[6] + GROUPTABLE[7] + GROUPTABLE[8])):
        hybr = 'd3sp3'
    elif qi == 4 and zi > 10 and di == 2:
        hybr = 'd2sp3'
    elif zi in TRANSITION_METALS:
        hybr = 'd2sp3'
    elif zi > 10 and zi not in [14, 15, 16, 34]:  # {Si, P, S, Se} = [14, 15, 16, 34]
        hybr = 'd2sp3' if qi > 4 else 'sp3'
    elif qi == 4 or (qi == 3 and di == 3):
        hybr = 'sp3'
    elif qi > 2 and zi in GROUPTABLE[6] + GROUPTABLE[7] + GROUPTABLE[8]:
        hybr = 'sp3'
    elif zi not in [6, 7, 8, 14, 15, 16, 34]:  # {C, N, O, Si, P, S, Se }
        hybr = 'sp3'
    else:
        hybr = '?'  # 默认未知
    return hybr


def calc_atom_weight(mol: Molecule, idx: int):
    """根据原子的连接方式计算原子权重 BFS
    """
    conn = len(mol.graph.adj[idx])
    weight = -20.0
    if conn > 3:
        return -20.0
    for key, value in WEIGHT_INFO:
        current = set([idx])
        visited = set()
        for k in key:
            matched = set(i for i in current if mol.elements[i].symbol == k)
            if not matched:
                break
            visited |= matched
            current = mol.getNeighbor(matched) - visited
        else:
            weight = value[conn]
            break
    return weight


def mwma(edges):
    """Max Weighted Matching
    """
    import networkx as nx
    from networkx.algorithms.matching import max_weight_matching
    G = nx.Graph()
    G.add_weighted_edges_from(edges)
    res = max_weight_matching(G)
    return sorted([sorted(item) for item in res])


GROUPTABLE = {
    1: [1, 3, 11, 19, 37, 55, 87],
    2: [4, 12, 20, 38, 56, 88],
    3: [5, 13, 31, 49, 81, 113],
    4: [6, 14, 32, 50, 82, 114],
    5: [7, 15, 33, 51, 83, 115],
    6: [8, 16, 34, 52, 84, 116],
    7: [9, 17, 35, 53, 85, 117],
    8: [2, 10, 18, 36, 54, 86, 118],
}

TRANSITION_METALS = list(range(21, 31)) + list(range(39, 49)) + list(range(72, 81)) + list(range(104, 109))

BONDLENGTH = {
    "C-C": 1.54, "C-N": 1.47, "C-O": 1.43, "C-Si": 1.86, "C-P": 1.85, "C-S": 1.75, "C-Se": 1.97,
    "N-N": 1.45, "N-O": 1.43, "N-Si": 1.75, "N-P": 1.68, "N-S": 1.76, "N-Se": 1.85,
    "O-O": 1.47, "O-Si": 1.63, "O-P": 1.57, "O-S": 1.57, "O-Se": 1.97,
    "Si-Si": 2.36, "Si-P": 2.26, "Si-S": 2.15, "Si-Se": 2.42,
    "P-P": 2.27, "P-S": 2.07, "P-Se": 2.27,
    "S-S": 2.05, "S-Se": 2.19,
    "Se-Se": 2.34
}

WEIGHT_INFO = [
    [["C", "O"], {1: 1.3, 2: 4.0, 3: 4.0}],
    [["C", "N"], {1: -6.9, 2: 4.0, 3: 4.0}],
    [["C"], {1: 0.0, 2: 4.0, 3: 4.0}],
    [["N", "C", "O"], {1: -2.4, 2: -0.8, 3: -7.0}],
    [["N", "C", "N"], {1: -1.4, 2: 1.3, 3: -3.0}],
    [["N"], {1: 1.2, 2: 1.2, 3: 0.0}],
    [["O", "C", "O"], {1: 4.2, 2: -8.1, 3: -20.0}],
    [["O", "C", "N"], {1: 4.2, 2: -8.1, 3: -20.0}],
    # [["O"], {1: 0.2, 2: -6.5, 3: -20.0}],
    [["O"], {1: 2.0, 2: -6.5, 3: -20.0}],  # 强制把末端 O 权重调高.

    [['Si', 'S'], {1: 1.2, 2: 3.9, 3: 3.9}],
    [['Si', 'P'], {1: -7.0, 2: 3.9, 3: 3.9}],
    [['Si'], {1: -0.1, 2: 3.9, 3: 3.9}],
    [['P', 'Si', 'S'], {1: -2.5, 2: -0.9, 3: -7.1}],
    [['P', 'Si', 'P'], {1: -1.5, 2: 1.2, 3: -3.1}],
    [['P'], {1: 1.1, 2: 1.1, 3: -0.1}],
    [['S', 'Si', 'S'], {1: 4.1, 2: -8.2, 3: -20.1}],
    [['S', 'Si', 'P'], {1: 4.1, 2: -8.2, 3: -20.1}],
    [['S'], {1: 0.1, 2: -6.6, 3: -20.1}]
]


if __name__ == '__main__':
    # python3 -m sdm.analysis.find_bonds_idatm
    import sdm
    st = sdm.open('/home/bochen.li/tmp/ACSALA.cif')
    d = st.toDict()
    print(find_idatm_bonds(d['atom_symbols'], d['coords']))
