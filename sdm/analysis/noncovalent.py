# -*- coding:utf-8 -*-
from __future__ import annotations
from typing import Optional
import networkx as nx


HDonor = {"N", "O"}
HAcceptor = {"B", "N", "O", "F", "S", "Cl"}
XAtom = {"I", "Br", "Cl"}
XAcceptor = {"N", "O", "S"}
HeterocyclicAtom = {"N", "O", "S", "P"}


def find_pipi_sites(
        mol,
        r_bonds: Optional[list] = None,
        break_conjugate: Optional[list] = None) -> list[list]:
    """分子中寻找pipi堆积位点

    r_bonds: 柔性键列表
    break_conjugate: 手动断开共价的原子
    """
    if r_bonds is None:
        r_bonds = mol.findRotatableBonds()
    d_bonds = [(i, j) for i, j, t in mol.bonds(data='type') if t in {2, 3, 4}]
    d_atoms = [i for b in d_bonds for i in b]
    G = mol.graph.copy()
    for i, e in mol.atoms(data='element'):
        if all(j in d_atoms for j in G.adj[i]) and e.symbol in HeterocyclicAtom:
            d_atoms.append(i)
    G.remove_edges_from(r_bonds)
    G.remove_nodes_from(break_conjugate)
    G = G.subgraph(d_atoms)
    aroma_groups = [list(g) for g in nx.connected_components(G)]
    return aroma_groups


def find_HB_sites(mol) -> tuple[list, list]:
    """分子中寻找氢键位点"""
    donors = []
    acceptors = []
    adj = mol.graph.adj
    for i, e in mol.atoms(data='element'):
        if e.symbol in HDonor:
            donors += [(i, h) for h in adj[i] if mol.atoms[h]['element'].symbol == "H"]
        if e.symbol in HAcceptor and len(adj[i]) < 4:
            acceptors.append(i)
    return donors, acceptors


def find_XB_sites(mol) -> tuple[list, list]:
    """分子中寻找卤键位点"""
    donors = []
    acceptors = []
    adj = mol.graph.adj
    for i, e in mol.atoms(data='element'):
        if e.symbol in XAtom:
            donors += [(d, i) for d in adj[i]]
        if e.symbol in XAcceptor:
            acceptors.append(i)
    return donors, acceptors
