# -*- coding: utf-8 -*-
from __future__ import annotations
from numpy import ndarray

import numpy as np
from ..models.space import Space


def kabsch(P: ndarray, Q: ndarray) -> ndarray:
    """Kabsch algorithm
    https://en.wikipedia.org/wiki/Kabsch_algorithm
    """
    assert P.shape == Q.shape, "P, Q must have same shape"
    H = np.dot(P.T, Q)
    U, S, V = np.linalg.svd(H)
    d = np.linalg.det(V) * np.linalg.det(U)
    if d < 0:
        U[:, -1] *= -1
    R = np.dot(U, V)
    return R


def kabsch_w(P: ndarray, Q: ndarray, w: ndarray) -> ndarray:
    """weighted-Kabsch algorithm
    https://en.wikipedia.org/wiki/Wahba%27s_problem
    """
    assert P.shape == Q.shape, "P, Q must have same shape"
    assert P.shape[0] == len(w), "P, Q, w must have same length"
    H = np.dot(P.T, Q * w)
    U, S, V = np.linalg.svd(H)
    d = np.linalg.det(V) * np.linalg.det(U)
    if d < 0:
        V[:, -1] *= -1
    R = np.dot(V, U.T)
    return R


def rmsd(V: ndarray, W: ndarray) -> float:
    """Root-mean-square deviation
    https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
    """
    assert V.shape == W.shape, "V, W must have same shape"
    return np.sqrt(np.sum((V - W) ** 2) / V.shape[0])


class MoleculeComparsion(Space):
    def __init__(self, entity1, entity2):
        self.entity1 = entity1
        self.entity2 = entity2
        self.rmsd = None
        super().__init__(np.vstack([entity1.coords, entity2.coords]))

    @classmethod
    def compare(cls, mol_a, mol_b, without_H=True):
        group = None
        if without_H:
            group = [i for i, el in enumerate(mol_a.elements) if el.name != "H"]
            coords_a = mol_a.coords[group]
            coords_b = mol_b.coords[group]
        else:
            coords_a = mol_a.coords
            coords_b = mol_b.coords
        cent_a = mol_a.getGeometryCenter(group)
        cent_b = mol_b.getGeometryCenter(group)
        rotate = kabsch(coords_a - cent_a, coords_b - cent_b)
        mol_b2 = mol_b.copy()
        mol_b2.coords = np.dot(mol_b2.coords - cent_b, rotate) + cent_a

        result = cls(mol_a, mol_b2)
        result.rmsd = rmsd(coords_a, mol_b2.coords[group] if without_H else mol_b2.coords)
        return result
