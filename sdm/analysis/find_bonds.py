# -*- coding: utf-8 -*-
import numpy as np
from ..utils import latt62Matrix


def find_empiric_bonds(coords, atoms_covalent_radius, detect_cut=1.25):
    """经验找键

    Args:
        coords:
        atoms_covalent_radius:
        detect_cut:

    Returns:
        bonds(numpy.array)
    """
    radius = atoms_covalent_radius * detect_cut
    col, _ = detectDistance(coords, coords, radius, radius)
    unique_col = np.unique(col[col[:, 0] < col[:, 1]], axis=1)
    bonds = [(int(i), int(j), 1) for i, j in unique_col]  # numpy.int64 -> int
    return bonds


def find_babel_bonds(xyz_string):
    """find bond by babel (openbabel 2.x+)

    Args:
        coords:原子坐标
        atom_symbols:原子符号

    Returns:
        numpy.dnarray: bonds

    """
    from comlink.processor.babel import Babel
    from ..io.mol import parser
    babel = Babel()
    mol_string = babel.convert(xyz_string, "xyz", "mol")
    info_dict = parser(mol_string)
    return info_dict["bonds"]


def find_empiric_bonds_pbc(fcoords, atoms_covalent_radius, lattice_param, affine_matrices, detect_cut=1.25):
    radius = atoms_covalent_radius * detect_cut
    col, vec = ddosap_py(fcoords, radius, lattice_param, affine_matrices)
    halfbonds = {}
    for (c1, c2), vec in zip(col, vec):
        if c1 in halfbonds:
            halfbonds[c1][c2] = vec
        else:
            halfbonds[c1] = {c2: vec}
    return halfbonds


def detectDistance(coords0, coords1, radius0, radius1, symprec=1e-4):
    """
    Args:
        coords0: (float[n * 3]) 原子坐标0
        coords1: (float[m * 3]) 原子坐标1
        radius0: (float[n]) 原子坐标0中, 原子的半径
        radius1: (float[m]) 原子坐标1中, 原子的半径
        symmprec: (float) 距离低于这个值判断为相同原子

    Return:
        col (int[a * 2]) 碰撞原子对
        vec (float[a * 3]) 碰撞原子对的向量0->1
    """
    radius2 = (radius0.reshape((-1, 1)) + radius1)**2
    dist2 = -2 * np.dot(coords0, coords1.T) \
        + np.sum(coords0**2, axis=1).reshape((-1, 1)) \
        + np.sum(coords1**2, axis=1)
    mask = (symprec**2 < dist2) * (dist2 < radius2)
    col = np.vstack(np.where(mask)).T
    vec = coords1[col[:, 1]] - coords0[col[:, 0]]
    return col, vec


def ddosap_py(fcoords, radius, latt6, affine_matrix, symprec=1e-4):
    """包含周期性对称性的原子碰撞检测
    Args:
        fcoords: (float[n * 3]) 非对称单元的原子的分数坐标
        radius: (float[n]) 非对称单元的原子的半径
        latt6: (float[6]) 晶胞参数
        affine_matrix: (float[m * 4 * 4]) 对称操作的仿射变换矩阵列表
        symmprec: (float) 距离低于这个值判断为相同原子

    Usage:

    def ddosap(struct, factor=1.25):
        return ddosap_py(
            struct.frac_coords,
            np.array([a.covelent_radius for a in struct.elements]) * factor,
            struct.lattice.parameters,
            [op.affine_matrix for op in struct.space_group.symmetry_ops],
        )
    """
    nsym = len(affine_matrix)
    natom = len(radius)
    assert natom == len(fcoords)
    matrix = latt62Matrix(latt6)
    matinv = np.linalg.inv(matrix)
    cube = np.array([matrix[0, 0], matrix[1, 1], matrix[2, 2]])

    # expand symmerty
    fc4 = np.hstack([fcoords, np.ones((natom, 1))])
    fcoords2 = np.vstack([np.inner(fc4, affine)[:, :3] for affine in affine_matrix])
    fcoords2 -= np.floor(fcoords2)

    # PBC limit
    r_max = np.max(radius) * 2
    vertex = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                       [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
    cvertex_i = vertex * cube
    cvertex_o = cvertex_i + r_max * vertex + r_max * (vertex - 1)
    fvertex = np.dot(cvertex_o, matinv)
    pbc_min = np.floor(np.min(fvertex, axis=0)).astype(int)
    pbc_max = np.ceil(np.max(fvertex, axis=0)).astype(int)
    ncell = (pbc_max - pbc_min) + 1
    ncell = ncell[0] * ncell[1] * ncell[2]

    # PBC expand
    foffset = np.array(np.meshgrid(np.arange(pbc_min[0], pbc_max[0] + 1),
                                   np.arange(pbc_min[1], pbc_max[1] + 1),
                                   np.arange(pbc_min[2], pbc_max[2] + 1))).T.reshape(-1, 3)
    coords3 = np.tile(np.dot(fcoords2, matrix), (ncell, 1))
    coords3 += np.tile(np.dot(foffset, matrix), nsym * natom).reshape((-1, 3))

    # coords3: atoms out cube, in cube+r_max
    mask3 = np.all((cvertex_o[0] <= coords3), axis=1) * np.all((coords3 < cvertex_o[-1]), axis=1)
    coords3 = coords3[mask3]
    index3 = np.tile(np.arange(natom), nsym * ncell)[mask3]
    radius3 = radius[index3]
    # coords2: atoms in cube
    mask2 = np.all((cvertex_i[0] <= coords3), axis=1) * np.all((coords3 < cvertex_i[-1]), axis=1)
    coords2 = coords3[mask2]
    index2 = index3[mask2]
    radius2 = radius3[mask2]
    assert coords2.shape[0] == nsym * natom
    # collision test
    idx, vec = detectDistance(coords2, coords3, radius2, radius3, symprec=symprec)
    idx[:, 0] = index2[idx[:, 0]]
    idx[:, 1] = index3[idx[:, 1]]
    return idx, vec
