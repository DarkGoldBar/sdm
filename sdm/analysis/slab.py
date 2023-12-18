# -*- coding: utf-8 -*-
# Author : bochen.li
# Date : 2020/11/19
from math import gcd
from itertools import product
import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.surface import SlabGenerator
import sdm


def cutAllSlabs(struct, max_miller=1, min_thickness=10, min_vacuum=30, ftol=0.05, tol=0.1,
                primitive=True, in_cell=False):
    slabs = []
    for h, k, l in product(*(range(max_miller + 1), ) * 3):
        if gcd(gcd(h, k), l) != 1:
            print('{} skiped'.format((h, k, l)))
            continue
        print('{} processing'.format((h, k, l)))
        hkl_list = [(h, k, l)]
        if h != 0:
            hkl_list += [(-h, k, l) for h, k, l in hkl_list]
        if k != 0:
            hkl_list += [(h, -k, l) for h, k, l in hkl_list]
        if l != 0:
            hkl_list += [(h, k, -l) for h, k, l in hkl_list]
        hkl_slab = []
        for hkl in hkl_list:
            hkl_slab += cutOneSlab(struct, hkl, min_thickness, min_vacuum,
                                   ftol=ftol, tol=tol,
                                   primitive=primitive, in_cell=in_cell)
        slabs += removeDuplicatedSlabs(hkl_slab)
    return slabs


def removeDuplicatedSlabs(slabs, tol=1e-2):
    remove = set()
    sm = StructureMatcher(ltol=tol, stol=tol, primitive_cell=False, scale=False)
    mg_slabs = [slab.toMgStruct() for slab in slabs]
    for i, si in enumerate(mg_slabs):
        if i in remove:
            continue
        for j, sj in enumerate(mg_slabs[i+1:]):
            if j+i+1 in remove:
                continue
            rms = sm.get_rms_dist(si, sj)
            if (rms is not None) and (rms[0] < tol):
                remove.add(j+i+1)
    return [s for i, s in enumerate(slabs) if i not in remove]


def cutOneSlab(struct, hkl, min_thickness, min_vacuum, ftol=0.05, tol=0.1,
               primitive=True, in_cell=False):
    mgst = struct.toMgStruct(False)
    slab_gen = SlabGenerator(mgst, hkl, 1, 0, primitive=False)
    # 初始切片
    slab0 = slab_gen.get_slab(0)
    st0 = sdm.open("cif", slab0.to("cif"))
    st0.setLatticeLowerTriangular()
    st0.moveAtomsByMolecule()
    st0.findBonds(None)
    st0.findBonds('empiric')
    st0.moveMolCenterInCell()

    # 计算 shift
    shifts = []
    cart_centers = [st0.getGeometryCenter(g) for g in st0.groups]
    frac_centers = st0.lattice.getFractionalCoords(cart_centers)
    last = -1.0
    for c in sorted(frac_centers[:, 2]):
        if c - last > ftol:
            shifts.append(c)
            last = c

    # 计算 slabs
    slabs = []
    cell2x = st0.lattice.matrix * np.array([[1, 1, 2]]).T
    for s in shifts:
        slab = st0.copy()
        shift = [0, 0, s + ftol / 10]
        slab.frac_coords = slab.frac_coords - shift
        slab.moveMolCenterInCell()
        slab.setLattice(cell2x)
        slabs.append(slab)

    # 去重
    slabs_nodap = removeDuplicatedSlabs(slabs, tol=tol)

    # 展开厚度
    hc = st0.lattice.height[2]
    ncell = int(min_thickness // hc) + 1
    nvacuum = int(min_vacuum // hc) + 1
    new_cell = st0.lattice.matrix * np.array([[1, 1, ncell + nvacuum]]).T
    vec_c = st0.lattice.matrix[2]
    for slab in slabs_nodap:
        g = slab.graph.copy()
        c = slab.coords.copy()
        for i in range(1, ncell):
            slab.addGraph(g)
            slab.addCoords(c + i * vec_c)
        slab.setLattice(new_cell)

    # 转化为原胞
    if primitive:
        slabs_prim = []
        for slab in slabs_nodap:
            prim = slab.getPrimitiveStructure()
            prim.space_group = sdm.models.SpaceGroup.fromIntNumber(1)
            prim.setLatticeLowerTriangular()
            prim.moveAtomsByMolecule()
            prim.findBonds(None)
            prim.findBonds('empiric')
            prim.moveMolCenterInCell()
            slabs_prim.append(prim)
    else:
        slabs_prim = slabs_nodap

    # 沿着z轴平移到底面
    if in_cell:
        for slab in slabs_prim:
            fcoords = slab.frac_coords
            shift = min(fcoords[:, 2])
            fcoords[:, 2] += (0.001 - shift)
            slab.frac_coords = fcoords

    # 标记
    hkl_string = str(tuple(hkl)).replace(' ', '')
    for i, slab in enumerate(slabs_prim):
        slab.hkl = tuple(hkl)
        slab.hkli = slab.hkl + (i, )
        slab.title = struct.title + f"-slab{hkl_string}({i})"
    return slabs_prim
