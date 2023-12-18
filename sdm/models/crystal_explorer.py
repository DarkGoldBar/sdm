# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2021/8/10
from dataclasses import dataclass, field
from typing import List
import numpy as np
import networkx as nx
from sdm.analysis.find_bonds import detectDistance
from sdm.models.molecule import Molecule
from sdm.models.structure import Structure
from sdm.models.symmetry import SpaceGroup, SymmOp

class CrystalExplorer(Molecule):
    """用于分析晶体非键作用的分子团簇, 复现CrystalExplorer的操作

    Atom Properties:
        cell_x
        cell_y
        cell_z
        cell_sym
        cell_eq
    """

    def __init__(self,
                 struct: Structure,
                 unitcell: Structure = None,
                 asymmetric_unit: Structure = None):
        self._graph = struct.graph.copy()
        self._coords = struct.coords.copy()
        self._lattice = struct.lattice.copy()
        self._unitcell = unitcell
        self._asymmetric_unit = asymmetric_unit
        self.title = struct.title

        if self._asymmetric_unit is None and self._unitcell is None:
            self.initialize(struct)

    @property
    def graph(self):
        return self._graph

    @property
    def coords(self):
        return self._coords

    @property
    def lattice(self):
        return self._lattice

    @property
    def unitcell(self):
        return self._unitcell

    @property
    def asymmetric_unit(self):
        return self._asymmetric_unit

    def initialize(self, struct: Structure):
        norm = struct.copy()
        for n in norm.graph.nodes:
            data = {'cell_x': 0, 'cell_y': 0, 'cell_z': 0, 'cell_sym': 0, 'cell_eq': n}
            norm.graph.nodes[n].update(data)
        norm.setLatticeLowerTriangular()
        norm.moveMolCenterInCell()
        norm.space_group.sort()
        self._graph = norm.graph
        self._coords = norm.coords
        self._lattice = norm.lattice
        self._asymmetric_unit = norm
        self._unitcell = self.generateUnitCell()

    def generateUnitCell(self):
        unitcell = self._asymmetric_unit.copy()
        unitgraph = unitcell.graph.copy()
        new_frac_coords = [unitcell.frac_coords]
        for i, op in enumerate(unitcell.space_group.symmetry_ops[1:]):
            new_graph = unitgraph.copy()
            new_fcoord = op.operate_multi(unitcell.frac_coords)
            new_frac_coords.append(new_fcoord)
            for n in new_graph.nodes:
                new_graph.nodes[n]['cell_sym'] = i+1
            unitcell._addGraph(new_graph)
        unitcell.frac_coords = np.vstack(new_frac_coords)
        unitcell.moveAtomsByMolecule()
        unitcell.moveMolCenterInCell()
        remove_index = unitcell.delDuplicatedAtoms()
        if len(remove_index) > 0:
            print('[getP1Structure] atoms removed:', len(remove_index))

        for n, fc in zip(unitcell.graph.nodes, unitcell.frac_coords):
            c = np.floor(fc).astype(int)
            unitcell.graph.nodes[n].update({'cell_x': -c[0], 'cell_y': -c[1], 'cell_z': -c[2]})
        unitcell.moveAtomsInCell()
        unitcell.findBonds(None)
        unitcell.findBonds()
        return unitcell

    def generateSupercell(self, index, index_low):
        unitcell = self.unitcell
        supercell = unitcell.copy()
        for i in range(index_low[0], index[0]+1):
            for j in range(index_low[1], index[1]+1):
                for k in range(index_low[2], index[2]+1):
                    if (i, j, k) != (0, 0, 0):
                        new_graph = unitcell.graph.copy()
                        for n in new_graph:
                            new_graph.nodes[n]['cell_x'] += i
                            new_graph.nodes[n]['cell_y'] += j
                            new_graph.nodes[n]['cell_z'] += k
                        supercell._addGraph(new_graph)
                        supercell._addCoords(unitcell.coords + np.dot([i, j, k], unitcell.lattice.matrix))
        supercell.space_group = SpaceGroup.fromIntNumber(1)
        supercell.findBonds()
        return supercell

    def generateAtomsWithinRadius(self, radius=3.80, source: List[int] = None):
        unitcell = self.unitcell
        bbox = self.getBorderBox(source, expand=radius)
        bbox_frac = unitcell.lattice.getFractionalCoords(bbox)
        bbox_max = [int(x) for x in np.floor(np.max(bbox_frac, axis=0))]
        bbox_min = [int(x) for x in np.floor(np.min(bbox_frac, axis=0))]
        supercell = self.generateSupercell(bbox_max, bbox_min)
        coords0 = self.coords if source is None else self.coords[source]
        coords1 = supercell.coords
        radius0 = np.ones(len(coords0)) * radius
        radius1 = np.zeros(len(coords1))
        col, vec = detectDistance(coords0, coords1, radius0, radius1)
        in_range = sorted(np.unique(col[:, 1]))
        new_struct = supercell.getAtoms(in_range)
        new_object = self.__class__(new_struct, unitcell, self.asymmetric_unit)
        return new_object

    def generateCompletedFragments(self, step=5.0, limit=10):
        cur = self
        for i in range(limit):
            bbox = cur.getBorderBox(expand=step)
            bbox_frac = cur.lattice.getFractionalCoords(bbox)
            bbox_max = [int(x) for x in np.floor(np.max(bbox_frac, axis=0))]
            bbox_min = [int(x) for x in np.floor(np.min(bbox_frac, axis=0))]
            supercell = cur.generateSupercell(bbox_max, bbox_min)
            connected = set()
            for i in cur.graph:
                for j in supercell.graph:
                    if cur.atoms[i] == supercell.atoms[j]:
                        c = nx.connected.node_connected_component(supercell.graph, j)
                        connected = connected | c
                        break
                else:
                    raise ValueError(f'node {i} find no match in new graph')
            if len(connected) == len(cur.graph):
                break
            new = supercell.getAtoms(sorted(connected))
            cur = self.__class__(new, self.unitcell, self.asymmetric_unit)
        return cur

    def generateDimers(self, groups=None):
        groups = groups if groups else self.asymmetric_unit.groups

        unit_map = {}
        for i, d in self.atoms.data():
            unit_key = d['cell_x'], d['cell_y'], d['cell_z'], d['cell_sym']
            atom_key = d['cell_eq']
            if unit_key not in unit_map:
                unit_map[unit_key] = {}
            if atom_key in unit_map[unit_key]:
                raise KeyError(f'equivalent atom ({atom_key}: {i}) exists in map {unit_key}:{unit_map[unit_key]}')
            unit_map[unit_key][atom_key] = i

        group_map = {}
        asym_groups = []
        for ukey, unit in unit_map.items():
            for ig, g in enumerate(groups):
                atoms = [unit.get(x) for x in g]
                if any(x is None for x in atoms): continue
                cent = self.getCentroid(atoms)
                rot = self.asymmetric_unit.space_group.symmetry_ops[ukey[3]].rotation_matrix
                op = SymmOp.from_rotation_and_translation(rot, cent)
                gkey = ukey + (ig, )
                if ukey == (0, 0, 0, 0):
                    asym_groups.append(gkey)
                group = {'atoms': atoms, 'op': op}
                group_map[gkey] = group

        dimers = []
        for agkey in asym_groups:
            agroup = group_map[agkey]
            mol1 = self.getAtoms(agroup['atoms'])
            mol1.title += str(agkey)
            for gkey, group in group_map.items():
                if gkey == agkey: continue
                mol2 = self.getAtoms(group['atoms'])
                mol2.title += str(gkey)
                dimer = CEDimer(
                    mol1=mol1,
                    mol2=mol2,
                    op1=agroup['op'],
                    op2=group['op'],
                    frag1=agkey[-1],
                    frag2=gkey[-1],
                    tag=gkey[:-1],
                )
                dimers.append(dimer)

        return dimers

    def getAtoms(self, groups):
        title = self.asymmetric_unit.title
        return Molecule(self.graph, self.coords, title=title).getAtoms(groups)

@dataclass
class CEDimer:
    mol1: Molecule
    mol2: Molecule
    op1: SymmOp
    op2: SymmOp
    frag1: int = 0
    frag2: int = 0
    tag: any = 'NoName'
    op12: SymmOp = field(init=False)

    def __repr__(self):
        return f'CEDimer(tag={self.tag}, frag1={self.frag1}, frag2={self.frag2})'

    def __post_init__(self):
        self.op12 = self.op1.inverse * self.op2

    def equal(self, other, precision=1e-2):
        if self.frag1 == other.frag1 and self.frag2 == other.frag2:
            if np.isclose(self.op12.affine_matrix, other.op12.affine_matrix, precision, precision).all():
                return True
        if self.frag2 == other.frag1 and self.frag1 == other.frag2:
            if np.isclose(self.op12.inverse.affine_matrix, other.op12.affine_matrix, precision, precision).all():
                return True
        return False

    def mindist(self):
        coords0 = self.mol1.coords
        coords1 = self.mol2.coords
        dist2 = -2 * np.dot(coords0, coords1.T) \
            + np.sum(coords0**2, axis=1).reshape((-1, 1)) \
            + np.sum(coords1**2, axis=1)
        min_dist = np.sqrt(np.min(dist2))
        return min_dist

    @staticmethod
    def filter_distance(dimer_list, min_dist=1.5, max_dist=5.0):
        return [d for d in dimer_list if min_dist < d.mindist() < max_dist]
            
    @staticmethod
    def filter_symmetry(dimer_list):
        passed = []
        grouped = {}
        for dimer in dimer_list:
            for i, d in enumerate(passed):
                if dimer.equal(d):
                    if i not in grouped:
                        grouped[i] = []
                    grouped[i].append(dimer)
                    break
            else:
                passed.append(dimer)
        return passed, grouped

    @staticmethod
    def plot(dimer_list: List["CEDimer"]):
        from sdm.analysis.plotly import Vis
        checked0 = CEDimer.filter_distance(dimer_list)
        checked, grouped = CEDimer.filter_symmetry(checked0)

        asymm = {}
        for dimer in dimer_list:
            if dimer.frag1 not in asymm:
                asymm[dimer.frag1] = dimer.mol1

        v = Vis()
        for i, mol in asymm.items():
            v.add_mol(mol)

        for i, d in enumerate(checked):
            if d.tag == (0,0,0,0): continue
            h = i / (len(checked)+1) * 300
            hsv = f'hsv({h},90%,50%)'    
            text =  f'{i}:{d.tag}:{{i}}:{{t}}'
            v.add_mol(d.mol2, atom_color=hsv, bond_color=hsv, text=text)
            if i in grouped:
                for dd in grouped[i]:
                    text = f'{i}:{dd.tag}:{{i}}:{{t}}'
                    v.add_mol(dd.mol2, atom_color=hsv, bond_color=hsv, text=text)

        v._plot()
        return v