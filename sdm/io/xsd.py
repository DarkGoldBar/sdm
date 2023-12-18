# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2020/9/23
from __future__ import absolute_import, unicode_literals, division, print_function
import numpy as np
from xml.etree import ElementTree as etree


def parser(string, height2D=30.0):
        root = etree.fromstring(string)
        system_node = root.find("AtomisticTreeRoot/SymmetrySystem")
        title = system_node.attrib['Name']

        data_node = system_node.find("MappingSet/MappingFamily/IdentityMapping")
        dimension = 3
        symmetry_grp = data_node.find("SpaceGroup")
        if symmetry_grp is None:
            dimension = 2
            symmetry_grp = data_node.find("PlaneGroup")
        if symmetry_grp is None:
            dimension = 0
            matrix = np.eye(3, dtype=float)
        else:
            matrix = np.array([
                np.array(symmetry_grp.attrib['AVector'].split(','), float),
                np.array(symmetry_grp.attrib['BVector'].split(','), float),
                np.array(symmetry_grp.attrib['CVector'].split(','), float)
            ])

        atom_list = []
        for node in data_node.findall('Atom3d'):
            atom_list.append([
                int(node.attrib["ID"]),
                node.attrib["Components"],
                np.array(node.attrib["XYZ"].split(','), float),
                node.attrib.get("Name"),
                node.attrib.get("ForcefieldType", "xx"),
                node.attrib.get("Charge", "0.0")
            ])
        atoms_prop = list(zip(*sorted(atom_list)))
        fcoords = np.vstack(atoms_prop[2])
        coords = np.dot(fcoords, matrix)
        if dimension == 2:
            matrix[2] *= height2D + np.max(fcoords[:, 2]) - np.min(fcoords[:, 2])

        return {"title": title,
                "coords": coords,
                "atom_symbols": atoms_prop[1],
                "atom_titles": atoms_prop[3],
                "hall_number": 1,
                "lattice": matrix}

def writer(info_dict):
    pass
