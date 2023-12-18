# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2021/08/20
from ..utils import Units


def parser(string) -> dict:
    ss = iter(string.splitlines())
    for line in ss:
        if '$cell vectors' in line:
            matrix = []
            for i in range(3):
                line = next(ss)
                matrix.append([float(x) / Units.bohr for x in line.split()])
        if '$coordinates' in line:
            coords = []
            atom_symbols = []
            line = next(ss)
            while '$end' not in line:
                lsp = line.split()
                atom_symbols.append(lsp[0])
                coords.append([float(x) / Units.bohr for x in lsp[1:]])
                line = next(ss)
    info_dict = {
        'int_number': 1,
        'lattice': matrix,
        'coords': coords,
        'atom_symbols': atom_symbols}
    return info_dict


def writer(info_dict) -> str:
    ss = []
    if 'lattice' in info_dict:
        ss += ['$cell vectors']
        for vec in info_dict['lattice']:
            v = [x * Units.bohr for x in vec]
            ss.append('{:30.14f}{:20.14f}{:20.14f}'.format(*v))
    ss += ['$coordinates']
    for at, coord in zip(info_dict['atom_symbols'], info_dict['coords']):
        v = [x * Units.bohr for x in coord]
        ss.append('{:10}{:20.14f}{:20.14f}{:20.14f}'.format(at, *v))
    ss += ['$end\n']
    string = '\n'.join(ss)
    return string
