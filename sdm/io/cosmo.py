# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2021/7/8
import numpy as np


def parser(string):
    lines = string.splitlines()
    max_line = len(lines)
    i = 0
    atom_symbols = []
    coords = []
    matrix = []
    while i < max_line:
        ln = lines[i]
        if ln.startswith("$cell vectors"):
            if "[au]" in ln:
                factor = 0.529177249
            else:
                raise Exception("$cell vectors 单位不是 [au], 修改代码")

            for i in range(i+1, i+4):
                ls = lines[i].split()
                matrix.append([float(x) for x in ls])
            matrix = np.array(matrix) * factor

        if ln.startswith("$coordinates xyz"):
            if "[au]" in ln:
                factor = 0.529177249
            else:
                raise Exception("$cell vectors 单位不是 [au], 修改代码")
            i += 1
            ln = lines[i]
            while not ln.startswith("$end"):
                ls = ln.split()
                atom_symbols.append(ls[0])
                coords.append([float(x) for x in ls[1:4]])
                i += 1
                ln = lines[i]
            coords = np.array(coords) * factor
        i += 1

    return {
        "atom_symbols": atom_symbols,
        "lattice": matrix,
        "coords": coords,
        "symmops_xyz": ['x, y, z'],
    }
