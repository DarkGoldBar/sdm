# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2020/9/23
from jinja2 import Template
import re
import numpy as np


def parser(string):
    string = re.sub(r'#.*', '', string.strip())
    if len(string) == 0:
        raise ValueError("Empty POSCAR")
    lines = re.split(r'\s*\n\s*', string)
    # Line 1 注释（或名称）
    title = lines[0]
    # Line 2 缩放因子
    scale = float(lines[1])
    # Line 3-5 晶胞矩阵
    lattice = np.array([line.split() for line in lines[2:5]], dtype=float)
    if scale < 0:
        # In vasp, a negative scale factor is treated as a volume. We need
        # to translate this to a proper lattice vector scaling.
        vol = abs(np.linalg.det(lattice))
        lattice *= (-scale / vol) ** (1 / 3)
    else:
        lattice *= scale

    # Line 6-7 元素类型(及数量)
    """
    Atoms and number of atoms in POSCAR written with vasp appear on
    multiple lines when atoms of the same type are not grouped together
    and more than 20 groups are then defined ...
    Example :
    Cr16 Fe35 Ni2
        1.00000000000000
            8.5415010000000002   -0.0077670000000000   -0.0007960000000000
        -0.0077730000000000    8.5224019999999996    0.0105580000000000
        -0.0007970000000000    0.0105720000000000    8.5356889999999996
        Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Ni   Fe   Cr   Fe   Cr
        Fe   Ni   Fe   Cr   Fe
            1   1   2   4   2   1   1   1     2     1     1     1     4     1     1     1     5     3     6     1
            2   1   3   2   5
    Direct
        ...
    """
    ipos = 5
    atom_types = []
    atom_qty = []
    while True:
        lsp = lines[ipos].split()
        if len(atom_qty) == 0 and all(x.isalpha() for x in lsp):
            atom_types += lsp
        elif all(x.isdigit() for x in lsp):
            atom_qty += [int(x) for x in lsp]
        elif len(atom_types) != len(atom_qty):
            raise IOError("POSCAR: number of atom Type and Quantity not match @ line {}".format(ipos))
        elif ipos > 20:
            raise IOError("POSCAR: loop too long @ line {}".format(ipos))
        else:
            break
        ipos += 1
    atomic_symbols = [x for tp, qty in zip(atom_types, atom_qty) for x in [tp]*qty]

    # Line 7(or 8) 坐标类型
    postype = lines[ipos].split()[0]
    sdynamics = False
    ipos += 1
    # Selective dynamics
    if postype[0] in "sS":
        sdynamics = True
        ipos += 1
        postype = lines[ipos].split()[0]

    cart = postype[0] in "cCkK"
    # print(f"{postype=}, {cart=}")
    nsites = sum(atom_qty)

    linesp = list(zip(*(ll.split() for ll in lines[ipos: ipos + nsites])))
    atom_titles = None
    atitle_idx = 6 if sdynamics else 3
    if len(linesp) > atitle_idx:
        atom_titles = linesp[atitle_idx]
    coords = np.array(linesp[:3], float).T
    if cart:
        coords = coords * scale
    else:
        coords = np.dot(coords, lattice)
    if sdynamics:
        selective_dynamics = np.array([[t.upper()[0] == "T" for t in tok] for tok in linesp[3:6]]).T
    else:
        selective_dynamics = None

    info_dict = {"hall_number": 1,
                 "lattice": lattice,
                 "atom_symbols": atomic_symbols,
                 "atom_titles": atom_titles,
                 "title": title,
                 "coords": coords,
                 "selective_dynamics": selective_dynamics}
    return info_dict


def __concentrate_symbol(atom_symbols):
    symbol_list = [atom_symbols[0]]
    count_list = [1]
    for c in atom_symbols[1:]:
        if c == symbol_list[-1]:
            count_list[-1] += 1
        else:
            symbol_list.append(c)
            count_list.append(1)
    return symbol_list, count_list


def writer(info_dict, direct=True, scaling_factor=1.0, selective_dynamics=None):
    """
    scaling_factor: 比例因子,用于缩放所有晶格矢量和所有原子坐标（此值为负时, 表示晶格总体积), 默认值为1.0
    selective_dynamics: 笛卡尔坐标, 并指定每个原子的坐标是否可以改变,
    """
    template = Template(POSCAR_TEMPLEATE)
    info_dict['lattice'] = [tuple(line) for line in info_dict['lattice']]
    if selective_dynamics is not None:
        selective_dynamics = [tuple(x) for x in selective_dynamics]
    symbol_list, count_list = __concentrate_symbol(info_dict['atom_symbols'])
    poscar_str = template.render(
        direct=direct,
        count_list=count_list,
        symbol_list=symbol_list,
        scaling_factor=scaling_factor,
        selective_dynamics=selective_dynamics,
        zip=zip,
        **info_dict
    )
    return poscar_str


POSCAR_TEMPLEATE = """{{title}}
{{scaling_factor}}
{%- for row in lattice %}
{{ "%12.6f%12.6f%12.6f" % row }}
{%- endfor %}
{% for symbol in symbol_list -%} {{ symbol }} {% endfor %}
{% for count in count_list -%} {{ count }} {% endfor %}
{%- if selective_dynamics %}
Selective dynamics
{%- endif %}
{%- if direct %}
Direct
{%- else %}
Cartesian
{%- endif %}
{% for cd, fcd, syb in zip(coords, frac_coords, atom_symbols) %}
{%- if direct %}
{%- for x in fcd %}{{ "%12.6f" % x }}{% endfor %}
{%- else %}
{%- for x in cd %}{{ "%12.6f" % x }}{% endfor %}
{%- endif -%}
{%- if selective_dynamics %}
{%- for f in selective_dynamics[loop.index-1] %}{% if f %} T{% else %} F{% endif %}{% endfor %}
{%- endif -%}
{{ " %s" % syb }}
{% endfor %}
"""
