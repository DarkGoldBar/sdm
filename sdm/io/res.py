# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2020/12/30
from jinja2 import Template
import re
import numpy as np


LATTLIST = 'xPIHFABC'


def parser(string):
    info_dict = {}
    info_dict['symmops_xyz_string'] = ['x,y,z']
    fcoords = []
    atom_titles = []
    atom_symbols = []
    atom_spacies = []
    for line in string.strip().splitlines():
        key_ma = re.match(r'([A-Z]{4}) ', line)
        if key_ma:
            key = key_ma.group(1)
            if key == 'TITL':
                info_dict['title'] = line.split(None, 1)[1]
            elif key == 'CELL':
                info_dict['latt6'] = [float(x) for x in line.split(None)[2:]]
            elif key == 'SYMM':
                info_dict['symmops_xyz_string'].append(line.split(None, 1)[1])
            elif key == 'SFAC':
                atom_spacies = line.split()[1:]
            elif key == 'LATT':
                res_latt = int(line.split()[1])
                if res_latt > 0:
                    info_dict['_symmops_centring'] = LATTLIST[res_latt]
        else:
            lsp = line.split()
            if len(lsp) == 5:
                atom_titles.append(lsp[0])
                atom_symbols.append(atom_spacies[int(lsp[1]) - 1])
                fcoords.append([float(x) for x in lsp[2:]])
            elif re.match(r'END', line):
                break
    info_dict['frac_coords'] = np.array(fcoords)
    info_dict['atom_titles'] = atom_titles
    info_dict['atom_symbols'] = atom_symbols
    return info_dict


def writer(info_dict):
    SYMM = info_dict['symmops_xyz_string']
    SYMM.pop(SYMM.index('x, y, z'))
    SFAC = list(set(info_dict['atom_symbols']))
    if 'H' in SFAC:
        SFAC.remove('H')
        SFAC.insert(0, 'H')
    if 'C' in SFAC:
        SFAC.remove('C')
        SFAC.insert(0, 'C')
    atomlines = zip(
        info_dict['atom_titles'],
        [SFAC.index(x)+1 for x in info_dict['atom_symbols']],
        info_dict['frac_coords'].tolist()
    )
    return Template(TEMPLATE).render(
        title=info_dict['title'],
        latt6=info_dict['latt6'],
        SYMM=SYMM,
        SFAC=SFAC,
        atomlines=atomlines
    )


TEMPLATE = """TITL {{ title }}
CELL 1.0{% for i in latt6 %}{{ "%10.4f" % i }}{% endfor %}
ZERR {{ SYMM|length + 1 }}  0.000 0.000 0.000 0.000 0.000 0.000
LATT -1
{%- for xyz_string in SYMM %}
SYMM {{ xyz_string.upper() }}
{%- endfor %}
SFAC {% for i in SFAC %} {{ i }}{% endfor %}
{%- for att, aid, acd in atomlines %}
{{ "%-5s %-2d" % (att, aid) }} {{ "%12.8f %12.8f %12.8f" % (acd[0], acd[1], acd[2]) }}
{%- endfor %}
END
"""
