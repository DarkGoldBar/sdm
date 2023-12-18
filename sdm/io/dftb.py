# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2020/9/23
import re
import numpy as np
from jinja2 import Template


def parser(string):
    title = 'Untitled'
    if string.startswith('#'):
        title = string.split('\n')[0][1:].strip()

    string = re.sub(r'#.*', '', string).strip()
    lines = re.split(r'\s*\n\s*', string)
    natom, isFC = lines[0].split()
    atom_type_map = lines[1].split()
    isFC = (isFC.upper() == 'F')
    natom = int(natom)
    atom_type_map = {str(k+1): v for k, v in enumerate(atom_type_map)}
    # getAtomCoords
    symbols = []
    coords = []
    for line in lines[2: 2 + natom]:
        sp = line.split()
        ele = atom_type_map.get(sp[1])
        symbols.append(ele)
        coords.append(np.array(sp[2:5], dtype=np.float))
    coords = np.array(coords)

    matrix = None
    matstr = [ll.split() for ll in lines[-3:]]
    if len(matstr[0]) == 3:
        matrix = np.array(matstr, dtype=np.float)
        if isFC:
            coords = np.dot(coords, matrix)

    info_dict = {"title": title,
                 "atom_symbols": symbols,
                 "coords": coords,
                 "hall_number": 1}
    if matrix is not None:
        info_dict["lattice"] = matrix
    return info_dict


def writer(info_dict, direct=True):
    if "lattice" not in info_dict:
        direct = False
    symbols = sorted(set(info_dict["atom_symbols"]))
    itypes = [1 + symbols.index(syb) for syb in info_dict["atom_symbols"]]
    template = Template(GEN_TEMPLATE)
    string = template.render(symbols=symbols,
                             itypes=itypes,
                             direct=direct,
                             zip=zip,
                             **info_dict)
    return string


GEN_TEMPLATE = """# {{ title }}
{{ coords|length }} {% if direct -%} F {%- else -%} C {%- endif %}
{% for symb in symbols -%}
    {{ "%-3s" % symb }}
{%- endfor %}
{% for itype in itypes -%}
    {{ "%4d%3d" % (loop.index, itype) }}
    {%- if direct -%}
    {% for c in frac_coords[loop.index-1] %}{{ "%14.10f " % c}}{% endfor %}
    {%- else -%}
    {% for c in coords[loop.index-1] %}{{ "%14.10f " % c}}{% endfor %}
    {%- endif %}
{% endfor %}
{%- if lattice is defined -%}
0.0 0.0 0.0
{%- for vec in lattice %}
{% for v in vec %}{{ "%14.10f " % v }}{% endfor %}
{%- endfor %}
{%- endif %}
"""
