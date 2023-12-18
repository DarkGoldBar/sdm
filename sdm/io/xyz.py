# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2020/9/24
from __future__ import absolute_import, division, print_function
import re
import numpy as np
from jinja2 import Template


def parser(string):
    lines = re.split(r"\s*\n\s*", string.strip())
    natom = int(lines[0].strip())
    title = lines[1]
    linesp = list(zip(*(l.strip().split() for l in lines[2: 2 + natom])))
    coords = np.array(linesp[1:4], float).T
    atom_titles = linesp[0]
    atom_symbols = [re.match(r'^[A-Z][a-z]?', tt).group() for tt in atom_titles]
    return {"title": title,
            "atom_symbols": atom_symbols,
            "atom_titles": atom_titles,
            "coords": coords}


def writer(info_dict, remove_titles=False):
    if remove_titles or info_dict.get("atom_titles") is None:
        info_dict['atom_titles'] = info_dict['atom_symbols']
    info_dict['coords'] = [tuple(x) for x in info_dict['coords']]
    template = Template(XYZ_TEMPLATE)
    string = template.render(zip=zip, **info_dict)
    return string


XYZ_TEMPLATE = """{{ coords|length }}
{{ title }}
{% for atom_title, coord in zip(atom_titles, coords) %}\
{{ "%-5s" % atom_title }}{{ "%14.8f%14.8f%14.8f" % coord }}
{% endfor %}
"""
