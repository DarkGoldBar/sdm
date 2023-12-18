# -*- coding: utf-8 -*-
from __future__ import absolute_import, unicode_literals, division, print_function
import numpy as np
from jinja2 import Template


def parser(string):
    struct_lines = string.split('\n')
    # line 1: dimensionality, centring, crystal system, and comment
    firstline = struct_lines[0].split(None, 3)
    dim_code, center_code, system_code = (int(x) for x in firstline[:3])
    title = 'CRYSTAL-Structure'
    if len(firstline) > 3:
        title = firstline[3]
    # line 2-4: cartesian components of direct lattice vectors
    matrix = np.array([x.split() for x in struct_lines[1:4]], dtype=float)
    inv_matrix = np.linalg.inv(matrix)
    # line 5: number of symmetry operators
    # line x + 1-3: symmetry operators matrices in cartesian coordinates
    # line x + 4: cartesian components of the translation
    nsymm = int(struct_lines[4].strip())
    nlatom = 5 + 4 * nsymm
    affine_matrix_raw = np.array(' '.join(struct_lines[5: nlatom]).split(), dtype=float).reshape((-1, 4, 3)).round(12)
    affine_matrix = []
    for mat in affine_matrix_raw:
        affine = np.eye(4)
        affine[:3, :3] = mat[:3, :3]
        affine[:3, 3] = np.dot(mat[3, :3], inv_matrix)
        affine_matrix.append(affine)
    # line n: number of atoms in the primitive cell (up to CRYSTAL14, irreducible atoms only)
    natom = int(struct_lines[nlatom].strip())
    atom_coords_str = ' '.join(struct_lines[nlatom + 1: nlatom + natom + 1])
    atom_coords = np.array(atom_coords_str.split(), dtype=float).reshape((-1, 4))
    return {
        "title": title,
        "atom_numbers": [int(round(x)) for x in atom_coords[:, 0]],
        "lattice": matrix,
        "coords": atom_coords[:, 1:],
        "affine_matrices": affine_matrix,
        "dim_code": dim_code,
        "lattice_centering": centering_code_map.get(center_code, 'P'),
        "system_code": system_code,
        "cartesian_symmetry_ops": True,
    }


def writer(info_dict):
    if "lattice" in info_dict:
        dim_code = 3
        matrix = info_dict['lattice']
        if np.all(np.isclose(matrix[2] == [0, 0, 500])):
            dim_code = 2
            if np.all(np.isclose(matrix[1] == [0, 500, 0])):
                dim_code = 1
        center_code = centering_code_map[info_dict['int_symbol'][0]]
        system_code = crystal_system_list.index(info_dict['lattice_system'])
        rot_trans = info_dict['CRYSTAL_cartops']
        render_data = dict(
            title=info_dict['title'],
            dim_code=dim_code,
            center_code=center_code,
            system_code=system_code,
            matrix=matrix,
            rot_trans=rot_trans,
            atoms_coords=list(zip(info_dict['atom_numbers'], info_dict['coords'])),
        )
    else:
        dim_code = 0
        matrix = np.eye(3) * 500
        render_data = dict(
            title=info_dict['title'],
            dim_code=dim_code,
            center_code=1,
            system_code=1,
            matrix=matrix,
            rot_trans=np.vstack(np.eye(3), [[0, 0, 0]]),
            atoms_coords=list(zip(info_dict['atom_numbers'], info_dict['coords'])),
        )
    template = Template(TEMPLATE)
    return template.render(**render_data)


TEMPLATE = """{{dim_code}} {{center_code}} {{system_code}} {{title}}
{%- for vec in matrix %}
{% for v in vec %}{{ "%20.12e " % v }}{% endfor %}
{%- endfor %}
{{rot_trans|length // 4}}
{%- for vec in rot_trans %}
{% for v in vec %}{{ "%20.12e " % v }}{% endfor %}
{%- endfor %}
{{atoms_coords|length}}
{%- for atom, coord in atoms_coords %}
{{"%4d " % atom.z}}{% for c in coord %}{{ "%20.12e " % c}}{% endfor %}
{%- endfor %}
"""

centering_code_map = {
    'P': 1,
    'R': 1,
    'A': 2,
    'B': 3,
    'C': 4,
    'F': 5,
    'I': 6,
    'H': 7,
    1: 'P',
    2: 'A',
    3: 'B',
    4: 'C',
    5: 'F',
    6: 'I',
    7: 'H',
}

crystal_system_list = [  # 编号从1到7
    'dummy'
    'Triclinic',
    'Monoclinic',
    'Orthorhombic',
    'Tetragonal',
    'Trigonal',
    'Hexagonal',
    'Cubic'
]
