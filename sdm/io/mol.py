# -*- coding:utf-8 -*-
# Author : bochen.li
# Data : 2020/9/22
import re
import os
from jinja2 import Template


def parser(mol_str):
    lines = re.split(r'\r?\n', mol_str)
    title = lines[0]  # 第一行为文件名
    atom_count = int(lines[3][0:3])  # 第三行第一项为原子个数
    bond_count = int(lines[3][3:6])  # 第三行第二项为化学键个数
    atom_symbols, coords = [], []
    for l in lines[4:4 + atom_count]:
        ll = l.split()
        atom_symbols.append(ll[3])  # 每行第四项为原子元素符
        if len(ll) == 16 or len(ll) == 5:
            coords.append(list(map(float, [ll[0], ll[1], ll[2]])))
        elif len(ll) == 15:
            x = ll[0]
            y = ll[1][:-10]
            z = ll[1][-10:]
            coords.append(list(map(float, [x, y, z])))
        elif len(ll) == 14:
            x, y, z = map(float, re.findall(r"(-?\d+.\d\d\d\d)\s*(-?\d+.\d\d\d\d)\s*(-?\d+.\d\d\d\d)", ll[0])[0])
            coords.append(list(map(float, [x, y, z])))
        else:
            print(l)
            raise Exception('lenght == %d' % len(ll))

    bonds = []
    for l in lines[4 + atom_count:4 + atom_count + bond_count]:  # 从 4+atom_count 到 4+atom_count+bond_count 行， 每行表示一个化学键
        # 第一，二项为原子序号，第三项为化学键类型
        a_idx, b_idx, bond_type = int(l[:3]) - 1, int(l[3:6]) - 1, int(l[6:9].strip())
        if a_idx < b_idx:
            bonds.append((a_idx, b_idx, bond_type))
        else:
            bonds.append((b_idx, a_idx, bond_type))

    charges = [0] * atom_count
    for l in lines[(4 + atom_count + bond_count):]:
        if len(l) < 6:
            continue
        if l[0] == 'M' and l[3:6] == 'CHG':
            i = 10
            while i < len(l):
                charges[int(l[i: i + 3]) - 1] = float(l[i + 4: i + 7])
                i += 8

    return {"atom_symbols": atom_symbols,
            "bonds": bonds,
            "title": title,
            "coords": coords,
            "charge": int(sum(charges)),
            "atom_charges": charges}


def writer(info_dict):
    info_dict['coords'] = [tuple(coord) for coord in info_dict['coords']]
    info_dict['bonds'] = [tuple((i+1, j+1, t, 0)) for i, j, t in info_dict.get('bonds', [])]
    # M..CHG..n..i..c  # n <= 8, 后接n组电荷信息
    atom_charges = info_dict.get('atom_charges', [])
    m_charge = [(idx+1, v) for idx, v in enumerate(atom_charges) if v != 0]
    m_charge = [m_charge[i*8:i*8+8] for i in range((len(m_charge)-1)//8+1)]
    info_dict['m_charge'] = m_charge
    template = Template(MOL_TEMPLATE)
    mol_str = template.render(zip=zip, **info_dict).lstrip()
    return mol_str


MOL_TEMPLATE = """{{title}}
SDM 3D

{{"%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (coords|length, bonds|length)}}
{%- for atom_symbol, coord in zip(atom_symbols, coords) %}
{{"%10.4f%10.4f%10.4f" % coord}}\
{{" %-2s  0  0  0  0  0  0  0  0  0  0  0  0" % atom_symbol}}
{%- endfor %}
{%- for bond in bonds %}
{{"%3d%3d%3d%3d" % bond}}
{%- endfor %}
{%- for line in m_charge %}
M  CHG{{"%3d" % line|length}}
{%- for x in line %}{{"%3d%3d" % x}}{% endfor %}
{%- endfor %}
M  END
"""


def rdkit_parse_mol(file_path):
    """
    解析.mol文件， 返回rdkit的mol对象
    Args:
        file_path:  .mol文件路径

    Returns:  rdkit.Chem.rdchem.Mol (参考https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol)

    """
    from rdkit import Chem
    if not os.path.exists(file_path):
        raise Exception("no such file {}".format(file_path))
    return Chem.MolFromMolFile(file_path, sanitize=False, removeHs=False)


def rdkit_write_mol(mol, out_file_path):
    """
    输入rdkit的mol对象， 生成.mol文件
    Args:
        mol: rdkit.Chem.rdchem.Mol (参考https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol)
        out_file_path: 生成的.mol文件路径

    """
    from rdkit import Chem
    Chem.MolToMolFile(mol, out_file_path)
