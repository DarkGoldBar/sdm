# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2020/9/20
import io
import os
import re
import importlib


format_name_dict = {
    'cif': 'cif',
    'mol': 'mol',
    'xyz': 'xyz',
    'res': 'res',
    'xsd': 'xsd',
    'gen': 'dftb',
    'vasp': 'vasp',
    'poscar': 'vasp',
    'contcar': 'vasp',
    '34': 'crystal',
    'd12': 'crystal',
    'cosmo': 'cosmo',
    'incoor': 'dmol',
}


def getWriter(filename):
    _ext = ext = os.path.basename(filename).split(".")[-1].lower()
    if re.match(r'opt[ac]\d+', ext):
        _ext = '34'
    format_name = format_name_dict.get(_ext)
    if format_name is None:
        raise IOError('no support for this file: {}'.format(ext))
    format_lib = importlib.import_module(".." + format_name, __name__)
    return format_lib.writer


def getParser(filename):
    _ext = ext = os.path.basename(filename).split(".")[-1].lower()
    if re.match(r'opt[ac]\d+', ext):
        _ext = '34'
    format_name = format_name_dict.get(_ext)
    if format_name is None:
        raise IOError('no support for this file: {}'.format(ext))
    format_lib = importlib.import_module(".." + format_name, __name__)
    return format_lib.parser


def open(filename, string=None, bond_method='empiric', astype=None):
    from ..models import Structure, Molecule  # 不写在这里的话 pip install 时会报错
    """open crystal or molecule file

    Args:
        filename: path to file
        string: string of data, overrides `filename`
        astype:
            "dict" -> dict
            "mole" -> Molecule
            "struct" -> Structure
    """
    astype = "" if astype is None else astype
    if string is None:
        with io.open(filename, encoding="utf-8") as f:
            string = f.read()
    info_dict = getParser(filename)(string)

    # 字典 Dict
    if astype is dict or astype.lower().startswith("d"):
        return info_dict

    if "lattice" in info_dict or "latt6" in info_dict:
        obj = Structure.fromDict(info_dict, bond_method=bond_method)
    else:
        obj = Molecule.fromDict(info_dict, bond_method=bond_method)

    if astype is Molecule or astype.lower().startswith("m"):
        obj = Molecule(obj.graph, obj.coords, title=obj.title)
    elif astype is Structure or astype.lower().startswith("s"):
        info_dict['latt6'] = (1, 1, 1, 90, 90, 90)
        info_dict['hall_number'] = 1
        obj = Structure.fromDict(obj.toDict())

    return obj
