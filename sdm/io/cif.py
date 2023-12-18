# -*- coding:utf-8 -*-
import re
import warnings
from collections import OrderedDict, deque
import numpy as np
from jinja2 import Template


def _str2float(text):
    """字符串转换为符点型"""
    return float(re.sub(r"\(.+\)", "", text))


def _remove_non_ascii(s):
    """
    Remove non-ascii characters in a file. Needed when support for non-ASCII
    is not available.

    Args:
        s (str): Input string

    Returns:
        String with all non-ascii characters removed.
    """
    return "".join(i for i in s if ord(i) < 128)


class CifBlock:
    """
    Object for storing cif data. All data is stored in a single dictionary.
    Data inside loops are stored in lists in the data dictionary, and
    information on which keys are grouped together are stored in the loops
    attribute.
    """

    maxlen = 70  # not quite 80 so we can deal with semicolons and things

    def __init__(self, data, loops, header):
        """
        Args:
            data: dict or OrderedDict of data to go into the cif. Values should
                    be convertible to string, or lists of these if the key is
                    in a loop
            loops: list of lists of keys, grouped by which loop they should
                    appear in
            header: name of the block (appears after the data_ on the first
                line)
        """
        self.loops = loops
        self.data = data
        # AJ says: CIF Block names cannot be more than 75 characters or you
        # get an Exception
        self.header = header[:74]

    def __eq__(self, other):
        return self.loops == other.loops and self.data == other.data and self.header == other.header

    def __getitem__(self, key):
        return self.data[key]

    def __str__(self):
        """
        Returns the cif string for the data block
        """
        s = ["data_{}".format(self.header)]
        keys = self.data.keys()
        written = []
        for k in keys:
            if k in written:
                continue
            for l in self.loops:
                # search for a corresponding loop
                if k in l:
                    s.append(self._loop_to_string(l))
                    written.extend(l)
                    break
            if k not in written:
                # k didn't belong to a loop
                v = self._format_field(self.data[k])
                if len(k) + len(v) + 3 < self.maxlen:
                    s.append("{}   {}".format(k, v))
                else:
                    s.extend([k, v])
        return "\n".join(s)

    def _loop_to_string(self, loop):
        s = "loop_"
        for l in loop:
            s += "\n " + l
        for fields in zip(*[self.data[k] for k in loop]):
            line = "\n"
            for val in map(self._format_field, fields):
                if val[0] == ";":
                    s += line + "\n" + val
                    line = "\n"
                elif len(line) + len(val) + 2 < self.maxlen:
                    line += "  " + val
                else:
                    s += line
                    line = "\n  " + val
            s += line
        return s

    def _format_field(self, v):
        v = v.__str__().strip()
        # 移除依赖 textwrap
        # if len(v) > self.maxlen:
        #     return ";\n" + textwrap.fill(v, self.maxlen) + "\n;"
        # add quotes if necessary
        if v == "":
            return '""'
        if (" " in v or v[0] == "_") and not (v[0] == "'" and v[-1] == "'") and not (v[0] == '"' and v[-1] == '"'):
            if "'" in v:
                q = '"'
            else:
                q = "'"
            v = q + v + q
        return v

    @classmethod
    def _process_string(cls, string):
        # remove comments
        string = re.sub(r"(\s|^)#.*$", "", string, flags=re.MULTILINE)
        # remove empty lines
        string = re.sub(r"^\s*\n", "", string, flags=re.MULTILINE)
        # remove non_ascii
        string = _remove_non_ascii(string)
        # since line breaks in .cif files are mostly meaningless,
        # break up into a stream of tokens to parse, rejoining multiline
        # strings (between semicolons)
        q = deque()
        multiline = False
        ml = []
        # this regex splits on spaces, except when in quotes.
        # starting quotes must not be preceded by non-whitespace
        # (these get eaten by the first expression)
        # ending quotes must not be followed by non-whitespace
        p = re.compile(r"""(data_.*|(?:[^'"\s][\S]*))|'(.*?)'(?!\S)|"(.*?)"(?!\S)""")
        for l in string.splitlines():
            if multiline:
                if l.startswith(";"):
                    multiline = False
                    q.append(("", "", "", " ".join(ml)))
                    ml = []
                    l = l[1:].strip()
                else:
                    ml.append(l)
                    continue
            if l.startswith(";"):
                multiline = True
                ml.append(l[1:].strip())
            else:
                for s in p.findall(l):
                    # s is tuple. location of the data in the tuple
                    # depends on whether it was quoted in the input
                    q.append(s)
        return q

    @classmethod
    def from_string(cls, string):
        """
        Reads CifBlock from string.

        :param string: String representation.
        :return: CifBlock
        """
        q = cls._process_string(string)
        header = q.popleft()[0][5:].strip("'").strip('"')
        data = OrderedDict()
        loops = []
        while q:
            s = q.popleft()
            # cif keys aren't in quotes, so show up in s[0]
            if s[0] == "_eof":
                break
            if s[0].startswith("_"):
                try:
                    data[s[0]] = "".join(q.popleft())
                except IndexError:
                    data[s[0]] = ""
            elif s[0].startswith("loop_"):
                columns = []
                items = []
                while q:
                    s = q[0]
                    if s[0].startswith("loop_") or not s[0].startswith("_"):
                        break
                    columns.append("".join(q.popleft()))
                    data[columns[-1]] = []
                while q:
                    s = q[0]
                    if s[0].startswith("loop_") or s[0].startswith("_"):
                        break
                    items.append("".join(q.popleft()))
                n = len(items) // len(columns)
                assert len(items) % n == 0
                loops.append(columns)
                for k, v in zip(columns * n, items):
                    data[k].append(v.strip())
            elif "".join(s).strip() != "":
                warnings.warn("Possible issue in cif file" " at line: {}".format("".join(s).strip()))
        return cls(data, loops, header)


def parser(cif_string):
    data = CifBlock.from_string(cif_string)
    title = data.header
    atom_symbols = [atom_symbol for atom_symbol in data["_atom_site_type_symbol"]]
    atom_titles = data.data.get("_atom_site_label")
    latt6 = [float(re.sub(r"\([0-9.]+\)", "", data[key])) \
             for key in ["_cell_length_a", "_cell_length_b", "_cell_length_c", \
                         "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]]
    frac_coords = np.array([[_str2float(i) for i in data["_atom_site_fract_x"]],
                            [_str2float(i) for i in data["_atom_site_fract_y"]],
                            [_str2float(i) for i in data["_atom_site_fract_z"]]]).T

    hall_number = int_number = symmops_xyz_string = None
    for symmetry_label in ["_symmetry_space_group_name_hall",
                           "_symmetry_space_group_name_hall_"]:
        if data.data.get(symmetry_label):
            hall_number = data.data.get(symmetry_label)
    for symmetry_label in ["_space_group_IT_number",
                           "_space_group_IT_number_",
                           "_symmetry_Int_Tables_number",
                           "_symmetry_Int_Tables_number_"]:
        if data.data.get(symmetry_label):
            int_number = data.data.get(symmetry_label)
    for symmetry_label in ["_symmetry_equiv_pos_as_xyz",
                           "_symmetry_equiv_pos_as_xyz_",
                           "_space_group_symop_operation_xyz",
                           "_space_group_symop_operation_xyz_"]:
        if data.data.get(symmetry_label):
            symmops_xyz = data.data.get(symmetry_label)

    return {
        "title": title,
        "atom_symbols": atom_symbols,
        "atom_titles": atom_titles,
        "latt6": latt6,
        "frac_coords": frac_coords,
        "hall_number": hall_number,
        "int_number": int_number,
        "symmops_xyz": symmops_xyz,
    }


def writer(cif_dict):
    """
    {
        "title",
        "spg_int_symbol",
        "spg_int_number",
        "latt6",
        "symmops_xyz_string",
        "atom_titles",
        "atom_symbols",
        "frac_coords",
    }
    """
    template = Template(CIF_TEMPLATE)
    cif_dict['frac_coords'] = [tuple(coord) for coord in cif_dict['frac_coords']]
    if ' ' in cif_dict['title']:
        cif_dict['title'] = "'" + cif_dict['title'] + "'"
    cif_string = template.render(zip=zip, **cif_dict)
    return cif_string


CIF_TEMPLATE = """data_{{ title }}

_symmetry_space_group_name_H-M '{{ spg_int_symbol }}'
_symmetry_Int_Tables_number {{ spg_int_number }}

_cell_length_a     {{ "%14.8f" % latt6[0] }}
_cell_length_b     {{ "%14.8f" % latt6[1] }}
_cell_length_c     {{ "%14.8f" % latt6[2] }}
_cell_angle_alpha  {{ "%14.8f" % latt6[3] }}
_cell_angle_beta   {{ "%14.8f" % latt6[4] }}
_cell_angle_gamma  {{ "%14.8f" % latt6[5] }}

loop_
_space_group_symop_operation_xyz
{%- for op_string in symmops_xyz_string %}
{{ op_string.replace(" ","") }}
{%- endfor %}

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
{%- for label, atom_symbol, coord in zip(atom_titles, atom_symbols, frac_coords) %}
{{"%5s" % label}}  {{"%-2s" % atom_symbol}}{{"%14.8f %14.8f %14.8f" % coord}}  1.0000
{%- endfor %}
"""
