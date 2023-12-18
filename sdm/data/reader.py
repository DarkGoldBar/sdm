# -*- coding:utf-8 -*-
import re
import json


def read_periodic_table(fp):
    with open(fp) as f:
        data = json.load(f)
    cols, lines = data[0], data[1:]
    table = {}
    for l in lines:
        table[l[0]] = dict(zip(cols, l))
    # 数据处理
    table_index = {d['Atomic no']: k for k, d in table.items()}
    return table, table_index


def read_space_group_table(fp):
    with open(fp) as f:
        data = json.load(f)
    hm_info = data['hm_info']
    num_info = {int(k): v for k, v in data['num_info'].items()}
    sign_to_hm = {v['sign']: k for k, v in hm_info.items()}
    hall_to_hm = {v['hall_num']: k for k, v in hm_info.items() if 'hall_num' in v}
    return hm_info, num_info, sign_to_hm, hall_to_hm


def read_switch_table(fp):
    with open(fp) as f:
        data = json.load(f)
    hm_switch = {}
    for k, v in data.items():
        k = tuple(k.split())
        vs = re.split(r'[,;] ?', v)
        v = ",".join(vs[i] if vs[i+3] == "0" else vs[i] + "+" + vs[i+3] for i in range(3)).replace("+-", "-")
        hm_switch[k] = v
    return hm_switch


def write_periodic_table(table):
    import pandas as pd
    df = table
    if isinstance(table, dict):
        df = pd.DataFrame(table, 'O').T.sort_values("Atomic no")
    col_width = [0] * (len(df.columns) + 1)
    header = ['symbol'] + df.columns.tolist()
    table = [header] + [[i] + v for i, v in zip(df.index.tolist(), df.values.tolist())]
    lines = [[json.dumps(k) for k in val] for val in table]
    for ln in lines:
        col_width = [max(x, len(k)) for x, k in zip(col_width, ln)]
    formatter = ", ".join("{:%d}"%x for x in col_width)
    formatter = "[" + formatter + "]"
    fstring = ",\n".join(formatter.format(*ln) for ln in lines)
    fstring = "[\n" + fstring + "]"
    return fstring


def write_space_group_table(fp):
    """
    Hall number 重复
    7fb43d887b23b623ebc0e15751991af1 324 322
    a46a67a0e6718a2cd937e6f29649d120 328 326
    8b52526eb48ae3c7a17fb9c04f130b72 332 330
    """
    import spglib
    from sdm.models import SpaceGroup, SymmOp
    with open(fp) as f:
        d = json.load(f)

    SIGN_TO_HM = {}
    for k, v in d['hm_info'].items():
        spg = SpaceGroup.fromXyzOps(v['symops'])
        sign = spg._get_sign()
        v['sign'] = sign
        SIGN_TO_HM[sign] = k

    for i in range(1, 531):
        symmetry_dataset = spglib.get_symmetry_from_database(i)
        spg = SpaceGroup([SymmOp.from_rotation_and_translation(r, t) for r, t in zip(symmetry_dataset['rotations'], symmetry_dataset['translations'])])
        sign = spg._get_sign()
        if sign in SIGN_TO_HM:
            hm = SIGN_TO_HM[sign]
            if 'hall_num' not in d['hm_info'][hm]:
                d['hm_info'][hm]['hall_num'] = i
            else:
                print('hall重复', d['hm_info'][hm]['hall_num'], i)
        else:
            print('hall丢失', i, sign)

    with open(fp, 'w') as f:
        json.dump(d, f)
