# -*- coding:utf-8 -*-
from pathlib import Path
from .reader import read_periodic_table, read_space_group_table, read_switch_table

__all__ = ["PERIODIC_TABLE", "PERIODIC_INDEX", "HM_INFO", "NUM_INFO", "SIGN_TO_HM", "HALL_TO_HM", "HM_SWITCH"]

current_path = Path(__file__).parent

PERIODIC_TABLE, PERIODIC_INDEX = read_periodic_table(current_path / "periodic_table.json")

HM_INFO, NUM_INFO, SIGN_TO_HM, HALL_TO_HM = read_space_group_table(current_path / "sg_info_db.json")

HM_SWITCH = read_switch_table(current_path / "symops_switch_db.json")
