# -*- coding:utf-8 -*-
from __future__ import annotations

commands = [
    "VS.structLoadMulti",
    "VS.structLoad",
    "VS.structSave",
    "VS.close",
    "VS.modify",
    "VS.stack",
    "VS.setAttr",
    "VS.setSelect",
    "VS.setCamera",
    "UI.updateStruct",
    "UI.updateSelect",
]


def get_hint(current: str, limit=5) -> list[str]:
    query = current.split()[-1]
    hits = [s for s in commands if s.startswith(query)]
    return hits[:limit]
