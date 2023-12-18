# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2021/7/12
import copy
from enum import Enum
from ..data import PERIODIC_TABLE, PERIODIC_INDEX

__all__ = ["ElementBase", "Element", "CustomizedElement", "PERIODIC_TABLE", "PERIODIC_INDEX"]


class ElementBase():
    def __init__(self, *args):
        self._data = {}

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name

    @property
    def data(self) -> dict:
        return self._data.copy()

    @property
    def Z(self) -> int:
        return int(self._data["Atomic no"])

    @property
    def number(self) -> int:
        return int(self._data["Atomic no"])

    @property
    def symbol(self) -> str:
        return str(self._data["symbol"])

    @property
    def long_name(self) -> str:
        return str(self._data["Name"])

    @property
    def X(self) -> float:
        return float(self._data.get("X", "nan"))

    @property
    def atomic_mass(self) -> float:
        return float(self._data.get("Atomic mass", "nan"))

    @property
    def atomic_radius(self) -> float:
        return float(self._data.get("Atomic radius", "nan"))

    @property
    def covelent_radius(self) -> float:
        return float(self._data.get("Covelent radius", "nan"))

    @property
    def van_der_waals_radius(self) -> float:
        return float(self._data.get("Van der waals radius", "nan"))

    @property
    def common_oxidation_states(self) -> list:
        return list(self._data.get("Common oxidation states", list()))

    @property
    def oxidation_states(self) -> list:
        return list(self._data.get("Oxidation states", list()))

    @property
    def vcolor(self) -> str:
        """GUI默认颜色"""
        return str(self._data.get("Visual color", "000fff"))

    @property
    def vradius(self) -> float:
        """GUI默认半径"""
        return float(self._data.get("Visual radius", 0.15))


class ElementEnum(ElementBase, Enum):
    """元素类, 使用枚举类型

    自定义元素时使用 CustomizedElement 类创建新元素

    判断类型时使用基类 ElementBase 判断是否为元素

    调用方法:
        >>> Element("C")
        <Element.C: 'C'>

        >>> Element(8)
        <Element.O: 'O'>
    """
    def __init__(self, *args):
        super().__init__(*args)
        self._data = {k: v for k, v in PERIODIC_TABLE[self._value_].items() if v is not None}
        self.__doc__ = ElementEnum.__doc__

    @classmethod
    def _missing_(cls, value):
        if value in PERIODIC_INDEX:
            return cls(PERIODIC_INDEX[value])


Element = ElementEnum('Element', {k: k for k in PERIODIC_TABLE})


class CustomizedElement(ElementBase):
    """自定义元素

    调用方法:
        >>> CustomizedElement("D", "H", {"Atomic mass": 2})  # 填写元素符号
        <Customized Element D>

    data = {
        'symbol': None,
        'Atomic mass': None,
        'Atomic no': None,
        'Atomic radius': None,
        'Common oxidation states': None,
        'Name': None,
        'Oxidation states': None,
        'Van der waals radius': None,
        'X': None,
        'Covelent radius': None,
        'Visual color': None,
        'Visual radius': None
    }
    """
    def __init__(self, name, copy_from=None, data=None):
        self.name = name
        self.value = copy_from
        self._data = copy.deepcopy(PERIODIC_TABLE.get(copy_from, {}))

        data["symbol"] = name
        self._data.update(data)

    def __repr__(self):
        return "<Customized Element {}>".format(self.name)
