# -*- coding: utf-8 -*-
import wx
import numpy as np
from itertools import product
from .components import MyID
from .layout import (MyMolCompareDialog,
                     MyFileConfigDialog,
                     MySymmetryDialog,
                     MyLatticeDialog,
                     MyRotateDialog,
                     MyTwistDialog,
                     MyBondDialog)
from sdm.models import SpaceGroup, Lattice
from . import logger as root_logger
logger = root_logger.getChild('dialog')


class IMyMolCompareDialog(MyMolCompareDialog):
    def __init__(self, parent, namedict=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.main = parent
        self.choice_map = choice_map = {}
        for i, (k, name) in enumerate(namedict.items()):
            self.check_list_box_1.Append(name)
            self.check_list_box_2.Append(name)
            choice_map[i] = k

    def onBtnOK(self, event):
        choice_ref = [self.choice_map[i] for i in self.check_list_box_1.GetCheckedItems()]
        choice_com = [self.choice_map[i] for i in self.check_list_box_2.GetCheckedItems()]
        without_H = self.checkbox_1.IsChecked()
        logger.info("ref={}, com={}".format(choice_ref, choice_com))
        for ref, com in product(choice_ref, choice_com):
            if ref != com:
                wx.CallAfter(self.main.canvas.structCompare, ref=ref, com=com, without_H=without_H)
        self.Destroy()

    def onBtnCancel(self, event):
        self.Destroy()


class IMyFileConfigDialog(MyFileConfigDialog):
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)
        self.text = [self.text_ctrl_1, self.text_ctrl_2]
        self.attrs = ["title", "symprec"]
        self.types = [str, float]
        self.struct = None
        self.load()

    def load(self, struct=None):
        self.struct = struct = self.Parent.canvas.vs.data if struct is None else struct
        for attr, obj in zip(self.attrs, self.text):
            if hasattr(struct, attr):
                obj.SetValue(str(getattr(struct, attr)))

    def apply(self):
        attr_dict = {attr: t(obj.GetValue()) for attr, t, obj in zip(self.attrs, self.types, self.text)}
        wx.CallAfter(self.Parent.canvas.structSetAttr, attr_dict=attr_dict)

    def OnOK(self, event):
        self.apply()
        self.Close()

    def OnCancel(self, event):
        self.Close()


class IMySymmetryDialog(MySymmetryDialog):
    def __init__(self, parent, **kwargs):
        MySymmetryDialog.__init__(self, parent, **kwargs)
        self.text = [
            self.text_ctrl_1,
            self.text_ctrl_2,
            self.text_ctrl_3,
            self.text_ctrl_4,
            self.text_ctrl_5
        ]
        self.spg = None
        self.struct = None
        self.load()

    def update(self, spg=None):
        if spg is None:
            spg = self.spg
        else:
            self.spg = spg
        data = [spg.int_symbol, spg.int_number, spg.hall_symbol, spg.hall_number, '\n'.join(spg.xyz_string)]
        for val, obj in zip(data, self.text):
            if val is not None:
                obj.SetValue(str(val))

    def load(self, struct=None):
        self.struct = struct = self.Parent.canvas.vs.data if struct is None else struct
        if hasattr(struct, 'space_group'):
            self.update(struct.space_group.copy())

    def apply(self):
        attr_dict = {"space_group": self.spg}
        wx.CallAfter(self.Parent.canvas.structSetAttr, attr_dict=attr_dict)

    def OnUpperTextEnter(self, event):  # wxGlade: MySymmetryDialog.<event_handler>
        widget = event.GetEventObject()
        wid = widget.GetId()
        val = widget.GetValue()
        if wid == MyID.DlgSymmTxt1:
            spg = SpaceGroup.fromIntSymbol(val)
        elif wid == MyID.DlgSymmTxt2:
            spg = SpaceGroup.fromIntNumber(int(val))
        elif wid == MyID.DlgSymmTxt3:
            spg = SpaceGroup.fromSymbol(val)
        elif wid == MyID.DlgSymmTxt4:
            spg = SpaceGroup.fromHallNumber(int(val))
        self.update(spg)

    def OnOK(self, event):
        self.apply()
        self.Close()

    def OnApply(self, event):
        self.apply()

    def OnCancel(self, event):
        self.Close()

    def OnRefresh(self, event):
        self.load()


class IMyLatticeDialog(MyLatticeDialog):
    def __init__(self, parent, **kwargs):
        MyLatticeDialog.__init__(self, parent, **kwargs)
        self.text = [
            self.text_ctrl_1,
            self.text_ctrl_2,
            self.text_ctrl_3,
            self.text_ctrl_4,
            self.text_ctrl_5,
            self.text_ctrl_6,
        ]
        self.lattice = None
        self.struct = None
        self.load()

    def load(self, struct=None):
        self.struct = struct = self.Parent.canvas.vs.data if struct is None else struct
        if hasattr(struct, 'lattice'):
            self.update(struct.lattice.copy())

    def update(self, lattice=None, lower=True, upper=True):
        if lattice is None:
            lattice = self.lattice
        else:
            self.lattice = lattice
        if upper:
            for val, ctrl in zip(lattice.parameters, self.text):
                ctrl.SetValue('{:.6f}'.format(val))
        if lower:
            mat_str = '\n'.join(('{:15.8f}'*3).format(*l) for l in lattice.matrix)
            self.text_ctrl_matrix.SetValue(mat_str)

    def apply(self):
        mat_str = self.text_ctrl_matrix.GetValue()
        matrix = np.array(mat_str.split(), dtype='f').reshape((3, 3))
        args = [matrix,
                self.checkbox_1.GetValue(),
                self.checkbox_2.GetValue(),
                self.checkbox_3.GetValue()]
        wx.CallAfter(self.Parent.canvas.sessionModify, attrib="setLattice", args=args, no_return=True)

    def OnUpperTextEnter(self, event):
        latt6 = [float(ctrl.GetValue()) for ctrl in self.text]
        lattice = Lattice.fromParameters(*latt6)
        self.update(lattice, upper=False)

    def OnOK(self, event):
        self.apply()
        self.Close()

    def OnCancel(self, event):
        self.Close()


class IMyBondDialog(MyBondDialog):
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)
        self.main = parent
        self.selection = parent.canvas.vs.selection[-2:]
        if len(self.selection) != 2:
            self.Destroy()
        struct = parent.canvas.vs.data
        length = struct.getDistance(*self.selection)
        bond = struct.graph.get_edge_data(*self.selection)
        self.follower = struct.GUIgetFollowerGroups(self.selection)

        self.config_d = {0: "TG", 1: "TA", 2: "FIX"}
        self.config = ["TG", "TG"]
        self.vrange = vrange = (0.01, 5.0)
        srange = self.slider_1.GetRange()
        self.slider_tickrate = (srange[1] - srange[0]) / (vrange[1] - vrange[0])

        bond_type = 0 if bond is None else bond["type"].value
        self.combo_box_1.Select(bond_type)
        self.text_ctrl_1.SetValue("{:.4f}".format(length))
        self.setSlider(length)

        wx.CallAfter(self.main.canvas.sessionChange, struct)

    def setSlider(self, value):
        """关联参数: self.slider_tickrate
        """
        val = int(round((value + self.vrange[0]) * self.slider_tickrate))
        self.slider_1.SetValue(val)

    def getSlider(self):
        """关联参数: self.slider_tickrate
        """
        return self.slider_1.GetValue() / self.slider_tickrate + self.vrange[0]

    def applyModify(self):
        value = float(self.text_ctrl_1.GetValue())
        attrib = "GUItranslate"
        args = [self.selection, value, self.follower, self.config]
        # print(args)
        wx.CallAfter(self.main.canvas.sessionModify, attrib=attrib, args=args, no_return=True, preview=True)

    def OnBondTypeSelect(self, event):
        bond_type_new = self.combo_box_1.GetSelection()
        if bond_type_new == 0:
            attrib = "delBonds"
            args = [[list(self.selection)]]
        else:
            attrib = "setBonds"
            args = [[list(self.selection) + [bond_type_new]]]
        wx.CallAfter(self.main.canvas.sessionModify, attrib=attrib, args=args, no_return=True, preview=True)

    def OnConfigSelect(self, event):
        self.config = [
            self.config_d[self.combo_box_2.GetSelection()],
            self.config_d[self.combo_box_3.GetSelection()],
        ]

    def OnScroll(self, event):
        value = self.getSlider()
        self.text_ctrl_1.SetValue("{:.4f}".format(value))
        self.applyModify()

    def OnSliderEnter(self, event):
        value = float(self.text_ctrl_1.GetValue())
        self.slider_1.SetValue(value)
        self.applyModify()

    def OnOK(self, event):
        self.Destroy()

    def OnCancel(self, event):
        wx.CallAfter(self.main.canvas.sessionChange, undo=True)
        self.Destroy()


class IMyRotateDialog(MyRotateDialog):
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)
        self.main = parent
        self.selection = parent.canvas.vs.selection[-3:]
        if len(self.selection) != 3:
            self.Destroy()
        struct = parent.canvas.vs.data
        angle = struct.getAngle(*self.selection)
        self.follower = struct.GUIgetFollowerGroups(self.selection)
        # [[RG|TG|RA|FIX], [TG|TA|FIX], [RG|TG|RA|FIX]]
        self.config_d13 = {0: "RG", 1: "TG", 2: "RA", 3: "FIX"}
        self.config_d2 = {0: "TG", 1: "TA", 2: "FIX"}
        self.OnConfigSelect(None)
        self.vrange = vrange = (0.0, 180)
        srange = self.slider_1.GetRange()
        self.slider_tickrate = (srange[1] - srange[0]) / (vrange[1] - vrange[0])

        self.text_ctrl_1.SetValue("{:.4f}".format(angle))
        self.setSlider(angle)

        wx.CallAfter(self.main.canvas.sessionChange, struct)

    def setSlider(self, value):
        """关联参数: self.slider_tickrate
        """
        val = int(round((value + self.vrange[0]) * self.slider_tickrate))
        self.slider_1.SetValue(val)

    def getSlider(self):
        """关联参数: self.slider_tickrate
        """
        return self.slider_1.GetValue() / self.slider_tickrate + self.vrange[0]

    def applyModify(self):
        value = float(self.text_ctrl_1.GetValue())
        attrib = "GUIrotate"
        args = [self.selection, value, self.follower, self.config]
        # print(args)
        wx.CallAfter(self.main.canvas.sessionModify, attrib=attrib, args=args, no_return=True, preview=True)

    def OnConfigSelect(self, event):
        self.config = [
            self.config_d13[self.combo_box_2.GetSelection()],
            self.config_d2[self.combo_box_3.GetSelection()],
            self.config_d13[self.combo_box_4.GetSelection()],
        ]

    def OnScroll(self, event):
        value = self.getSlider()
        self.text_ctrl_1.SetValue("{:.4f}".format(value))
        self.applyModify()

    def OnSliderEnter(self, event):
        value = float(self.text_ctrl_1.GetValue())
        self.slider_1.SetValue(value)
        self.applyModify()

    def OnOK(self, event):
        self.Destroy()

    def OnCancel(self, event):
        wx.CallAfter(self.main.canvas.sessionChange, undo=True)
        self.Destroy()


class IMyTwistDialog(MyTwistDialog):
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)
        self.main = parent
        self.selection = parent.canvas.vs.selection[-4:]
        if len(self.selection) != 4:
            self.Destroy()
        struct = parent.canvas.vs.data
        angle = struct.getDihedral(*self.selection)
        self.follower = struct.GUIgetFollowerGroups(self.selection)
        # [[RS|RG|RA|FIX], [RS|RG|RA|FIX]]
        self.config_d = {0: "RS", 1: "RG", 2: "RA", 3: "FIX"}
        self.OnConfigSelect(None)
        self.vrange = vrange = (-180.0, 180.0)
        srange = self.slider_1.GetRange()
        self.slider_tickrate = (srange[1] - srange[0]) / (vrange[1] - vrange[0])

        self.text_ctrl_1.SetValue("{:.4f}".format(angle))
        self.setSlider(angle)

        wx.CallAfter(self.main.canvas.sessionChange, struct)

    def setSlider(self, value):
        """关联参数: self.slider_tickrate
        """
        val = int(round((value - self.vrange[0]) * self.slider_tickrate))
        self.slider_1.SetValue(val)

    def getSlider(self):
        """关联参数: self.slider_tickrate
        """
        return self.slider_1.GetValue() / self.slider_tickrate + self.vrange[0]

    def applyModify(self):
        value = float(self.text_ctrl_1.GetValue())
        attrib = "GUItwist"
        args = [self.selection, value, self.follower, self.config]
        # print(args)
        wx.CallAfter(self.main.canvas.sessionModify, attrib=attrib, args=args, no_return=True, preview=True)

    def OnConfigSelect(self, event):
        self.config = [
            self.config_d[self.combo_box_2.GetSelection()],
            self.config_d[self.combo_box_3.GetSelection()],
        ]

    def OnScroll(self, event):
        value = self.getSlider()
        self.text_ctrl_1.SetValue("{:.4f}".format(value))
        self.applyModify()

    def OnSliderEnter(self, event):
        value = float(self.text_ctrl_1.GetValue())
        self.setSlider(value)
        self.applyModify()

    def OnOK(self, event):
        self.Destroy()

    def OnCancel(self, event):
        wx.CallAfter(self.main.canvas.sessionChange, undo=True)
        self.Destroy()
