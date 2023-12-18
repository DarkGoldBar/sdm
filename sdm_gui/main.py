# -*- coding: utf-8 -*-
import os
import re
import wx
from .components import MyID
from .console import get_hint
from .layout import MyMainFrame
from .dialog import (IMyMolCompareDialog,
                     IMyFileConfigDialog,
                     IMySymmetryDialog,
                     IMyLatticeDialog,
                     IMyBondDialog,
                     IMyRotateDialog,
                     IMyTwistDialog)
from . import config
from . import logger as root_logger
logger = root_logger.getChild('main')


class IMyMainFrame(MyMainFrame):
    def __init__(self, *args, init_struct=None, **kwargs):
        MyMainFrame.__init__(self, *args, **kwargs)
        self.canvas.main = self

        self.Bind(wx.EVT_MENU, self.onMBEDelete, id=MyID.MEDelete)
        self.Bind(wx.EVT_MENU, self.onMBVPageup, id=MyID.MVPageup)
        self.Bind(wx.EVT_MENU, self.onMBVPagedown, id=MyID.MVPagedown)

        entries = [
            wx.AcceleratorEntry(wx.ACCEL_NORMAL, wx.WXK_DELETE, MyID.MEDelete),
            wx.AcceleratorEntry(wx.ACCEL_NORMAL, wx.WXK_PAGEUP, MyID.MVPageup),
            wx.AcceleratorEntry(wx.ACCEL_NORMAL, wx.WXK_PAGEDOWN, MyID.MVPagedown),
        ]
        self.SetAcceleratorTable(wx.AcceleratorTable(entries))

    def onChangeSession(self, info):
        self.text_ctrl_detail.SetValue(info)

    def onChangeSelect(self, info):
        self.statusbar.SetStatusText(info[0], 0)
        self.statusbar.SetStatusText(info[1], 1)
        self.statusbar.SetStatusText(info[2], 2)

    def onMBFOpen(self, event):
        wildcard = (
            "CIF files (*.cif)|*.cif|"
            "XYZ files (*.xyz)|*.xyz|"
            "RES files (*.res)|*.res|"
            "MOL files (*.mol)|*.mol|"
            "*.*|*.*"
        )
        with wx.FileDialog(self, "Open file", wildcard=wildcard, defaultDir=config.PREFER['last'],
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            filenames = fileDialog.GetPaths()
            config.PREFER['last'] = os.path.dirname(filenames[0])
            wx.CallAfter(self.canvas.structLoadMulti, filenames=filenames)

    def onMBFSave(self, event):
        wildcard = (
            "CIF files (*.cif)|*.cif|"
            "XYZ files (*.xyz)|*.xyz|"
            "RES files (*.res)|*.res|"
            "MOL files (*.mol)|*.mol|"
            "*.*|*.*"
        )
        with wx.FileDialog(self, "Save file", wildcard=wildcard, defaultDir=config.PREFER['last'],
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            filename = fileDialog.GetPath()
            config.PREFER['last'] = os.path.dirname(filename)
            wx.CallAfter(self.canvas.structSave, filename=filename)

    def onMBFClose(self, event):
        wx.CallAfter(self.canvas.sessionClose)

    def onMBFExportVASP(self, event):
        logger.error("Event handler 'onMBFExportVASP' not implemented!")

    def onMBFExportCRYSTAL(self, event):
        logger.error("Event handler 'onMBFExportCRYSTAL' not implemented!")

    def onMBFExportDFTB(self, event):
        logger.error("Event handler 'onMBFExportDFTB' not implemented!")

    def onMBEUndo(self, event):
        wx.CallAfter(self.canvas.sessionChange, undo=True)

    def onMBERedo(self, event):
        wx.CallAfter(self.canvas.sessionChange, redo=True)

    def onMBEReload(self, event):
        wx.CallAfter(self.canvas.sessionChange, reset=True)

    def onMBEFindBondsN(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="findBonds", kwargs={"method": None})

    def onMBEFindBondsE(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="findBonds", kwargs={"method": 'e'})

    def onMBEFindBondsEP(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="findBonds", kwargs={"method": 'ep'})

    def onMBEFindBondsI(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="findBonds", kwargs={"method": 'i'})

    def onMBEConvertP1(self, event):
        wx.CallAfter(self.canvas.sessionModify, attrib="getP1Structure")

    def onMBEConvertPrim(self, event):
        wx.CallAfter(self.canvas.sessionModify, attrib="getPrimitiveStructure")

    def onMBEConvertAsymm(self, event):
        wx.CallAfter(self.canvas.sessionModify, attrib="getAsymmetricStructure")

    def onMBEMoveInCell(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="moveAtomsInCell")

    def onMBEMoveInCellM(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="moveMolCenterInCell")

    def onMBEMoveInMole(self, event):
        wx.CallAfter(self.canvas.sessionModify, no_return=True, attrib="moveAtomsByMolecule")

    def onMBEModify(self, event):
        case = len(self.canvas.vs.selection)
        if case == 2:
            dlg = IMyBondDialog(self)
            dlg.Show()
        elif case == 3:
            dlg = IMyRotateDialog(self)
            dlg.Show()
        elif case == 4:
            dlg = IMyTwistDialog(self)
            dlg.Show()

    def onMBVAlong(self, event):
        winid = event.GetId()
        along = {
            MyID.MVAA: 'a',
            MyID.MVAB: 'b',
            MyID.MVAC: 'c',
            MyID.MVAAs: 'a*',
            MyID.MVABs: 'b*',
            MyID.MVACs: 'c*',
        }
        self.canvas.setCamera(along=along.get(winid))
        self.canvas.Refresh()

    def onMBVSetCenter(self, event):
        self.canvas.setCamera(target='select')
        self.canvas.Refresh()

    def onMBVSync(self, event):
        menu = event.GetEventObject()
        obj = menu.FindItemById(MyID.MVSync)
        if obj.IsChecked():
            self.canvas.camera_sync = True
        else:
            self.canvas.camera_sync = False

    def onMBAMolCompare(self, event):
        raise NotImplementedError
        # namedict = {k: v.data.title for k, v in enumerate(self.canvas.vslist) if v._class == "ViewStruct"}
        # dlg = IMyMolCompareDialog(self, namedict=namedict)
        # dlg.Show()

    def onMBEDelete(self, event):
        self.canvas.onKeyDelete(event)

    def onMBVPageup(self, event):
        self.canvas.onKeyPageup(event)

    def onMBVPagedown(self, event):
        self.canvas.onKeyPagedown(event)

    def onTBCloudopen(self, event):
        from .backend_cli import get_s3_file_string
        with wx.TextEntryDialog(self, 'Open Cloud File', 'address', style=wx.OK | wx.CANCEL) as dlg:
            retcode = dlg.ShowModal()
            if retcode == wx.ID_OK:
                rfp = dlg.GetValue()
                if 'fusial-backend-dev' in rfp:
                    profile = 'fusial-dev'
                elif 'asajj-us-east-1' in rfp:
                    profile = 'fusial-test'
                else:
                    profile = 'fusial-dev'
                string = get_s3_file_string(rfp, profile)
                wx.CallAfter(self.canvas.structLoad, rfp, string, rfp)
        event.Skip()

    def onTBSelect(self, event):
        if self.canvas.vs is None:
            return
        selection = self.statusbar.GetStatusText(0)
        with wx.TextEntryDialog(self, 'Selection', 'Select', selection, wx.OK | wx.CANCEL) as dlg:
            retcode = dlg.ShowModal()
            if retcode == wx.ID_OK:
                selection = [int(x) for x in re.findall(r'\d+', dlg.GetValue())]
                self.canvas.sessionSelect(selection)
        event.Skip()

    def onTBPacking(self, event):
        if self.canvas.vs is None:
            return
        selection = self.canvas.vs.selection
        with wx.TextEntryDialog(self, 'Packing', 'Packing', '1 1 1', wx.OK | wx.CANCEL) as dlg:
            retcode = dlg.ShowModal()
            if retcode == wx.ID_OK:
                packing = [int(x) for x in re.findall(r'\d+', dlg.GetValue())]
                group = selection if len(selection) > 0 else None
                kwargs = {'index': packing, 'group': group}
                wx.CallAfter(self.canvas.sessionModify, attrib="getSupercell", kwargs=kwargs)
        event.Skip()

    def onTBConfig(self, event):
        dlg = IMyFileConfigDialog(self)
        dlg.Show()

    def onTBLattice(self, event):
        dlg = IMyLatticeDialog(self)
        dlg.Show()

    def onTBSymmetry(self, event):
        dlg = IMySymmetryDialog(self)
        dlg.Show()
