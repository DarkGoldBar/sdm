# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
import wx
from wx import glcanvas
from .glcore import GLCore
from .session import ViewSession, ViewStruct
from .components import MyID, FileDropTarget, IMyPopMenu
from . import logger as root_logger
logger = root_logger.getChild('canvas')


class MyGLCanvas(glcanvas.GLCanvas, GLCore):
    """WxPython GLCanvas 窗口模块
    """
    def __init__(self, parent, ID):
        glcanvas.GLCanvas.__init__(self, parent, ID)
        GLCore.__init__(self)

        self.main = None
        self._vsid = 0  # 展示中的结构编号
        self.vslist = []  # 已打开的结构列表
        self.status = 0  # 画布状态 0=未初始化 1=运行中
        self.camera_sync = False  # 是否同步不同结构的镜头
        self.mouse_draged = False
        self.mouse_xy = (0, 0)

        self.popup = IMyPopMenu(self)  # (WxPython) 右键菜单
        self.GLContext = glcanvas.GLContext(self)  # (WxPython) 设置画布焦点为当前对象
        self.SetDropTarget(FileDropTarget(canvas=self))  # (WxPython) 设置鼠标拖拽事件

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.onEraseBackground)
        self.Bind(wx.EVT_SIZE, self.onSize)
        self.Bind(wx.EVT_PAINT, self.onPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.onLeftDown)
        self.Bind(wx.EVT_LEFT_UP, self.onLeftUp)
        self.Bind(wx.EVT_LEFT_DCLICK, self.onLeftDClick)
        self.Bind(wx.EVT_MIDDLE_DOWN, self.onMouseDown)
        self.Bind(wx.EVT_MIDDLE_UP, self.onMouseUp)
        self.Bind(wx.EVT_RIGHT_DOWN, self.onMouseDown)
        self.Bind(wx.EVT_RIGHT_UP, self.onMouseUp)
        self.Bind(wx.EVT_MOTION, self.onMouseMotion)
        self.Bind(wx.EVT_MOUSEWHEEL, self.onMouseWheel)
        self.Bind(wx.EVT_MENU, self.onStyle, id=172, id2=178)

    @property
    def vs(self) -> ViewSession:
        if 0 <= self.vsid < len(self.vslist):
            return self.vslist[self.vsid]

    @property
    def vsid(self) -> int:
        return self._vsid

    @vsid.setter
    def vsid(self, val):
        # 保存当前场景
        if self.vs and not self.camera_sync:
            self.vs.scene['view'] = self.view
            self.vs.scene['scale'] = self.scale
            self.vs.scene['camera'] = self.camera
        # 切换session
        if val is None:
            self._vsid = 0
        elif val < 1:
            self._vsid = 1
        elif val >= len(self.vslist):
            self._vsid = len(self.vslist) - 1
        else:
            self._vsid = val
        # 读取场景
        if not self.camera_sync and self._vsid > 0:
            self.view = self.vs.scene['view']
            self.scale = self.vs.scene['scale']
            self.camera = self.vs.scene['camera']
        # 更新状态栏
        info = self.getSessionInfo() if self._vsid > 0 else ''
        wx.CallAfter(self.main.onChangeSession, info=info)

    def append(self, sess: ViewSession):
        self.vslist.append(sess)
        self.vsid = len(self.vslist) - 1
        self.setCamera()
        self.Refresh()

    def sessionChange(self, data=None, undo=False, redo=False, reset=False, preview=False):
        if preview:
            self.vs.dataPush(data, preview=preview)
        elif undo:
            self.vs.dataUndo()
        elif redo:
            self.vs.dataRedo()
        elif reset:
            self.vs.dataReset()
        else:
            self.vs.dataPush(data)
        self.Refresh()
        wx.CallAfter(self.main.onChangeSession, info=self.getSessionInfo())

    def sessionClose(self):
        if self.vsid <= 0:
            return
        self.vslist.pop(self.vsid)
        self.vsid = self.vsid
        self.Refresh()

    def sessionModify(self, attrib, args=[], kwargs={}, no_return=False, preview=False):
        if no_return:
            new_data = self.vs.data.copy()
            getattr(new_data, attrib)(*args, **kwargs)
        else:
            new_data = getattr(self.vs.data, attrib)(*args, **kwargs)
        self.sessionChange(new_data, preview=preview)

    def structLoad(self, filename, string=None, source='', update=True):
        vs = ViewStruct.load_string(
            filename,
            string,
            info={'source': source},
            style=self.vs.style.copy() if self.vs else dict(),
        )
        vs.makeArray()
        vs.setCamera()
        logger.info('Struct Loaded: {}'.format(filename))
        return self.append(vs) if update else vs

    def structLoadMulti(self, filenames):
        new_vslist = [self.structLoad(fp, source=fp, update=False) for fp in filenames]
        self.vslist += new_vslist[:-1]
        self.append(new_vslist[-1])

    def structSave(self, filename, **kwargs):
        self.vs.data.save(filename, **kwargs)

    def structSetAttr(self, attr_dict):
        new_struct = self.vs.data.copy()
        for attr, value in attr_dict.items():
            setattr(new_struct, attr, value)
        self.sessionChange(new_struct)

    def structCompare(self, ref, com, without_H=True):
        raise NotImplementedError
    #     compare = MoleculeComparsion.compare(self.vslist[ref].data, self.vslist[com].data, without_H=without_H)
    #     vs = ViewCompare()
    #     vs.updateCompare(compare, source="")
    #     vs.style.setStyle("Wireframe")
    #     vs.makeArray()
    #     self.append(vs)

    def sessionSelect(self, selection=None):
        if selection is not None:
            self.vs.selection = selection
        info = self.vs.getSelectInfo()
        wx.CallAfter(self.main.onChangeSelect, info=info)
        self.vs.makeSelectArray()
        self.Refresh()

    def getAtomIdxByMouse(self, x, y) -> int:
        idx = self.vs.calcSelectByRay(*self.getMouseRay(x, y))
        return idx

    def setCamera(self, along=None, target=None):
        vs = self.vs
        if target == 'select':
            if vs.selection:
                target = np.average(vs.coords[vs.selection])
            else:
                target = None
        vs.setCamera(along=along, target=target, old_camera=self.camera)
        self.camera = vs.scene['camera']
        self.view = vs.scene['view']

    def getSessionInfo(self) -> str:
        info_string = '[{}/{}]\n{}\n{}'.format(
            self.vsid,
            len(self.vslist)-1,
            self.vs.info['source'],
            self.vs.getDataInfo()
        )
        return info_string

    def onStyle(self, event):
        """
        MVStyle = 171
        MVSBall = 172
        MVSWire = 173
        MVSStick = 174
        MVSFill = 175
        MVTheme = 176
        MVTDark = 177
        MVTLight = 178
        """
        idmap = {
            MyID.MVSBall: 'Ball-Stick',
            MyID.MVSWire: 'Wireframe',
            MyID.MVSStick: 'Stick',
            MyID.MVSFill: 'Spacefill',
            MyID.MVTDark: 'Dark',
            MyID.MVTLight: 'Light',
        }
        self.vs.setStyle(idmap.get(event.GetId()))
        self.vs.makeArray()
        self.glClearColor(*self.vs.style['C_background'])
        self.Refresh()

    def onEraseBackground(self, event):
        pass  # Do nothing, to avoid flashing on MSW.

    def onSize(self, event):
        if self.IsShown():
            size = self.GetClientSize()
            self.win_wh = (size.width, size.height)
            self.Refresh()

    def onPaint(self, event):
        if self.status == 0:
            self.SetCurrent(self.GLContext)
            self.InitGL()
            self.vslist = [ViewSession()]
            self.status = 1
        self.Draw()
        self.SwapBuffers()

    def onKeyDelete(self, event):
        if len(self.vs.selection) > 0:
            attrib = "delAtoms"
            args = [self.vs.selection]
            wx.CallAfter(self.sessionModify, no_return=True, attrib=attrib, args=args)
        self.sessionSelect([])

    def onKeyPageup(self, event):
        self.vsid -= 1
        self.Refresh()

    def onKeyPagedown(self, event):
        self.vsid += 1
        self.Refresh()

    def onLeftDown(self, event):
        self.SetFocus()
        if len(self.vslist) <= 1:
            return
        self.mouse_xy = event.GetPosition()
        sel_idx = self.getAtomIdxByMouse(*self.mouse_xy)
        self.vs.select(sel_idx, shift=event.ShiftDown())
        self.sessionSelect()

    def onLeftUp(self, event):
        pass

    def onLeftDClick(self, event):
        self.SetFocus()
        if len(self.vslist) <= 1:
            return
        self.mouse_xy = event.GetPosition()
        sel_idx = self.getAtomIdxByMouse(*self.mouse_xy)
        self.vs.select(sel_idx, shift=event.ShiftDown(), dclick=True)
        self.sessionSelect()

    def onMouseDown(self, event):
        self.mouse_draged = False
        if self.vsid <= 0:
            return
        if not self.HasCapture():
            self.CaptureMouse()
        self.SetFocus()

    def onMouseUp(self, event):
        if self.HasCapture():
            self.ReleaseMouse()
        if event.RightUp() and not self.mouse_draged:
            self.PopupMenu(self.popup, event.GetPosition())

    def onMouseMotion(self, event):
        self.mouse_draged = True
        x, y = event.GetPosition()
        if event.Dragging():
            dx = self.mouse_xy[0] - x
            dy = self.mouse_xy[1] - y
            if event.RightIsDown():
                self.RotateCamera(dx, -dy)
            if event.MiddleIsDown():
                self.TranslateCamera(dx, -dy)
            self.Refresh()
        self.mouse_xy = x, y

    def onMouseWheel(self, evt):
        if self.vsid <= 0:
            return
        delta = evt.GetWheelRotation() / (120 * 20)
        self.scale *= (1 + delta)
        self.Refresh()
