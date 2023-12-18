# -*- coding:utf-8 -*-
import wx
import numpy as np
import yaml
from pathlib import Path


RESOURCE = Path(__file__).parent / 'resource'


class Model3D:
    def loadOBJ(modelname: str, swapyz=True) -> np.ndarray:
        """load .obj model file

        Return:
            (np.ndarray (n * 6)) : x, y, z, norm_x, norm_y, norm_z
        """
        vertices = []
        normals = []
        faces = []
        filepath = RESOURCE / (modelname + '.obj')
        for line in open(filepath, "r", encoding='utf-8'):
            if line.startswith('#'):
                continue
            values = line.split()
            if not values:
                continue
            if values[0] == 'v':
                v = [float(x) for x in values[1:4]]
                if swapyz:
                    v = v[0], v[2], v[1]
                vertices.append(v)
            elif values[0] == 'vn':
                v = [float(x) for x in values[1:4]]
                if swapyz:
                    v = v[0], v[2], v[1]
                normals.append(v)
            elif values[0] == 'f':
                face = []
                norms = []
                for v in values[1:]:
                    w = v.split('/')
                    face.append(int(w[0]))
                    if len(w) >= 3 and len(w[2]) > 0:
                        norms.append(int(w[2]))
                    else:
                        norms.append(0)
                faces.append((face, norms))
        index_array = np.array(faces, int)
        nvert = index_array.shape[0] * index_array.shape[2]
        vertex_normal = np.zeros((nvert, 9), 'f')
        vertex_normal[:, :3] = np.array(vertices, 'f')[index_array[:, 0, :].flatten() - 1]
        vertex_normal[:, 3:6] = np.array(normals, 'f')[index_array[:, 1, :].flatten() - 1]
        return vertex_normal

    ball = loadOBJ('ball')
    stick = loadOBJ('stick')


class FileDropTarget(wx.FileDropTarget):
    def __init__(self, canvas):
        super().__init__()
        self.canvas = canvas

    def OnDropFiles(self, x, y, filenames):
        wx.CallAfter(
            self.canvas.structLoadMulti,
            filenames=filenames
        )
        return True


class IMyPopMenu(wx.Menu):
    UIstring = """
        - text: Style
          item:
            - {id: 172, text: Ball-stick}
            - {id: 173, text: Wireframe}
            - {id: 174, text: Stick}
            - {id: 175, text: Spacefill}
            - {}
            - {id: 177, text: Dark}
            - {id: 176, text: Light}
        - text: View
          item:
            - {id: 161, text: Reset view}
            - {id: 163, text: Set view centre}
            - {}
            - {id: 165, text: Along a axis}
            - {id: 166, text: Along b axis}
            - {id: 167, text: Along c axis}
            - {id: 168, text: Along a* axis}
            - {id: 169, text: Along b* axis}
            - {id: 170, text: Along c* axis}
        - {id: 134, text: Modify}
        """
    UIdata = yaml.safe_load(UIstring)

    def __init__(self, parent):
        super(IMyPopMenu, self).__init__()
        [self.createMenuItem(self, d) for d in self.UIdata]

    def createMenuItem(self, parent: wx.Menu, data: dict):
        if not data:
            parent.AppendSeparator()
        elif "item" in data:
            menu = wx.Menu()
            [self.createMenuItem(menu, d) for d in data["item"]]
            parent.AppendSubMenu(menu, data["text"])
        else:
            parent.Append(data.get("id", -1), item=data["text"])


class MyID:
    MFile = 101
    MFOpen = 102
    MFSave = 103
    MFClose = 104
    MFExport = 105
    MFVasp = 106
    MFCrystal = 107
    MFDftb = 108
    MEdit = 130
    MEUndo = 131
    MERedo = 132
    MEReload = 133
    MEModify = 134
    MEFindbonds = 135
    MEFNone = 136
    MEFEmpiric = 137
    MEFEmpiricpbc = 138
    MEFBabel = 139
    MEFIdatm = 149
    MEConvert = 140
    MECP1 = 141
    MECPrim = 142
    MECAsym = 143
    MEMoveatom = 144
    MEMCell = 145
    MEMCellmol = 146
    MEMMol = 147
    MEDelete = 148
    MView = 160
    MVReset = 161
    MVSync = 162
    MVCenter = 163
    MVAlong = 164
    MVAA = 165
    MVAB = 166
    MVAC = 167
    MVAAs = 168
    MVABs = 169
    MVACs = 170
    MVStyle = 171
    MVSBall = 172
    MVSWire = 173
    MVSStick = 174
    MVSFill = 175
    MVTheme = 176
    MVTDark = 177
    MVTLight = 178
    MVPageup = 179
    MVPagedown = 180
    MAnalysis = 190
    MAMolecompare = 191

    TBSelect = 301
    TBConfig = 302
    TBLattice = 303
    TBSymmetry = 304
    TBPacking = 305
    TBCloudopen = 306

    DlgSymmTxt1 = 201
    DlgSymmTxt2 = 202
    DlgSymmTxt3 = 203
    DlgSymmTxt4 = 204
    DlgLattTxt1 = 301
    DlgLattTxt2 = 302
    DlgLattTxt3 = 303
    DlgLattTxt4 = 304
    DlgLattTxt5 = 305
    DlgLattTxt6 = 306

    ConsoleInput = 99
    ConsoleText = 98
    ConsolePromote = 97
