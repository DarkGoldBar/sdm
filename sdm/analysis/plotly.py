# coding: utf-8
# plotly == 5.11.0
# 2022/10/31

from itertools import product
from copy import deepcopy
import numpy as np
import plotly.graph_objs as go
import plotly.offline as offline


class Vis:
    _init_notebook_mode = False

    def __init__(self):
        self.data = []
        self.layout = deepcopy(_template_layout)
        self._atom_scale = 7.0
        if not self._init_notebook_mode:
            offline.init_notebook_mode()
            self._init_notebook_mode = True

    @classmethod
    def plot(cls, object):
        o = cls()
        o.add_mol(object)
        if hasattr(object, 'lattice'):
            o.add_cell(object.lattice.matrix)
        o._plot()
        return o

    def _plot(self):
        fig = go.Figure(data=self.data, layout=self.layout)
        offline.iplot(fig)

    def add_cell(self, matrix: np.ndarray):
        def get_edge(x):
            return [x[0], x[1], x[3], x[2], x[0], x[4], x[5], x[7], x[6], x[4], None, x[1], x[5], None, x[2], x[6], None, x[3], x[7]]
        mp = np.array(list(product([0, 1], [0, 1], [0, 1]))) @ matrix
        mx, my, mz = map(get_edge, zip(*mp))
        data_cell = go.Scatter3d(x=mx, y=my, z=mz, mode='lines',
                                 line={'color': '#202020', 'width': 2},
                                 hoverinfo='skip')
        self.data.append(data_cell)

    def add_mol(self, mol, without_H=False, bond_width=3, text='{i}:{t}',
                atom_color=None, atom_size=None, bond_color=None):
        species = {}
        for i in mol.atoms:
            e = mol.atoms[i]['element']
            t = mol.atoms[i]['title']
            p = mol.coords[i]
            if without_H and e.symbol == 'H':
                continue
            if e.symbol not in species:
                species[e.symbol] = {
                    'coords': [],
                    'text': [],
                    'color': atom_color if atom_color else '#' + e.vcolor,
                    'size': atom_size if atom_size else e.vradius,
                }
            species[e.symbol]['coords'].append(p)
            species[e.symbol]['text'].append(text.format(i=i, t=t))

        for sp in species.values():
            self.add_points(**sp)

        bonds = list(mol.bonds)
        self.add_lines(mol.coords, bonds, bond_color, width=bond_width)

    def add_points(self, coords, color='#FFAAFF', size=None, marker={}, **kwargs):
        x, y, z = np.array(coords).T
        size = 1 if size is None else size
        size *= self._atom_scale
        points = deepcopy(_template_atom)
        points['marker'].update(color=color, size=size, **marker)
        points.update(x=x, y=y, z=z, **kwargs)
        self.data.append(go.Scatter3d(points))

    def add_lines(self, coords, edges, color='#808080', width=2):
        x = [x for i, j in edges for x in (coords[i, 0], coords[j, 0], None)][:-1]
        y = [x for i, j in edges for x in (coords[i, 1], coords[j, 1], None)][:-1]
        z = [x for i, j in edges for x in (coords[i, 2], coords[j, 2], None)][:-1]
        lines = deepcopy(_template_line)
        lines['line'].update(color=color, width=width)
        lines.update(x=x, y=y, z=z)
        self.data.append(go.Scatter3d(lines))

_template_axis = {'showbackground': False,
                  'showgrid': True,
                  'showline': True,
                  'zeroline': True,
                  'showticklabels': False,
                  'rangemode': 'normal',
                  'tickmode': 'linear',
                  'dtick': 1,
                  'title': ''}

_template_atom = {'hoverinfo': 'x+y+z+text',
                  'mode': 'markers',
                  'text': [],
                  'x': [],
                  'y': [],
                  'z': [],
                  'marker': {'color': '#F8F',
                             'colorscale': 'Rainbow',
                             'size': 2,
                             'symbol': 'circle',
                             'line': {'color': '#444',
                                      'width': 0.5}}}

_template_line = {'hoverinfo': 'skip',
                  'mode': 'lines',
                  'x': [],
                  'y': [],
                  'z': [],
                  'line': {'color': '#888',
                           'colorscale': 'Rainbow',
                           'width': 2}}

_template_layout = {'title': '',
                    'width': 800,
                    'height': 600,
                    'showlegend': False,
                    'scene': {'xaxis': _template_axis.copy(),
                              'yaxis': _template_axis.copy(),
                              'zaxis': _template_axis.copy(),
                              'camera': {'projection': {'type': 'orthographic'}},
                              'aspectmode': 'data'}}
