# -*- coding:utf-8 -*-
# Author : bochen.li
# Date : 2022/11/09



from sdm.models.molecule import Molecule


class Dimer(Molecule):
    def __init__(self, graph=None, coords=None, title='Untitled', main=None):
        super().__init__(graph, coords, title)
        self.main = main
        if self.main is None:
            self.main = self.findMainGroup()

    def findMainGroup(self):
        