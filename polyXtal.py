# -*- coding: utf-8 -*-
# @Date    : 2017-12-26 10:18:07
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging


import lattice as LTTC
import xray as XR
import myfunctools as FT

class PolyXtal(object):
	'''
		Class for poly crystal
		Attributes:
			lattice: lattice information of this poly crystal
			reciprocal_matrix: a copy of reciprocal_matrix attribute in Lattice

		Method:
			d_spacing: d_spacing for specific hkl
	'''
	def __init__(self, lttc):
		super(PolyXtal, self).__init__()

		self.lattice = lttc
		self.reciprocal_matrix = self.lattice.reciprocal_matrix

	def d_spacing(self, hklfamilies):
		if isinstance(hklfamilies, LTTC.Familyindex):
			hklfamilies = np.array((hklfamilies))[None] ## add a dimension to make it iterable

		return self.lattice.d_spacing(hklfamilies)

	def D_spacing(self, hklfamilies):
		return FT.tolist(self.d_spacing(hklfamilies)) 

if __name__ == '__main__':
	Cu = LTTC.Lattice(material = 'Cu')
	Mg = LTTC.Lattice(material = 'Mg')
	X = PolyXtal(Mg)
	hkl = LTTC.Familyindex((1,0,0)), LTTC.Familyindex((0,0,2)), LTTC.Familyindex((1,0,1))
	print(X.D_spacing(hkl))
