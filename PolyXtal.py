# -*- coding: utf-8 -*-
# @Date    : 2017-12-26 10:18:07
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging
import Lattice

class PolyXtal(object):
	'''
		Class for poly crystal
		Attributes:
			lattice: lattice information of this poly crystal

		Method:
			d_spacing: d_spacing for specific hkl
	'''
	def __init__(self, lattice):
		super(PolyXtal, self).__init__()

		self.lattice = lattice

	def d_spacing(self, hklfamilies):
		if isinstance(hklfamilies, Lattice.Familyindex):
			hklfamilies = np.array((hklfamilies))[None] ## add a dimension to make it iterable

		results = self.lattice.d_spacing([hklfamily.sons()[0] for hklfamily in hklfamilies])
		logging.debug('d_spacing is %s'%(str(results)))
		return results



if __name__ == '__main__':
	Cu = Lattice.Lattice(material = 'Cu')
	X = PolyXtal(Cu)
	hkl = Lattice.Familyindex((2,0,0)), Lattice.Familyindex((3,3,1)), Lattice.Familyindex((3,1,1))
	print(X.d_spacing(hkl))
