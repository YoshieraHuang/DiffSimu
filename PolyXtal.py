# -*- coding: utf-8 -*-
# @Date    : 2017-12-26 10:18:07
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging
import Lattice
import HKL

class PolyXtal(object):
	'''
		Class for poly crystal
	'''
	def __init__(self, lattice):
		super(PolyXtal, self).__init__()

		self.lattice = lattice

	def d_spacing(self, i):
		if not isinstance(i, (HKL.index, HKL.Familyindex)):
			logging.error('Unexpected Class!')

		if isinstance(i, HKL.index):
			logging.warning('For d_spacing of PolyXtal, input index should be Familyindex!')

		return self.lattice.d_spacing(i.sons()[0])
