# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 19:21:58
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging
import matplotlib.pyplot as plt


import lattice as LTTC
import polyXtal as PX
import xray as XR

class DiffSimu_1D(object):
	'''
		1D Diffraction Simulation for poly crystal 
	'''
	def __init__(self):
		super(DiffSimu_1D, self).__init__()

	@property
	def xtal(self):
		if not hasattr(self, '_xtal'):
			raise ValueError('there are no lattice information')
		return self._xtal

	@xtal.setter
	def xtal(self, ls):
		if isinstance(ls, LTTC.Lattice):
			px = PX.PolyXtal(ls)
		elif isinstance(ls, PX.PolyXtal):
			px = ls
		else:
			raise TypeError('Arguments must be Lattice or PolyXtal')

		self._xtal = lttc

	@property
	def xray(self):
		if not hasattr(self, '_xray'):
			raise ValueError('there are no x-ray information')
		return self._xray

	@xray.setter
	def xray(self, xr):
		if not isinstance(xr, XR.Xray):
			raise TypeError('Arguments must be Xray')

		self._xray = xr

	def calc(self, *, precision = 0.01, norm = True):
		for wl, intn in self.xray.spectrum:
			for hkl in LTTC.Gen_hklfamilies():
				


	def show(self):
		plt.plot(self.tth, self.intensity)
		plt.title('1D simulation')
		plt.xlabel('2theta (degrees)')
		plt.ylabel('Intensity (a. u.)')
		plt.show()

if __name__ == '__main__':
	pass
