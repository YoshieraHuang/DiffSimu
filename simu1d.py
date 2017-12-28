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

	def calc(self, *, precision = 0.01, norm = True, range_2th = (0, np.pi/2)):
		profile = Profile1d()
		for wl, intn_wl in self.xray.spectrum:
			hkls = LTTC.Gen_hklfamilies(lattice = self.xtal.lattice)
			peaks = self.bragg(hkls, wl, intn)
			for peak in peaks:
				if peak.tth < range_2th[0]:
					continue
				if peak.tth > range_2th[1]:
					break
				profile.Add_peak(peak)
			if norm:
				profile.norm()

		self.profile = profile

	def bragg(self, hkls, wl, intn, peak_width = 0.1, peak_shape = 'Gauss'):
		for hkl, d_spacing, factor in zip(hkls, self.xtal.d_spacings(hkls), self.xtal.factor(hkls)):
			if factor is None:
				continue
			peak = Peak(peak_width = peak_width, peak_shape = peak_shape)
			tth = np.arcsin(wl/(2*d_spacing))
			peak.tth =tth
			peak.intn = intn*factor
			yield peak

	

def Profile1d(object):
	'''
		1D profile of diffraction
	'''
	def __init__(self):
		super().__init__()


if __name__ == '__main__':
	pass
