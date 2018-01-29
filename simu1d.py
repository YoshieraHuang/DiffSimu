# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 19:21:58
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging; logging.basicConfig(level = logging.INFO)
import matplotlib.pyplot as plt
from multiprocessing.dummy import Pool 
import functools

import lattice as LTTC
import polyXtal as PX
import xray as XR
import peak as PK
import myfunctools as FT
import LUT

class Profile1D(object):
	'''
		1D Diffraction Simulation for poly crystal 
	'''
	def __init__(self):
		super(Profile1D, self).__init__()

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

		self._xtal = px

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

	def Add_peaks(self, peaks):
		if not hasattr(self, 'peaks'):
			self.peaks = peaks
		else:
			self.peaks = np.concatenate((self.peaks, peaks))


	def Calc_peak(self, range_2th = (0, 90)):
		hkls = FT.tolist(LTTC.Gen_hklfamilies(lattice = self.xtal.lattice))
		pool = Pool(4)
		peaks = pool.map(lambda x: FT.tolist(self.bragg(x, hkls, range_2th = range_2th)), self.xray.spectrum)
		pool.close()
		self.Add_peaks(np.concatenate(peaks))

	def Gen_profile(self, peak_shape = 'Gaussian1d', peak_args = (np.deg2rad(0.1),)):
		p = []
		for peak in self.peaks:
			peak.shape = peak_shape
			peak.args = peak_args
			p.append(peak.frame)

		self.profile = p

	def Calc_profile(self, start = 0, end = 90, precision = 0.1):
		self.tth = np.deg2rad(np.arange(start, end, precision))

		intn = 0
		for p in self.profile:
			intn += p(self.tth)

		self.intn = intn

	def calc_peak_intn(self, peak, wl, wlintn):
		def lp(tth):
			return (1 + np.cos(tth)**2)/(np.cos(tth/2)*(np.sin(tth/2)**2))

		sc_factor = peak.index.multiplicity * lp(peak.tth) * np.abs(self.xtal.lattice.sc_factor(peak.index, peak.tth, wl))**2
		return wlintn * sc_factor

	def bragg(self, w, hkls, range_2th = (0,90)):
		def where(arr, value):
			idx = np.argwhere(arr == value)
			if idx.shape == (0,1):
				return 0
			return idx[0,0] + 1

		wl, wlintn = w
		r_tth = np.deg2rad(range_2th)
		tag = where(self.xray.tag_wl, wl)
		if not hasattr(self, 'tags'):
			self.tags = []
		for hkl, d_spacing in zip(hkls, self.xtal.d_spacing(hkls)):
			sinth = wl/(2*d_spacing)
			if sinth < -1 or sinth > 1:
				continue
			tth = 2*np.arcsin(sinth)
			if tth < r_tth[0] or tth > r_tth[1]:
				continue
			peak = PK.Peak()
			peak.tth = tth
			peak.index = hkl
			peak.d_spacing = d_spacing
			peak.intn = self.calc_peak_intn(peak, wl, wlintn)
			peak.wl = wl
			if tag:
				self.tags.append((peak,tag))
			logging.debug('PEAK info:index: %r, d: %f, lambda:%f, tth: %f, intn: %f, tag: %r'%(hkl, d_spacing, wl, tth, peak.intn, tag))
			yield peak

	def show(self, scale_tag_intn = 1):
		fig = plt.figure()
		ax = fig.gca()
		tth = self.tth if min(self.tth) > 4 else np.rad2deg(self.tth) # radian or degree
		ax.plot(tth, self.intn, '-')
		ax.set(xlabel = r'$2\theta (^\circ )$', ylabel = 'Intensity (a. u.)', title = '1D Profile')
		for peak,tag in self.tags:
			colors = LUT.colors
			tth = np.rad2deg(peak.tth)
			intn = peak.intn * scale_tag_intn
			ax.plot((tth, tth), (0, intn), '-' + colors[tag])
			ax.annotate(peak.index.str, xy = (tth, intn), ha = 'center', va = 'bottom', size = 8)
		plt.show()

	def save_profile(self, filename):
		np.savetxt(filename, np.concatenate((self.tth, self.intn)))

if __name__ == '__main__':
	p = Profile1D()
	l = LTTC.Lattice(material = 'Cu')
	crystal = PX.PolyXtal(l)
	# xr = XR.Xray(wavelength = np.arange(0.4,0.43,0.0005))
	# xr = XR.Xray(wavelength = 0.4)
	# xr = XR.Xray(wavelength = 0.5)
	xr = XR.Xray(filename = 'u18_gap12mm.txt', islambda = False, EPS = 1e13)
	xr.Gen_tag()
	# xr.Add_wavelength(0.41)
	# xr.Gen_tag(range_wl = (0.3,0.1))
	xr.show()
	p.xray = xr
	# xr.show()
	p.xtal = crystal

	# p.xray = XR.Xray(filename = 'u18_gap12mm.txt', islambda = False, EPS = 1e12)
	p.Calc_peak(range_2th = (10,30))
	p.Gen_profile()
	p.Calc_profile(start = 10, end = 30, precision = 0.001)
	p.show(scale_tag_intn = 12)
	# p.save_profile('peaks.dat')