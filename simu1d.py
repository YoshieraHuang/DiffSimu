# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 19:21:58
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging; logging.basicConfig(level = logging.DEBUG)
import matplotlib.pyplot as plt
from multiprocessing.dummy import Pool 
import functools

import lattice as LTTC
import polyXtal as PX
import xray as XR
import peak as PK
import myfunctools as FT
import LUT
from Vec import Vector

class Profile1D(object):
	'''
		1D Diffraction Simulation for poly crystal 
	'''
	def __init__(self, xray = None, px = None):
		super(Profile1D, self).__init__()

		if not xray is None:
			self._xray = xray

		if not px is None:
			self._px = px

	@property
	def px(self):
		if not hasattr(self, '_px'):
			raise ValueError('there are no lattice information')
		return self._px

	@px.setter
	def px(self, ls):
		if isinstance(ls, LTTC.Lattice):
			px = PX.PolyXtal(ls)
		elif isinstance(ls, PX.PolyXtal):
			px = ls
		else:
			raise TypeError('Arguments must be Lattice or PolyXtal')

		self._px = px

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
			self.peaks = np.concatenate((self.peaks,peaks))
			
		# print(self.peaks)


	def Calc_peak(self, range_2th = (0, 90)):
		self.range_2th = range_2th
		hkls = FT.tolist(LTTC.Gen_hklfamilies(lattice = self.px.lattice))
		pool = Pool(4)
		peaks = pool.map(lambda x: FT.tolist(self.bragg(x, hkls, range_2th = self.range_2th)), self.xray.spectrum)
		pool.close()
		if len(peaks) == 0:
			return
		if len(peaks) > 1:
			peaks = np.concatenate(peaks)
		else:
			peaks = np.array(peaks[0])
		# print(peaks)
		self.Add_peaks(peaks)

	def Gen_profile(self, peak_shape = 'Gaussian1d', peak_args = (np.deg2rad(0.1),)):
		p = []
		for peak in self.peaks:
			peak.shape = peak_shape
			peak.args = peak_args
			p.append(peak.frame)

		self.profile = p

	def Calc_profile(self, range_2th = None, precision = 0.1):
		if range_2th is None:
			range_2th = self.range_2th

		self.tth = np.deg2rad(np.arange(range_2th[0], range_2th[1], precision))

		intn = 0
		for p in self.profile:
			intn += p(self.tth)

		self.intn = intn

	def calc_peak_intn(self, peak):
		def lp(tth):
			return (1 + np.cos(tth)**2)/(np.cos(tth/2)*(np.sin(tth/2)**2))

		sc_factor = peak.index.multiplicity * lp(peak.tth) * np.abs(self.px.lattice.sc_factor(peak.index, peak.tth, peak.wl))**2
		return peak.wlintn * sc_factor

	def bragg(self, w, hkls, range_2th = (0,90)):

		wl, wlintn = w
		r_tth = np.deg2rad(range_2th)
		tag = FT.where(self.xray.tag_wl, wl)
		if not hasattr(self, 'tags'):
			self.tags = []
		for hkl, d_spacing in zip(hkls, self.px.d_spacing(hkls)):
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
			peak.wl = wl
			peak.wlintn = wlintn
			peak.intn = self.calc_peak_intn(peak)
			if tag:
				self.tags.append((peak,tag))
			logging.debug('PEAK info:index: %r, d: %f, lambda:%f, tth: %f, intn: %f, tag: %r'%(hkl, d_spacing, wl, np.rad2deg(tth), peak.intn, tag))
			yield peak

	def Calc(self, range_2th = (0,90), precision = 0.1, peak_shape = 'Gaussian1d' , peak_args = (np.deg2rad(0.1),)):
		self.Calc_peak(range_2th = range_2th)
		self.Gen_profile(peak_shape = peak_shape, peak_args= peak_args)
		self.Calc_profile(range_2th = range_2th, precision = precision)

	def Add_geometry(self, inc = Vector(0,0,-1), x = Vector(1,0,0)):
		self.inc = inc
		self.vx = x

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
		np.savetxt(filename, np.vstack((self.tth, self.intn)).T)

if __name__ == '__main__':
	Ta = LTTC.Lattice(material = 'Ta')
	Cu = LTTC.Lattice(material = 'Cu')
	Mg = LTTC.Lattice(material = 'Mg')
	# sx = SX.SingleXtal(Ta, z = (1,1,1), x = (-1,1,0))
	# sx.rotate_by(axis = (1,0,0), degrees = -6)
	# sx.strain1d(axis = (1,0,0), ratio = -0.03)
	px = PX.PolyXtal(Mg)
	# xr = XR.Xray(filename = 'u18_gap12mm.txt', islambda = False, EPS = 1e13)
	# xr.show()
	xr = XR.Xray(wavelength = 0.52)
	inc = Vector(0,0,-1)
	# inc = inc.rotate_by(axis = (1,0,0), degree = 10)
	# p = S2D.Pattern2d(sx = sx, xray = xr, inc = Vector(0,0,-1))
	# p.Calc(hklrange = (5,5,5))
	p = Profile1D(xray = xr, px = px)
	p.Calc(range_2th = (10, 20), precision = 0.001)
	p.show()