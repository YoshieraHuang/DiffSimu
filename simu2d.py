# -*- coding: utf-8 -*-
# @Date    : 2018-01-29 15:24:25
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging; logging.basicConfig(level = logging.INFO)
import matplotlib.pyplot as plt

import singleXtal as SX
import lattice as LTTC
import xray as XR
import peak as PK
import LUT
import myfunctools as FT
from Vec import Vector

class Pattern2d(object):

	def __init__(self, sx = None, xray = None, inc = None, x = None):
		super(Pattern2d, self).__init__()
		if not sx is None:
			self._sx = sx
		if not inc is None:
			self._inc = inc
		if not xray is None:
			self._xray = xray
		if not x is None:
			self._vx = x

	@property
	def xray(self):
		return self._xray

	@xray.setter
	def xray(self, xr):
		if not isinstance(xr, XR.xray):
			raise ValueError('XRAY!')
		self._xray = xr

	@property
	def sx(self):
		return self._sx

	@sx.setter
	def sx(self, sxtal):
		if not isinstance(sxtal, SX.SingleXtal):
			raise ValueError('SINGLEXTAL!')
		self._sx = sxtal

	@property
	def inc(self):
		return self._inc

	@inc.setter
	def inc(self, v):
		if not isinstance(v, (Vector, LTTC.index)):
			raise ValueError('VECTOR or INDEX')

		if isinstance(v, Vector):
			self._inc = v.norm

		if isinstance(v, LTTC.index):
			self._inc = self.sx.Vec_in_sx(v).norm

	@property
	def vx(self):
		if not hasattr(self, '_vx'):
			if self.inc[1] == 0 and self.inc[0] == 0:
				self._vx = Vector(1,0,0)
				return self._vx
			self._vx = Vector((self.inc[1], -self.inc[0], 0))
			if self.inc[1] < 0:
				self._vx = -self._vx

		return self._vx

	@vx.setter
	def vx(self, v):
		self._vx = v

	def Calc_peak_intn(self, peak):
		def lp(tth):
			return (1 + np.cos(tth)**2)/(np.cos(tth/2)*(np.sin(tth/2)**2))

		sc_factor = lp(peak.tth) * np.abs(self.sx.lattice.sc_factor(peak.index, peak.tth, peak.wl))**2
		return peak.wlintn * sc_factor

	def Calc(self, hklrange = (10,10,10)):
		if not hasattr(self, 'sx'):
			raise ValueError('No singleXtal')
		if not hasattr(self, 'xray'):
			raise ValueError('No Xray')
		if not hasattr(self, 'inc'):
			raise ValueError('No Incidence')

		self.peaks = list(self.Calc_peak(hklrange = hklrange))

	def Calc_peak(self, hklrange = (10,10,10)):
		sx = self.sx
		inc = self.inc
		xr = self.xray
		hkls = FT.tolist(LTTC.Gen_hkls(hklrange = hklrange, lattice = sx.lattice))
		for hkl, vec in zip(hkls, sx.vec_in_sx(hkls)):
			d_spacing = 1 / vec.length
			tth = 2*vec.betweenangle_rad(inc) - np.pi
			if tth <= 0 or tth >= np.pi:
				continue
			wl = 2 * d_spacing * np.sin(tth/2)
			wlintn = xr.intn_wl(wl)
			if wlintn == 0:
				continue
			k_inc = inc / wl
			k_out = k_inc + vec
			peak = PK.Peak()
			peak.index = hkl
			peak.tth, peak.gamma = sx.tthgam_rad(k_out, z = -inc, x = self.vx)
			peak.wl = wl
			peak.wlintn = wlintn
			peak.d_spacing = d_spacing
			peak.intn = self.Calc_peak_intn(peak)
			logging.debug('PEAK info:index: %r, d: %f, lambda:%f, tth: %f, intn: %f'%(hkl, d_spacing, wl, np.rad2deg(tth), peak.intn))
			yield peak

	def rotate_sx(self, q):
		self.inc = self.inc.rotate_by(q.inverse)

	def show(self):
		fig = plt.figure()
		ax = fig.gca()
		ax.set(xlabel = r'$2\theta (^\circ)$', ylabel = r'$\gamma (^\circ)$', title = '2D Pattern')
		tth = [np.rad2deg(peak.tth) for peak in self.peaks]
		gamma = [np.rad2deg(peak.gamma) for peak in self.peaks]
		intn = [peak.intn  for peak in self.peaks]
		size = intn / np.linalg.norm(intn) * 200
		for peak in self.peaks:
			print(peak.index, np.rad2deg(peak.tth), np.rad2deg(peak.gamma), peak.wl, peak.intn)
		ax.scatter(tth, gamma, size)
		plt.show()


if __name__ == '__main__':
	lttc = LTTC.Lattice(material = 'Ta')
	sx = SX.SingleXtal(lttc, x = (-1,1,0), z = (1,1,1))
	sx.rotate_by(axis = (1,0,0), degrees = -6)
	xr = XR.Xray(wavelength = np.linspace(0.5, 0.6, 1000))
	inc = Vector(0,0,-1)
	p = Pattern2d(sx = sx, xray = xr, inc = inc)
	p.Calc(hklrange = (5,5,5))
	p.show()


