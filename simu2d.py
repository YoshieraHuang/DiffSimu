# -*- coding: utf-8 -*-
# @Date    : 2018-01-29 15:24:25
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging; logging.basicConfig(level = logging.INFO)
import matplotlib.pyplot as plt
from collections import Iterable
import copy
from pathos.multiprocessing import ProcessingPool as Pool

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
			if not isinstance(sx, Iterable):
				sx = np.array(sx)[None]
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
		if not isinstance(sxtal, Iterable):
			sxtal = np.array(sxtal)[None]

		for sxt in sxtal:
			if not isinstance(sxt, SX.SingleXtal):
				raise ValueError('SINGLE XTAL!')

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
			self._inc = self.sx.vec_in_rcp(v).norm

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

		sc_factor = lp(peak.tth) * np.abs(peak.grain.lattice.sc_factor(peak.index, peak.tth, peak.wl))**2
		return peak.wlintn * sc_factor

	def Calc(self, hklrange = (10,10,10), EPS = 0):
		if not hasattr(self, 'sx'):
			raise ValueError('No singleXtal')
		if not hasattr(self, 'xray'):
			raise ValueError('No Xray')
		if not hasattr(self, 'inc'):
			raise ValueError('No Incidence')

		import time
		start = time.clock()
		self.num_tag(self.sx)
		# hkls = list(LTTC.Create_hkls(hklrange = hklrange, lattice = self.sx[0].lattice))
		# print("Calc time 1 : %f s"%(time.clock() - start))
		# self.peaks = self.Calc_peak(hkls, EPS = EPS)
		# print("Calc time 2 : %f s"%(time.clock() - start))
		for sx in self.sx:
			sx.Calc_rcp_space(hklrange)
		print("Calc time 1 : %f s"%(time.clock() - start))
		self.peaks = self.Calc_peak_para()
		print("Calc time 2 : %f s"%(time.clock() - start))

	def Calc_peak_para(self, EPS = 0):
		def give_peak(args):
			hkl, vec = args
			d_spacing = 1 / vec.length
			tth = 2 * vec.betweenangle_rad(inc) - np.pi
			if tth <= np.deg2rad(1) or tth >= np.deg2rad(179):
				return None
			wl = 2 * d_spacing * np.sin(tth/2)
			wlintn = xr.intn_wl(wl)
			if wlintn == 0:
				return None
			k_inc = inc / wl
			k_out = k_inc + vec
			peak = PK.Peak()
			peak.index = hkl
			peak.vec = vec 
			peak.tth, peak.gamma = sx.tthgam_rad(k_out, z = inc, x = self.vx)
			peak.wl = wl
			peak.energy = xr.energy_from_lambda(wl)
			peak.wlintn = wlintn
			peak.d_spacing = d_spacing
			peak.grain = sx
			peak.intn = self.Calc_peak_intn(peak)
			if peak.intn < EPS:
				return None
			return peak

		def isdeplicate(peak, peaks, EPS = np.deg2rad(0.05)):
			if len(peaks) == 0:
				return False

			tth, gamma = peak.tth, peak.gamma
			for i,p in enumerate(peaks):
				if FT.equal(p.tth, tth, EPS) and FT.equal(p.gamma, gamma, EPS):
					if peak.intn > p.intn:
						logging.debug('%s is excluded because of %s'%(peak.index.str, p.index.str))
						# print('%s is excluded because of %s'%(p.index.str, peak.index.str))	
						del peaks[i]
					else:
						logging.debug('%s is excluded because of %s'%(peak.index.str, p.index.str))
						# print('%s is excluded because of %s'%(peak.index.str, p.index.str))				
						return True

			return False

		pool = Pool()
		all_peak = []
		xr = self.xray
		inc = self.inc
		start = time.clock()
		for sx in self.sx:
			peaks = pool.map(give_peak, zip(sx.rcp_space_hkls, sx.rcp_space_vec))
			print("Calc time 3 : %f s"%(time.clock() - start))
			for peak in peaks:
				if (not peak is None) and not isdeplicate(peak,all_peak):
					all_peak.append(peak)
			print("Calc time 4 : %f s"%(time.clock() - start))
		return np.array(all_peak)

	def Calc_peak(self, hkls, EPS = 0):

		def isdeplicate(peak, peaks, EPS = np.deg2rad(0.05)):
			if len(peaks) == 0:
				return False

			tth, gamma = peak.tth, peak.gamma
			for i,p in enumerate(peaks):
				if FT.equal(p.tth, tth, EPS) and FT.equal(p.gamma, gamma, EPS):
					if peak.intn > p.intn:
						logging.debug('%s is excluded because of %s'%(peak.index.str, p.index.str))
						# print('%s is excluded because of %s'%(p.index.str, peak.index.str))	
						del peaks[i]
					else:
						logging.debug('%s is excluded because of %s'%(peak.index.str, p.index.str))
						# print('%s is excluded because of %s'%(peak.index.str, p.index.str))				
						return True

			return False

		peaks = []
		for sx in self.sx:
			inc = self.inc
			xr = self.xray
			for hkl in hkls:
				vec = sx.vec_in_rcp(hkl)
				d_spacing = 1 / vec.length
				tth = 2*vec.betweenangle_rad(inc) - np.pi
				if tth <= np.deg2rad(1) or tth >= np.deg2rad(179):
					continue
				wl = 2 * d_spacing * np.sin(tth/2)
				wlintn = xr.intn_wl(wl)
				if wlintn == 0:
					continue
				k_inc = inc / wl
				k_out = k_inc + vec
				peak = PK.Peak()
				peak.index = hkl
				peak.vec = vec 
				peak.tth, peak.gamma = sx.tthgam_rad(k_out, z = inc, x = self.vx)
				peak.wl = wl
				peak.energy = xr.energy_from_lambda(wl)
				peak.wlintn = wlintn
				peak.d_spacing = d_spacing
				peak.grain = sx
				peak.intn = self.Calc_peak_intn(peak)
				if peak.intn > EPS and not isdeplicate(peak,peaks):
					logging.debug('PEAK info:index: %r, d: %f, lambda:%f, tth: %f, intn: %f'%(hkl, d_spacing, wl, np.rad2deg(tth), peak.intn))
					peaks.append(peak)
		return np.array(peaks)

	def rotate_sx(self, q):
		self.inc = self.inc.rotate_by(q.inverse)

	def num_tag(self, sxs):
		for i,sx in enumerate(sxs):
			sx.number = i+1

	def show(self):
		fig = plt.figure()
		ax = fig.gca()
		ax.set(xlabel = r'$2\theta (^\circ)$', ylabel = r'$\gamma (^\circ)$', title = '2D Pattern')
		tth = [np.rad2deg(peak.tth) for peak in self.peaks]
		gamma = [np.rad2deg(peak.gamma) for peak in self.peaks]
		intn = [peak.intn  for peak in self.peaks]
		size = intn / np.linalg.norm(intn) * 200
		for peak in self.peaks:
			print(peak.index, peak.grain.number, np.rad2deg(peak.tth), np.rad2deg(peak.gamma), peak.vec, peak.wl, peak.intn)
		ax.scatter(tth, gamma, size)
		plt.show()

	def save(self, filename):
		if not hasattr(self, 'peaks'):
			print('No peaks!')
			return

		with open(filename, 'w') as f:
			f.writelines('%d\n'%(self.peaks.shape[0]))
			f.writelines('index grain tth gamma kx ky kz d_spacing intn lambda\n')
			for p in self.peaks:
				f.writelines('%s %d %f %f %f %f %f %f %f %f\n'%(p.index.str, p.grain.number, np.rad2deg(p.tth), np.rad2deg(p.gamma), p.vec[0], p.vec[1], p.vec[2], p.d_spacing, p.intn, p.wl))

	def copy(self):
		return copy.copy(self)

	def __add__(self,other):
		if not isinstance(other, Pattern2d):
			raise ValueError('Only \'Pattern2d\' can be added!')

		if not hasattr(other, 'peaks'):
			raise ValueError('No peaks information!')

		if not self.xray == other.xray:
			raise ValueError('XRAY should be identical!')

		new = self.copy()
		new.sx = np.hstack((self.sx, other.sx))
		new.num_tag(self.sx)
		new.peaks = np.hstack((self.peaks, other.peaks))

		return new



if __name__ == '__main__':
	import time
	start = time.clock()
	lttc = LTTC.Lattice(material = 'Si')
	sx = SX.SingleXtal(lttc, x = (1,0,0), z = (0,0,1))
	print('time1: %f s'%(time.clock()-start))
	# sx.rotate_by(axis = (1,0,0), degree = -6)
	# xr = XR.Xray(wavelength = np.linspace(0.5, 0.6, 1000))
	xr = XR.Xray('White')
	inc = Vector(0,0,1)
	p = Pattern2d(sx = sx, xray = xr, inc = inc)
	print('time2: %f s'%(time.clock()-start))
	p.Calc(hklrange = (15,15,15))
	print('time3: %f s'%(time.clock()-start))
	# p.show()
	p.save('peaks.txt')
	print('time: %f s'%(time.clock()-start))

