# -*- coding: utf-8 -*-
# @Date    : 2017-12-28 21:01:36
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging

import lattice as LTTC

class Peak(object):
	'''
		diffraction maximum
	'''
	@property
	def tth(self):
		return self._tth

	@tth.setter
	def tth(self, v):
		self._tth = v

	@property
	def gamma(self):
		return self._gamma

	@gamma.setter
	def gamma(self, v):
		self._gamma = v

	@property
	def wl(self):
		return self._wl

	@wl.setter
	def wl(self, v):
		self._wl = v

	@property
	def intn(self):
		return self._intn

	@intn.setter
	def intn(self, v):
		self._intn = v

	@property
	def args(self):
		return self._args

	@args.setter
	def args(self, v):
		self._args = v
		self.Calc_shape()

	@property
	def shape_name(self):
		return self._shape_name

	@shape_name.setter
	def shape_name(self, v):
		self._shape_name = v
		self.Calc_shape()

	@property
	def index(self):
		return self._index

	@index.setter
	def index(self, v):
		if not isinstance(v, LTTC.index):
			raise TypeError('Must be index!')
		self._index = v

	@property
	def d_spacing(self):
		return self._d_spacing

	@d_spacing.setter
	def d_spacing(self,v):
		if not type(v) is float:
			raise TypeError('Must be float!')
		self._d_spacing = v

	def __init__(self, *, peak_shape_name = 'Gauss1d', peak_args = (0.1,), tth = None, gamma = None, intn = 1):
		super().__init__()

		self._args = peak_args
		self._shape_name = peak_shape_name
		self._tth = tth
		self._gamma = gamma
		self._intn = intn

		self.Calc_shape()

	def Gaussian1d(self, width = 0.1):
		a = self.tth
		c2 = width**2 / (8 * np.log(2))
		return lambda x: np.exp(-(x - a)**2/(2*c2))

	def Lorentzian1d(self, width = 0.1):
		a = self.tth
		g = width/2
		return lambda x: g/(np.pi* ((x - a)**2 + g**2))/(1/(np.pi*g))

	def Calc_shape(self):
		if self.shape_name == 'Gaussian1d':
			self.peak_shape = self.Gaussian1d(*self.args)
		if self.shape_name == 'Lorentzian1d':
			self.peak_shape = self.Lorentzian1d(*self.args)

	def profile(self, x):
		return self.intn * self.peak_shape(x)

	def show(self, start = 0, end = 90, precision = 0.01):
		import matplotlib.pyplot as plt

		fig = plt.figure()
		ax = fig.gca()

		x = np.linspace(start, end, (end - start)/ precision)
		line, = ax.plot(x, self.profile(x))
		ax.set(xlim = [start, end], title = 'peak', xlabel = '2theta (degree)', ylabel = 'Intensity (a. u.)')
		plt.show()


if __name__ == '__main__':
	p = Peak(tth = 30, peak_shape_name = 'Gaussian1d', peak_args = (2,))
	p.show(20,50, 0.1)
