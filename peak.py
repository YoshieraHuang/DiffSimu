# -*- coding: utf-8 -*-
# @Date    : 2017-12-28 21:01:36
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import functools

class Peak(object):
	'''
		diffraction maximum
	'''
	def __init__(self, *, peak_shape = 'Gauss1d', peak_args = (0.1,), tth = None, gamma = None, intn = None):
		super().__init__()

		self.args = peak_args
		self.shape = peak_shape
		self.tth = tth
		self.gamma = gamma
		self.intn = intn
		self.dim = dimension

		self.Gen_peak()

	@property
	def tth(self):
		return self._tth

	@tth.setter
	def tth(self, v):
		if not self.tth is None:
			logging.warning('2theta of %s will be changed!'%(self.__name__))
		self._tth = v

	@property
	def gamma(self):
		return self._gamma

	@gamma.setter
	def gamma(self, v):
		if not self.gamma is None:
			logging.warning('gamma of %s will be changed!'%(self.__name__))
		self._gamma = v

	@property
	def intn(self):
		return self._intn

	@intn.setter
	def intn(self, v):
		if not self.intn is None:
			logging.warning('intensity of %s will be changed!'%(self.__name__))
		self._intn = v

	@property
	def args(self):
		return self._args

	@args.setter
	def args(self, v):
		self._args = v
		self.Gen_peak()

	@property
	def shape(self):
		return self._shape

	@shape.setter
	def shape(self, v):
		self._shape = v
		self.Gen_peak()

	def Guass1d(self, width = 0.1):
		a = self.tth
		r = 8*np.log(2)
		c2 = width**2 / r
		return lambda x: np.exp(-(x - a)**2/(2*c2))

	def Lorentz1d(self, width = 0.1):
		a = self.tth
		g = width
		return lambda x: 1/(np.pi* g * ((x - a) / g)**2)

	def gen_peak(self):
		if self.shape == 'Gauss1d':
			self.peak_frame = self.Gauss1d(*self.args)
		if self.shape == 'Lorentz1d':
			self.peak_frame = self.Lorentz1d(*self.args)

	def frame(self, x):
		return self.intn * self.peak_frame(x)

if __name__ == '__main__':
	pass
