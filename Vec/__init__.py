# -*- coding: utf-8 -*-
# @Date    : 2017-12-15 09:34:47
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import warnings

EPS = 1e-5
class Vectors(np.ndarray):
	"""3-dimensional vectors"""
	def __new__(cls, input_array):
		return np.asarray(input_array).view(cls)

	@property
	def length(self):
		return float(np.sqrt(np.sum(self **2)))

	@property
	def norm(self):
		if self.shape == (3,):
			if self.length == 0:
				raise ValueError("Zero Vectors!")
			return self / self.length
		elif self.shape == (2,):
			ele3 = 1 - self.length
			if ele3 < 0:
				raise ValueError("Can't convert to a normalized vector: Too long!")
			v_new = np.append(self,ele3)
			return Vectors(v_new)
		else:
			raise ValueError("Can't convert to a normalized vector: Wrong items!")

	def betweenangle_rad(self,other):
		vother = Vectors(other)

		return float(np.arccos(np.dot(self,other)/(self.length * vother.length)))

	def betweenangle_deg(self,other):
		return np.rad2deg(self.betweenangle_rad(other))

	def dot(self,other):
		return float(np.sum((s*o for s,o in zip(self,other))))

	def cross(self,other):
		if self.parallel(other):
			raise ValueError("Two vectors are parallel")

		return Vectors(np.cross(self,other))

	def parallel(self,other):
		'''
			When two Vectors is not parallel, 0 will be returned
			When two Vectors is parallel and have same direction, 1 will be returned.
			When two Vectors is reverse-parallel, 2 will be returned
			When directions do not matter, this function can also work 
		'''
		angle = self.betweenangle_rad(other)
		if angle < EPS:
			return 1

		if np.pi - angle < EPS:
			return 2

		return 0

	def perpendicular(self, other):
		if np.fabs(self.betweenangle_rad(other) - np.pi/2) < EPS:
			return True
		return False

	def tthgam_rad(self,*,z = (0,0,1), x = (1,0,0)):

		z = Vectors(z); x = Vectors(x)

		if x.betweenangle_rad(z) != 0:
			warnings.warn("x is not perpendicular to z, result may be weird!",DeprecationWarning)

		v = self.norm
		tth = v.betweenangle_rad(z)
		v_p2z =  v - np.cos(tth)*z
		gamma = v_p2z.betweenangle_rad(x)
		return (tth, gamma)

	def tthgam_deg(self,**kw):
		return np.rad2deg(self.tthgam_rad(**kw))