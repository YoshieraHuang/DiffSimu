# -*- coding: utf-8 -*-
# @Date    : 2017-12-15 09:34:47
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
from pyquaternion import Quaternion
import warnings

EPS = 1e-5
class Vector(np.ndarray):
	"""3-dimensional Vector"""
	def __new__(cls, *args, tth = None, gamma = None, z = (0,0,1), x = (1,0,0)):
		if args:
			if len(args) == 1:
				input_array = args[0]
			elif len(args) in (2,3):
				if all([isinstance(arg, (int, float)) for arg in args]):
					input_array =  args
				else:
					raise ValueError('Wrong Parameters')
			else:
				raise ValueError('Too many Parameters')
			return np.asarray(input_array).view(cls)

		if tth is None or gamma is None:
			raise ValueError('Empty Parameters')

		if np.dot(x, z) != 0:
			raise ValueError('x and z is not perpendicular')
		y = np.cross(z, x)
		frame = np.vstack((x,y,z))
		v0 = np.array((np.sin(tth)*np.cos(gamma), np.sin(tth)*np.sin(gamma), np.cos(tth)))
		v = np.dot(v0, frame)
		return v.view(cls)

	def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
		args = []
		in_no = []
		for i,input_ in enumerate(inputs):
			if isinstance(input_, Vector):
				args.append(input_.view(np.ndarray))
			else:
				args.append(input_)

		outputs = kwargs.pop('out', None)
		out_no = []
		if outputs:
			out_args = []
			for j, output in enumerate(outputs):
				if isinstance(output, Vector):
					out_args.append(output.view(np.ndarray))
				else:
					out_args.append(output)
			kwargs['out'] = tuple(out_args)
		else:
			outputs = (None,)* ufunc.nout

		results = super(Vector, self).__array_ufunc__(ufunc, method, *args, **kwargs)

		if results is NotImplemented:
			return NotImplemented
		if method == 'at':
			return

		if ufunc.nout == 1:
			results = (results,)

		results = tuple((np.asarray(result).view(Vector)
							if output is None else output)
							for result, output in zip(results, outputs))
		 
		return results[0] if len(results) == 1 else results

	@property
	def length(self):
		return float(np.sqrt(np.sum(self **2)))

	@property
	def norm(self):
		if self.shape == (3,):
			if self.length == 0:
				raise ValueError("Zero Vector!")
			return self / self.length
		elif self.shape == (2,):
			ele3 = 1 - self.length
			if ele3 < 0:
				raise ValueError("Can't convert to a normalized vector: Too long!")
			v_new = np.append(self, ele3)
			return Vector(v_new)
		else:
			raise ValueError("Can't convert to a normalized vector: Unknown items!")

	def betweenangle_rad(self,other):
		vother = Vector(other)
		if self.length == 0 or other.length == 0:
			raise ValueError("Zero Vector!")

		return float(np.arccos(self.dot(vother)/(self.length * vother.length)))

	def betweenangle_deg(self,other):
		return np.rad2deg(self.betweenangle_rad(other))

	def dot(self,other):
		return float(np.sum(np.multiply(self, other)))

	def cross(self,other):
		if self.isparallel(other):
			raise ValueError("Two Vector are parallel")

		return Vector(np.cross(self,other))

	def isparallel(self,other):
		if self.betweenangle_rad(other) < EPS:
			return True
		return False

	def isperpendicular(self, other):
		if np.fabs(self.betweenangle_rad(other) - np.pi/2) < EPS:
			return True
		return False

	def tthgam_rad(self, z = (0,0,1), x = (1,0,0)):
		def mygamma(x,y):
			if x == 0:
				if y > 0:
					return np.pi/2
				else:
					return 3*np.pi/2

			gamma = np.arctan(y/x)

			if (y >= 0) and (x >= 0):
				return gamma
			if (y >=0) and (x < 0):
				return np.pi + gamma
			if (x >= 0) and (y < 0):
				return 2*np.pi + gamma
			if (x < 0) and (y < 0):
				return np.pi + gamma

		z = Vector(z).norm; x = Vector(x).norm

		if not x.isperpendicular(z):
			warnings.warn("x is not perpendicular to z, result may be weird!",DeprecationWarning)

		y = z.cross(x)
		v = self.norm
		tth = v.betweenangle_rad(z)
		vx = v.dot(x); vy = v.dot(y)
		gamma = mygamma(vx,vy)
		return (tth, gamma)

	def tthgam_deg(self,**kw):
		return np.rad2deg(self.tthgam_rad(**kw))

	def rotate_by(self, q = None, axis = None, degree = None):
		if not q is None:
			if not isinstance(q, Quaternion):
				raise TypeError("parameter should be Quaternion from pyquaternion")
		else:
			q = Quaternion(axis = axis, degrees = degree)
		return Vector(q.rotate(self))