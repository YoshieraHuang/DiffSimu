# -*- coding: utf-8 -*-
# @Date    : 2017-12-26 19:15:18
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging; logging.basicConfig(level = logging.DEBUG)
import itertools


from pyquaternion import Quaternion
import lattice as LTTC
from Vec import Vector
EPS = 1e-5
np.set_printoptions(precision = 5, suppress = True)

'''
	CLASS DEFINITION
'''

class SingleXtal(object):
	'''
		class for single crystal
		Attributes:
			lattice: lattice information
			indexx, indexy, indexz: index at the direction in x-, y-, z-axis of lab frame;
			vecx, vecy, vecz: Vector in crystal at the direction in x-, y- and z-axis of lab frame;
				z is reverse-parallel to incident direction of x-ray; y is horizontal and x is verticle
				These Vectors are able to be set by value, which, however, are not recommanded
			reciprocal_matrix: reciprocal matrix of this single crystal, including the information of orientation 
		Method:
			check_orientation(): check whether the orientation of single crystal is decided.
			direction(index): the direction of given index in lab coordinates

	'''
	def __init__(self, lttc, *, x = None, y = None, z = None):
		super().__init__()

		self.lattice = lttc
		if (not x is None) and (not isinstance(x, LTTC.index)):
			self.indexx = LTTC.index(x)
		else:
			self.indexx = x

		if (not y is None) and (not isinstance(y, LTTC.index)):
			self.indexy = LTTC.index(y)
		else:
			self.indexy = y

		if (not z is None) and (not isinstance(z, LTTC.index)):
			self.indexz = LTTC.index(z)
		else:
			self.indexz = z

		self.orientation_decided = False
		self.reciprocal_matrix = self.lattice.reciprocal_matrix

		if self.check_orientation():
			logging.info('orientation is decided!')
			self.Cal_rcp_matrix()
		else:
			logging.info('orientation is not decided!')

	@property
	def indexx(self):
		return self._indexx

	@indexx.setter
	def indexx(self, v):
		if v is None:
			self._indexx = None
			return
		if not isinstance(v, LTTC.index):
			v = index(v)
		self._indexx = v

	@property
	def indexy(self):
		return self._indexy

	@indexy.setter
	def indexy(self,v):
		if v is None:
			self._indexy = None
			return
		if not isinstance(v, LTTC.index):
			v = index(v)
		self._indexy = v

	@property
	def indexz(self):
		return self._indexz

	@indexz.setter
	def indexz(self,v):
		if v is None:
			self._indexz = None
			return
		if not isinstance(v, LTTC.index):
			v = index(v)
		self._indexz = v

	@property
	def vecx(self):
		if hasattr(self, '_vecx'):
			return self._vecx
		else:
			if self.indexx is None:
				return None
			self._vecx = self.lattice.Vec_in_lattice(self.indexx)
			return self._vecx

	@vecx.setter
	def vecx(self, v):
		logging.warning('vecx is changed!')
		self._vecx = v

	@property
	def vecy(self):
		if hasattr(self,'_vecy'):
			return self._vecy
		else:
			if self.indexy is None:
				return None
			self._vecy = self.lattice.Vec_in_lattice(self.indexy)
			return self._vecy

	@vecy.setter
	def vecy(self, v):
		logging.warning('vecy is changed!')
		self._vecy = v

	@property
	def vecz(self):
		if hasattr(self,'_vecz'):
			return self._vecz
		else:
			if self.indexz is None:
				return None
			self._vecz = self.lattice.Vec_in_lattice(self.indexz)
			return self._vecz

	@vecz.setter
	def vecz(self, v):
		logging.warning('vecz is changed!')
		self._vecz = v

	def check_orientation(self):
		'''
			make sure that the orientation of single crystal is decided.
			if decided, unknown direction will be calculated and True will be returned. if not, return False
		'''
		if self.orientation_decided:
			return True 

		vecs = (self.vecx, self.vecy, self.vecz)
		number_of_None = 0
		for vec in vecs:
			if vec is None:
				number_of_None += 1
		logging.debug("number of None is %d"%(number_of_None))
		if number_of_None in (2,3):
			self.orientation_decided = False
			return False

		if number_of_None == 0:
			for pair in itertools.permutation(vecs,2):
				if not pair[0].isperpendicular(pair[1]):
					self.orientation_decided = False
					logging.error('indexx, indexy, indexz is not perpendicular!')
					return False

			v3 = vecs[0].cross(vecs[1])
			if vecs[2].isparallel(v3) != 1:
				logging.warning('three Vector may not be right-handed!')
				return False

			return True

		if number_of_None == 1:
			i_None, i1, i2 = -1,-1,-1
			for i,vec in enumerate(vecs):
				if vec is None:
					i_None = i
				elif i1 == -1:
					i1 = i
				else:
					i2 = i

			print(vecs)
			vec_unknown = vecs[i1].cross(vecs[i2])
			if i_None == 0:
				self.vecx = vec_unknown
			elif i_None == 1:
				self.vecy = vec_unknown
			else:
				self.vecz = vec_unknown

			if not vecs[i1].isperpendicular(vecs[i2]):
				logging.warning('Existing two Vector are not perpendicular, one Vector is adjusted!')
				vec_new = vecs[i2].cross(vec_unknown)
				if i1 == 0:
					self.vecx = vec_new
				elif i1 == 1:
					self.vecy = vec_new
				else:
					self.vecz = vec_new
				self.orientation_decided = True
			else:
				self.orientation_decided = True

			return True

	def Cal_rcp_matrix(self):
		qua = Guess_rotation((Vector((0,0,1)), Vector((1,0,0))), (self.vecz, self.vecx))
		self.reciprocal_matrix =  Rotate_vectors_by_qua(self.lattice.reciprocal_matrix, qua.inverse)
		logging.debug('reciprocal matrix is\n %s'%(self.reciprocal_matrix))

	def direction(self, hkls):
		return self.lattice.Vec_in_lattice(hkls, rcp_matrix = self.reciprocal_matrix)


'''
	FUNCTION DEFINITION
'''

def Rotate_vectors_by_qua(Vector, qua):
	return np.dot(Vector, qua.rotation_matrix.T)

def Guess_rotation(vecBefore, vecAfter):
	'''
		Give rotation quaternion by two changes of Vector
		Arguments:
			vecBefore: array contains two Vector before the rotation
			vecAfter: array contains two Vector after the rotation

		Returns:
			Quaternion
	'''
	logging.debug('vecBefore is %s, vecAfter is %s'%(str(vecBefore), str(vecAfter)))
	vec1Before, vec1After = vecBefore[0], vecAfter[0]
	if (Vector(vec1After - vec1Before).length < EPS):
		qua1 = Quaternion()
	else:
		axis1, angle1 = vec1Before.cross(vec1After), vec1Before.betweenangle_rad(vec1After)
		qua1 = Quaternion(axis = axis1, radians = angle1).unit

	logging.debug('axis1 = %s, angle1 = %f'%(str(qua1.axis), qua1.degrees))

	vec2AfterR1, vec2After = Rotate_vectors_by_qua(vecBefore[1], qua1), vecAfter[1]
	if (Vector(vec2AfterR1 - vec2After).length < EPS):
		qua2 = Quaternion()
	else:
		axis2 = vec1After
		angle2 = vec2After.cross(axis2).betweenangle_rad(vec2AfterR1.cross(axis2))
		if vec2AfterR1.cross(vec2After).dot(axis2) < 0:
			angle2 = -angle2
		qua2 = Quaternion(axis = axis2, radians = angle2).unit

	logging.debug('axis2 = %s, angle2 = %f'%(str(qua2.axis), qua2.degrees))

	qua = qua2 * qua1
	logging.debug('axis = %s, angle2 = %f'%(str(qua.axis), qua.degrees))
	logging.debug('Residual:\n %s'%(str(Rotate_vectors_by_qua(vecBefore,qua) - vecAfter)))
	return qua

'''
	TEST
'''
if __name__ == '__main__':
	l = LTTC.Lattice(material = 'Cu')
	sx = SingleXtal(l, z = LTTC.index((1,1,0)), x = LTTC.index((0,0,1)))
	hkls = LTTC.index((0,1,0)), LTTC.index((1,1,0))
	print(sx.direction(hkls))