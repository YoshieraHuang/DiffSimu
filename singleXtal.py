# -*- coding: utf-8 -*-
# @Date    : 2017-12-26 19:15:18
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging; logging.basicConfig(level = logging.INFO)
import itertools
from scipy import linalg
import copy

from pyquaternion import Quaternion
import lattice as LTTC
from Vec import Vector
import myfunctools as FT
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
			rcp_matrix: reciprocal matrix of this single crystal, including the information of orientation 
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
		self.rcp_matrix = self.lattice.rcp_matrix

		if self.check_orientation():
			logging.info('orientation is decided!')
			# self.Guess_rotation((Vector(0,0,1), Vector(1,0,0)), (self.vecz, self.vecx))
			self.Guess_rotation((self.vecz, self.vecx), (Vector(0,0,1), Vector(1,0,0)))
			self.Calc_rcp_matrix()
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
			self._vecx = self.vec_in_rcp(self.indexx).norm
			return self._vecx

	@vecx.setter
	def vecx(self, v):
		logging.debug('vecx is set as %s'%(str(v)))
		self._vecx = v

	@property
	def vecy(self):
		if hasattr(self,'_vecy'):
			return self._vecy
		else:
			if self.indexy is None:
				return None
			self._vecy = self.vec_in_rcp(self.indexy).norm
			return self._vecy

	@vecy.setter
	def vecy(self, v):
		logging.debug('vecy is set as %s'%(str(v)))
		self._vecy = v

	@property
	def vecz(self):
		if hasattr(self,'_vecz'):
			return self._vecz
		else:
			if self.indexz is None:
				return None
			self._vecz = self.vec_in_rcp(self.indexz).norm
			return self._vecz

	@vecz.setter
	def vecz(self, v):
		logging.debug('vecz is set as %s'%(str(v)))
		self._vecz = v

	@property
	def index(self):
		return (self.indexx, self.indexy, self.indexz)

	@property
	def vecs(self):
		return (self.vecx, self.vecy, self.vecz)

	def check_orientation(self):
		'''
			make sure that the orientation of single crystal is decided.
			if decided, unknown direction will be calculated and True will be returned. if not, return False
		'''

		vecs = self.vecs
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

			vec_unknown = vecs[i1].cross(vecs[i2]).norm

			if i_None == 0:
				self.vecx = vec_unknown
			elif i_None == 1:
				self.vecy = - vec_unknown
			else:
				self.vecz = vec_unknown

			if not vecs[i1].isperpendicular(vecs[i2]):
				logging.warning('Vector %s and %s are not perpendicular, one Vector is adjusted!'%(str(vecs[i1]), str(vecs[i2])))
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

	def Calc_rcp_matrix(self):
		self.rcp_matrix = Rotate_vectors_by_qua(self.lattice.rcp_matrix, self.R)
		logging.debug('reciprocal matrix is\n %s'%(self.rcp_matrix))

	@property
	def rcp_matrix(self):
		if not hasattr(self, '_rcp_matrix'):
			self._rcp_matrix = Calc_reciprocal_vectors(self.direct_matrix)

		return self._rcp_matrix

	@rcp_matrix.setter
	def rcp_matrix(self, m):
		self._rcp_matrix = m
		logging.debug('reciprocal matrix is set as %s'%(str(m)))
		if hasattr(self, '_direct_matrix'):
			del self._direct_matrix

	@property
	def direct_matrix(self):
		if not hasattr(self, '_direct_matrix'):
			self._direct_matrix = Calc_reciprocal_vectors(self.rcp_matrix)

		return self._direct_matrix

	@direct_matrix.setter
	def direct_matrix(self, m):
		self._direct_matrix = m
		logging.debug('direct matrix is set as %s'%(str(m)))
		if hasattr(self, '_rcp_matrix'):
			del self._rcp_matrix

	def Gen_vec_in_rcp(self, hkls):
		return self.lattice.Gen_vec_in_rcp(hkls, rcp_matrix = self.rcp_matrix)

	def vec_in_rcp(self, hkl):
		return self.lattice.vec_in_rcp(hkl, rcp_matrix = self.rcp_matrix)

	def tthgam_rad(self, vec, z = None, x = None):
		return vec.tthgam_rad(z = z, x = x)

	def tthgam_deg(self, vec, *args):
		return vec.tthgam_deg(*args)

	def strain1d(self, axis = (1,0,0), index = None, ratio = 0):
		if not index is None:
			axis = self.vec_in_rcp(index).norm
		else:
			axis = Vector(axis).norm

		vecs = []
		for vec in self.direct_matrix:
			d = vec.dot(axis)
			vec_ = vec + axis * ratio * d
			vecs.append(vec_)

		self.lattice.Add_latticeparameters(LTTC.LatticeParameter(np.array(vecs)))
		self.direct_matrix = vecs

	def Guess_rotation(self, vsbefore, vsafter):
		vbefore1, vafter1 = vsbefore[0].norm, vsafter[0].norm
		if ((vafter1 - vbefore1).length < EPS):
			axis1, angle1 = Vector((1,0,0)), 0
		else:
			axis1, angle1 = vbefore1.cross(vafter1), vbefore1.betweenangle_rad(vafter1)
		qua1 = Quaternion(axis = axis1, angle = angle1).unit
	
		logging.debug("axis1 = %s, angle1 = %f"%(str(qua1.axis), qua1.degrees))
		logging.debug("Residual = %s"%(str(vbefore1.rotate_by(qua1) - vafter1)))
	
		vbefore2, vafter2 = vsbefore[1].norm, vsafter[1].norm
		vafterq1 = vbefore2.rotate_by(qua1)
		if ((vafterq1 - vafter2).length < EPS):
			axis2, angle2 = Vector((1,0,0)), 0
		else:
			axis2 = vafter1
			angle2 = vafter2.cross(axis2).betweenangle_rad(vafterq1.cross(axis2))
			if vafterq1.cross(vafter2).dot(axis2) < 0:
				angle2 = -angle2
		qua2 = Quaternion(axis = axis2, angle = angle2).unit
	
		logging.debug("axis2 = %s, angle2 = %f"%(str(qua2.axis), qua2.degrees))
		logging.debug("Residual = %s"%(str(vafterq1.rotate_by(qua2) - vafter2)))
	
		qua = qua2 * qua1
	
		logging.debug("axis = %s, angle = %f"%(str(qua.axis), qua.degrees))
		
		self.R = qua

	def rotate_by(self, q = None, axis = None, degree = None):
		if q is None:
			if axis is None or degree is None:
				raise ValueError('Empty Parameters')
			q = Quaternion(axis = axis, degrees = degree)
		logging.debug('rotation: axis = %s, degree = %f'%(str(q.axis), q.degrees))
		self.vecx = self.vecx.rotate_by(q)
		self.vecy = self.vecy.rotate_by(q)
		self.vecz = self.vecz.rotate_by(q)
		self.R = q * self.R
		self.rcp_matrix = Rotate_vectors_by_qua(self.rcp_matrix, q)

	def Calc_rcp_space(self, hklrange, *, hkls = None):
		if not hkls:
			hkls = list(LTTC.Gen_hkls(hklrange = hklrange, lattice = self.lattice))
		self.rcp_space_hkls = hkls
		self.rcp_space_vec = list(self.Gen_vec_in_rcp(hkls))
		self.rcp_space_num = len(hkls)

	def Save_rcp_space(self, filename):
		if not hasattr(self, 'rcp_space_num'):
			raise ValueError('No rcp_space')

		with open(filename, 'w') as f:
			f.writelines('%d\n'%(self.rcp_space_num))
			f.writelines('h k l kx ky kz\n')
			for hkl, vec in zip(self.rcp_space_hkls, self.rcp_space_vec):
				f.writelines('%s %f %f %f\n'%(hkl.str, vec[0], vec[1], vec[2]))

	def Project_rcp_space(self, proj_vec = Vector(0,0,1), vecx = Vector(1,0,0), origin = Vector(0,0,0)):
		if not hasattr(self, 'rcp_space_num'):
			raise ValueError('No rcp_space')

		if not isinstance(proj_vec, Vector):
			proj_vec = Vector(proj_vec)

		if not isinstance(origin, Vector):
			origin = Vector(origin)

		if not isinstance(vecx, Vector):
			vecx = Vector(vecx)

		self.project_rcp_space_vec = list(self.Gen_project_vec(self.rcp_space_vec, proj_vec = proj_vec.norm, vecx = vecx.norm, origin = origin))

	def Gen_project_vec(self, vecs, proj_vec = Vector(0,0,1), vecx = Vector(1,0,0), origin = Vector(0,0,0)):
		for vec in vecs:
			yield self.project_vec(vec, proj_vec = proj_vec.norm, vecx = vecx.norm, origin = origin)

	def project_vec(self, vec, proj_vec = Vector(0,0,1), vecx = Vector(1,0,0), origin = Vector(0,0,0)):
		proj_vec = proj_vec.norm
		vecx = vecx.norm
		dis = (vec - origin).dot(proj_vec) ## distance from initial point to projection plane
		projected = vec - origin - dis*proj_vec ## projection vector on projection plane in 3D 
		vecy = proj_vec.cross(vecx) ## y axis in projection plane
		proj_x, proj_y = projected.dot(vecx), projected.dot(vecy) ## coordinates in projection coordinates system
		return (proj_x, proj_y, dis)

	def Save_proj_rcp_space(self, filename):
		if not hasattr(self, 'project_rcp_space_vec'):
			raise ValueError('No projected reciprocal space')

		with open(filename,'w') as f:
			f.writelines('%d\n'%(self.rcp_space_num))
			f.writelines('h k l vx vy distance\n')
			for hkl, vec in zip(self.rcp_space_hkls, self.project_rcp_space_vec):
				f.writelines('%s %f %f %f\n'%(hkl.str, vec[0], vec[1], vec[2]))

	def copy(self):
		return copy.copy(self)

'''
	FUNCTION DEFINITION
'''
def Calc_reciprocal_vectors(v):
	return np.vstack((np.cross(v[1],v[2]),np.cross(v[2],v[0]),np.cross(v[0],v[1])))/linalg.det(v)

def Rotate_vectors_by_qua(vs, qua):
	return np.vstack(map(qua.rotate, vs))

'''
	TEST
'''
if __name__ == '__main__':
	l = LTTC.Lattice(material = 'Cu')
	sx = SingleXtal(l, z = (0,0,1), x = (1,0,0))
	sx.Calc_rcp_space((5,5,5))
	sx.Save_rcp_space('rcp2.dat')
	sx.Project_rcp_space(proj_vec = (1,1,1), vecx = (1,-1,0))
	sx.Save_proj_rcp_space('rcp_proj.dat')
	# print(sx.vecs,sx.R,sx.rcp_matrix)
	# sx.strain1d(direction = (1,1,1), ratio = -0.04)
	# sx.lattice.LP.report()
	# sx.strain1d(direction = (-1,1,0), ratio = 0.02)
	# sx.strain1d(direction = (-1,-1,2), ratio = 0.02)
	# sx.lattice.LP.report()
	# print(sx.Vec_in_sx((LTTC.index((1,1,1)), LTTC.index((1,-1,-1)), LTTC.index((-1,-1,-1)))))