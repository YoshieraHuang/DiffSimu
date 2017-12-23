# -*- coding: utf-8 -*-
# @Date    : 2017-12-22 16:30:40
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
from scipy import linalg
import sys
import logging
import LUT
from Vec import Vectors

logging.basicConfig(level = logging.INFO)

class LatticeParameter(object):
	"""include lattice parameters"""

	def __init__(self, arr):
		super(LatticeParameter, self).__init__()
		if len(arr) == 6:
			self.a, self.b, self.c, self.alpha, self.beta, self.gamma = arr[0], arr[1], arr[2], np.deg2rad(arr[3]), np.deg2rad(arr[4]), np.deg2rad(arr[5])
		elif arr.shape == (3,3):
			va, vb, vc = arr
			va, vb, vc = Vectors(va), Vectors(vb), Vectors(vc)
			self.a, self.b, self.c = va.length, vb.length, vc.length
			self.alpha = vb.betweenangle_rad(vc)
			self.beta = vc.betweenangle_rad(va)
			self.gamma = va.betweenangle_rad(vb)

		logging.info('lattice parameters = (%f %f %f %f %f %f)'%(self.a, self.b, self.c, np.rad2deg(self.alpha), np.rad2deg(self.beta), np.rad2deg(self.gamma)))
		if not self.isvalid():
			logging.error('lattice parameter is invalid!')

	def isvalid(self):
		if self.a <= 0 or self.b <= 0 or self.c <= 0:
			return False

		for angle in (self.alpha, self.beta, self.gamma):
			if angle >= np.pi or angle <= 0:
				return False
		return True

	@property
	def direct_matrix(self):
		v1 = self.a*Vectors((1,0,0))
		v2 = self.b*Vectors((np.cos(self.gamma),np.sin(self.gamma),0))
		v3 = self.c*Vectors((np.cos(self.beta),(np.cos(self.alpha)-np.cos(self.beta)*np.cos(self.gamma))/np.sin(self.gamma))).norm
		m = np.vstack((v1,v2,v3))

		logging.info('Direct matrix is\t, (%s)'%(str(m)))
		return m

	@property
	def reciprocal_matrix(self):
		def Calc_reciprocal_vectors(v):
			return np.vstack((np.cross(v[1],v[2]),np.cross(v[2],v[0]),np.cross(v[0],v[1])))/linalg.det(v)
		m = Calc_reciprocal_vectors(self.direct_matrix)

		logging.info('Reciprocal matrix is\t, (%s)'%(str(m)))
		return m

	def report(self, fout = None):
		if fout != None:
			printfile = functools.partial(print, file = fout)
		else:
			printfile = print

		printfile('\ta = ', self.a, '\n\tb = ', self.b, '\n\tc = ',self.c)
		printfile('\talpha = ', np.rad2deg(self.alpha), '\n\tbeta = ', np.rad2deg(self.beta), '\n\tgamma = ', np.rad2deg(self.gamma))


class Atom(object):
	def __init__(self, element, atomcoor):
		super(Atom, self).__init__()

		if type(element) is Element:
			self.element = element
		else:
			self.element = Element(element)
		if type(atomcoor) is FracCoor:
			self.coor = atomcoor
		else:
			self.coor = FracCoor(atomcoor)

class FracCoor(object):
	def __init__(self, coors):
		super(FracCoor, self).__init__()

		if len(coors) == 3:
			pass
		else:
			logging.error('Coordinates are not 3-dimensional')
		
		for coor in coors:
			if coor < 0 or coor >1:
				logging.error('Fractional coordinates exceeds 0 ~ 1!')

		self.iter = np.array(coors)

	@property
	def x(self):
		return self.iter[0]

	@property
	def y(self):
		return self.iter[1]

	@property
	def z(self):
		return self.iter[2]

class Element(object):
	def __init__(self, name):
		super(Element, self).__init__()

		self.name, self.mass, self.scatter_factor = LUT.dict_element[name]

class Lattice(object):
	def __init__(self, *, material = None, structure = None, args = None):
		super(Lattice, self).__init__()

		self.LP = None
		self.atoms = []

		def Replace_placeholder(data, args, diction = {}, i = 0, top = True):
			if type(data) in (list,tuple):
				# if data is a list, recurse with the data
				results = []
				for datum in data:
					datum_done, diction, i = Replace_placeholder(datum, args, diction, i, False)
					results.append(datum_done)

				if top:
					return results
				else:
					return results, diction, i
			elif type(data) is str and data[0] == '$':
				# if data is a place holder, replace this
				if data in diction:
					if top:
						return diction[data]
					else:
						return diction[data], diction, i

				if type(args) in (list,tuple):
					diction[data] = args[i]
				else:
					diction[data] = args

				i += 1

				if top :
					return diction[data]
				else:
					return diction[data], diction , i
			else:

				if top:
					return data
				else:
					return data, diction, i


		if not material is None or (not structure is None and not args is None):
			if not material is None:
				# material is input
				structure, args = LUT.dict_material[material]
				latp, atomsymbols = Replace_placeholder(LUT.dict_structure[structure], args)
			else:
				# structure and args is input
				latp, atomsymbols = Replace_placeholder(LUT.dict_structure[structure], args)

			print(latp, atomsymbols)
			self.Add_latticeparameters(LatticeParameter(latp))
			for atomsymbol in atomsymbols[1:]:
				self.Add_atom(Atom(*atomsymbol))

			logging.info('init of lattice is done!')
		else:
			logging.info('No args input and empty lattice is initialized!')

	def Add_latticeparameters(self, latticeparam):
		if not type(latticeparam) is LatticeParameter:
			logging.error('argument should be type \'LatticeParameter\'')

		if not self.LP is None:
			logging.warning('the existing Lattice Parameters will be replaced!!')
		self.LP = latticeparam
		logging.info('Appending of lattice parameters is done!')

	def Add_atom(self, atom):
		if not type(atom) is Atom:
			logging.error('argument should be type \'Atom\'')

		self.atoms.append(atom)
		logging.info('Appending of atom %s at %s is done!'% (atom.element.name, str(atom.coor.iter)))

	def report(self):
		self.LP.report()
		for atom in self.atoms:
			print("atom %s at %s"%(atom.element.name, str(atom.coor.iter)))

