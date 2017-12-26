# -*- coding: utf-8 -*-
# @Date    : 2017-12-22 16:30:40
# @Author  : J. W. Huang (huangjasper@126.com)

# Public modules
import numpy as np
from scipy import linalg
import sys
import logging; logging.basicConfig(level = logging.INFO)
import itertools

# My modules
import LUT
from Vec import Vectors
EPS = 1e-5
EPS_sf = 1e-5

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


class LatticeParameter(object):
	'''
		Lattice Paramter
		Attributes:
			a,b,c,alpha,beta,gamma: six lattice parameters

		Methods:
			direct_matrix(): the direct matrix of lattice
			reciprocal_matrix(): the reciprocal matrix of lattice
			report(): show the informations of LatticeParamter
	'''

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

		logging.debug('lattice parameters = (%f %f %f %f %f %f)'%(self.a, self.b, self.c, np.rad2deg(self.alpha), np.rad2deg(self.beta), np.rad2deg(self.gamma)))
		if not self.isvalid():
			logging.error('lattice parameter is invalid!')
			raise ValueError

	def isvalid(self):
		if self.a <= 0 or self.b <= 0 or self.c <= 0:
			return False

		for angle in (self.alpha, self.beta, self.gamma):
			if angle >= np.pi or angle <= 0:
				return False
		return True

	def direct_matrix(self):
		if hasattr(self, '_direct_matrix'):
			return self._direct_matrix

		v1 = self.a*Vectors((1,0,0))
		v2 = self.b*Vectors((np.cos(self.gamma),np.sin(self.gamma),0))
		v3 = self.c*Vectors((np.cos(self.beta),(np.cos(self.alpha)-np.cos(self.beta)*np.cos(self.gamma))/np.sin(self.gamma))).norm
		m = np.vstack((v1,v2,v3))

		logging.debug('Direct matrix is\n %s'%(str(m)))
		self._direct_matrix = m
		return m

	def reciprocal_matrix(self):
		if hasattr(self, '_reciprocal_matrix'):
			return self._reciprocal_matrix

		def Calc_reciprocal_vectors(v):
			return np.vstack((np.cross(v[1],v[2]),np.cross(v[2],v[0]),np.cross(v[0],v[1])))/linalg.det(v)
		m = Calc_reciprocal_vectors(self.direct_matrix())

		logging.debug('Reciprocal matrix is\n %s'%(str(m)))
		self._reciprocal_matrix = m
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

class FracCoor(Vectors):
	def __new__(cls, input_array):
		if len(input_array) != 3:
			logging.error('Coordinates must be 3-dimensional!')
			raise ValueError
		if any((i < 0 or i > 1 for i in input_array)):
			logging.error('Fractional coordinates must be 0~1!')
			raise ValueError

		return np.asarray(input_array).view(cls)

	@property
	def x(self):
		return self[0]

	@property
	def y(self):
		return self[1]

	@property
	def z(self):
		return self[2]

class Element(object):
	def __init__(self, name):
		super(Element, self).__init__()

		self.name, self.mass, self.scatter_factor = LUT.dict_element[name]

class Lattice(object):
	def __init__(self, *, material = None, structure = None, args = None):
		super(Lattice, self).__init__()

		self.LP = None
		self.atoms = []

		if not material is None or (not structure is None and not args is None):
			if not material is None:
				# material is input
				structure, args = LUT.dict_material[material]
				latp, atomsymbols = Replace_placeholder(LUT.dict_structure[structure], args)
			else:
				# structure and args is input
				latp, atomsymbols = Replace_placeholder(LUT.dict_structure[structure], args)

			self.Add_latticeparameters(LatticeParameter(latp))
			for atomsymbol in atomsymbols[1:]:
				self.Add_atom(Atom(*atomsymbol))

			logging.info('init of lattice is done!')
		else:
			logging.info('No args input and empty lattice is initialized!')

	def Add_latticeparameters(self, latticeparam):
		if not type(latticeparam) is LatticeParameter:
			logging.error('argument should be type \'LatticeParameter\'')
			raise TypeError

		if not self.LP is None:
			logging.warning('the existing Lattice Parameters will be replaced!!')

		self.LP = latticeparam
		logging.info('Appending of lattice parameters is done!')

	def Add_atom(self, atom):
		if not type(atom) is Atom:
			logging.error('argument should be type \'Atom\'')
			raise TypeError

		self.atoms.append(atom)
		logging.debug('Appending of atom %s at %s is done!'% (atom.element.name, str(atom.coor)))

	def report(self):
		self.LP.report()
		for atom in self.atoms:
			print("atom %s at %s"%(atom.element.name, str(atom.coor)))

	def factor(self, hkls):
		'''
		 	Structure factor for one or a group of hkl
		'''
		if len(self.atoms) == 0:
			logging.error('No atoms in Lattice!')
			raise ValueError

		if isinstance(hkls, index):
			hkls = hkls[None,:]

		results = []
		for hkl in hkls:
			sf = 0
			for atom in self.atoms:
				sf0 = np.sum(atom.coor * hkl)
				sf += complex(atom.element.scatter_factor * np.exp(2j*np.pi*sf0))
			results.append(sf)
			logging.debug('The Structure Factor of %s is %r\n'%(str(hkl), sf))

		return np.array(results)

	def isextinct(self, hkls):
		'''
			Return True when this hkl is structurally extinct
		'''
		if isinstance(hkls, index):
			hkls = hkls[None,:]
		results = [np.absolute(self.factor(hkl)) < EPS_sf for hkl in hkls]
		return np.array(results)

	def d_spacing(self, hkls):
		'''
			return the d-spacing of one or a group of given hkls
			if hkl is extinct, None will be returned 
		'''
		
		results = []
		if isinstance(hkls, index):
			hkls = hkls[None,:] # add a dimensional to make it iterable

		for hkl in hkls:
			if self.isextinct(hkl):
				results.append(None)
			else:
				k_vector = Vectors(self.LP.reciprocal_matrix() * hkl)
				results.append(1 / k_vector.length)

		return np.array(results)


def Gen_hklfamilies(hklrange = (5,5,5), * , lattice = None):
	'''
		give an array of hklfamilies from given range of index
		NOTE:
		This is a GENERATOR
	'''
	for h in range(hklrange[0] + 1):
		for k in range(min(h, hklrange[1]) + 1):
			for l in range(min(k, hklrange[2]) + 1):
				if h == 0 and k == 0 and  l == 0:
					continue
				hkl = Familyindex((h,k,l))
				if (lattice is None) or (not lattice.isextinct(hkl)):
					yield hkl


def Gen_hkls(hklfamilies):
	'''
		give an array of hkl from given hklfamilies
		NOTE:
		This is a GENERATOR
	'''

	if isinstance(hklfamilies, Familyindex):
		hklfamilies = hklfamilies[None,:]

	for hklfamily in hklfamilies:
		for hkl in hklfamily.sons():
			yield hkl

class index(Vectors):
	'''
		index of lattice plane in a Xtal
	'''
	def __new__(cls, input_array):

		if len(input_array) != 3:
			logging.error('index must be 3-dimensional!')
			raise TypeError
		if any((np.fabs(e - np.floor(e)) > EPS for e in input_array)):
			logging.warning('Non-integer is input to hkl')

		return np.array(input_array).view(cls)


class Familyindex(index):
	'''
		index of a family of lattice plane in a Xtal
	'''
	def __new__(cls, input_array):
		return super().__new__(cls,input_array)

	def sons(self):
		'''
			give the sons of this family of lattice plane
		'''
		store = set()
		perms = list(itertools.permutations(self, 3))
		signs = list(itertools.product((1,-1),(1,-1),(1,-1)))
		for perm in perms:
			for sign in signs:
				i = [p*s for p,s in zip(perm, sign)]
				store.add(tuple(i))
		return [index(s) for s in store]

if __name__ == '__main__':
	l = Lattice(material = 'Cu')
	l.report()
