# -*- coding: utf-8 -*-
# @Date    : 2017-12-22 16:30:40
# @Author  : J. W. Huang (huangjasper@126.com)

# Public modules
import numpy as np
from scipy import linalg
import sys
import logging; logging.basicConfig(level = logging.INFO)
import itertools
import re

# My modules
import LUT
from Vec import Vector
import myfunctools as FT

# Parameters
EPS = 1e-5
EPS_sf = 1e-5
re_elesym = re.compile(r'^\s*([A-Za-z]+)([+-]*|\d[+-]+)\s*$') # Regex for element name

'''
	CLASS DEFINITION
'''

class LatticeParameter(object):
	'''
		Lattice Paramter
		Attributes:
			a,b,c,alpha,beta,gamma: six lattice parameters
			direct_matrix: the direct matrix of lattice
			reciprocal_matrix: the reciprocal matrix of lattice

		Methods:
			report(): show the informations of LatticeParamter
	'''

	def __init__(self, arr):
		super(LatticeParameter, self).__init__()
		if len(arr) == 6:
			self.a, self.b, self.c, self.alpha, self.beta, self.gamma = arr[0], arr[1], arr[2], np.deg2rad(arr[3]), np.deg2rad(arr[4]), np.deg2rad(arr[5])
		elif arr.shape == (3,3):
			va, vb, vc = arr
			va, vb, vc = Vector(va), Vector(vb), Vector(vc)
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

	@property
	def direct_matrix(self):
		if hasattr(self, '_direct_matrix'):
			return self._direct_matrix

		v1 = self.a*Vector((1,0,0))
		v2 = self.b*Vector((np.cos(self.gamma),np.sin(self.gamma),0))
		v3 = self.c*Vector((np.cos(self.beta),(np.cos(self.alpha)-np.cos(self.beta)*np.cos(self.gamma))/np.sin(self.gamma))).norm
		m = np.vstack((v1,v2,v3))

		logging.debug('Direct matrix is\n %s'%(str(m)))
		self._direct_matrix = m
		return m

	@property
	def reciprocal_matrix(self):
		if hasattr(self, '_reciprocal_matrix'):
			return self._reciprocal_matrix

		def Calc_reciprocal_vectors(v):
			return np.vstack((np.cross(v[1],v[2]),np.cross(v[2],v[0]),np.cross(v[0],v[1])))/linalg.det(v)
		m = Calc_reciprocal_vectors(self.direct_matrix)

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






class FracCoor(Vector):
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

		self.name = self.reg(name)
		self.sc_factor_coeff = [np.array(coeff) for coeff in LUT.dict_element[self.name]]

	def reg(self, name):
		re_match = re_elesym.match(name)
		if re_match is None:
			raise ValueError('Strange element name')
		reg_name = re_match.group(1).capitalize()
		if re_match.group(2) in ('+', '-'):
			reg_name = reg_name + '1' + re_match.group(2)
		else:
			reg_name = reg_name + re_match.group(2)
		return reg_name

	def sc_factor(self, tth, wl):
		a, b, c = self.sc_factor_coeff
		return np.sum(a * np.exp( b * ( - (np.sin(tth/2)/wl)**2 ) ) ) + c






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

	def strct_factor(self, hkls):
		'''
		 	Structure factor for one or a group of hkl
		 	Note: scattering factors of all elements is considered as 1
		 		if real scattering factors should be considered, please refer to function sc_factor() 
		'''
		if len(self.atoms) == 0:
			logging.error('No atoms in Lattice!')
			raise ValueError

		if isinstance(hkls, index):
			hkls = hkls[None,:]

		for hkl in hkls:
			sf = 0
			for atom in self.atoms:
				kr = float(np.sum(atom.coor * hkl))
				sf0 = np.exp(2j*np.pi*kr)
				sf += sf0

			if np.abs(sf) < EPS_sf:
				sf = None
			
			logging.debug('The Structure Factor of %s is %r\n'%(str(hkl), sf))
			yield sf
	
	def sc_factor(self, hkl, tth, wl):
		'''
			Only 1 hkl can be input because always called for calculation of 1 hkl
		'''
		if len(self.atoms) == 0:
			logging.error('No atoms in Lattice!')
			raise ValueError

		sf = 0
		for atom in self.atoms:
			kr = float(np.sum(atom.coor * hkl))
			sf0 = atom.element.sc_factor(tth, wl)* np.exp(2j*np.pi*kr)
			sf += sf0

		if np.abs(sf) < EPS_sf:
			sf = None

		logging.debug('The Structure Factor of %s is %r\n'%(str(hkl), sf))

		return sf 


	@property
	def reciprocal_matrix(self):
		return self.LP.reciprocal_matrix

	def isextinct(self, hkls):
		'''
			Return True when this hkl is structurally extinct
		'''
		for sta in self.strct_factor(hkls):
			if sta is None:
				yield True
			else:
				yield False

	def d_spacing(self, hkls):
		'''
			return the d-spacing of one or a group of given hkls
			if hkl is extinct, None will be returned 
		'''

		for vec in self.vec_in_lattice(hkls):
			yield 1/vec.length

	def D_spacing(self, hkls):
		return FT.tolist(self.d_spacing(hkls))

	def vec_in_lattice(self , hkls, rcp_matrix = None):
		'''
			transform index to Vector in lattice orthogonal frame 
		'''

		if isinstance(hkls, index):
			hkls = hkls[None,:]

		if rcp_matrix is None:
			rcp_matrix = self.reciprocal_matrix

		for hkl in hkls:
			yield Vector(np.dot(hkl,rcp_matrix))

	def Vec_in_lattice(self, hkls, rcp_matrix = None):
		return FT.tolist(self.vec_in_lattice(hkls, rcp_matrix))

	def isperpendicular(self, hkl1, hkl2):
		'''
			if these two hkls are perpendicular, return True
			else return False 
		'''

		vec1, vec2 = self.Vec_in_lattice(hkl1), self.Vec_in_lattice(hkl2)
		return vec1.isperpendicular(vec2)

	def isparallel(self, hkl1, hkl2):
		'''
			if these two hkls are parallel and have same direction, return 1;
			if reverse-parallel, return 2;
			else return 0
		'''
		vec1, vec2 = self.Vec_in_lattice(hkl1), self.Vec_in_lattice(hkl2)
		return vec1.isparallel(vec2)


class index(Vector):
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

	@property
	def str(self):
		return str(self[0]) + str(self[1]) + str(self[2])


class Familyindex(index):
	'''
		index of a family of lattice plane in a Xtal
	'''

	def sons(self):
		'''
			give the sons of this family of lattice plane
		'''
		if hasattr(self, '_sons'):
			return self._sons

		store = set()
		perms = list(itertools.permutations(self, 3))
		signs = list(itertools.product((1,-1),(1,-1),(1,-1)))
		for perm in perms:
			for sign in signs:
				i = [p*s for p,s in zip(perm, sign)]
				store.add(tuple(i))
		self._sons = [index(s) for s in store]
		return self._sons

	@property
	def multiplicity(self):
		return len(self.sons())


'''
	FUNCTION DEFINITION
'''


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

def Gen_hklfamilies(hklrange = (10,10,10), * , lattice = None):
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
				if (lattice is None) or (not all(lattice.isextinct(hkl))):
					yield hkl


def Gen_hkls(hklfamilies):
	'''
		give an array of hkl from given hklfamilies
		NOTE:
		This is a GENERATOR
	'''

	for hklfamily in hklfamilies:
		for hkl in hklfamily.sons():
			yield hkl



'''
	TEST
'''

if __name__ == '__main__':
	# l = Lattice(material = 'Cu')
	# hkls = FT.tolist(Gen_hklfamilies((3,3,3), lattice = l))
	# ds = l.d_spacing(hkls)
	# fs = l.strct_factor(hkls)
	# for hkl,d,f in zip(hkls,ds,fs):
	# 	print(hkl,hkl.multiplicity, d,f)

	# ele = Element('Cu')
	# print(ele.sc_factor(np.deg2rad(30), 0.5))
	print(index((0,0,-1)).str)