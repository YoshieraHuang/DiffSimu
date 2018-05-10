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
import types
from collections import Iterable

# My modules
import LUT
from Vec import Vector
import myfunctools as FT

# Parameters
EPS = 1e-5
EPS_sf = 1e-7
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
			rcp_matrix: the reciprocal matrix of lattice

		Methods:
			report(): show the informations of LatticeParamter
	'''

	def __init__(self, arr):
		super(LatticeParameter, self).__init__()
		logging.debug('input arr is %s'%(str(arr)))
		if len(arr) == 6:
			self.a, self.b, self.c, self.alpha, self.beta, self.gamma = arr[0], arr[1], arr[2], np.deg2rad(arr[3]), np.deg2rad(arr[4]), np.deg2rad(arr[5])
		elif np.array(arr).shape == (3,3):
			logging.debug('input matrix is %s'%(str(arr)))
			va, vb, vc = arr
			va, vb, vc = Vector(va), Vector(vb), Vector(vc)
			self.a, self.b, self.c = va.length, vb.length, vc.length
			self.alpha = vb.betweenangle_rad(vc)
			self.beta = vc.betweenangle_rad(va)
			self.gamma = va.betweenangle_rad(vb)
		else:
			raise ValueError('Unknown Parameters!')

		logging.debug('lattice parameters = (%f %f %f %f %f %f)'%(self.a, self.b, self.c, np.rad2deg(self.alpha), np.rad2deg(self.beta), np.rad2deg(self.gamma)))
		if not self.isvalid():
			raise Error('lattice parameter is invalid!')

	def isvalid(self):
		if self.a <= 0 or self.b <= 0 or self.c <= 0:
			return False

		for angle in (self.alpha, self.beta, self.gamma):
			if angle >= np.pi or angle <= 0:
				return False
		return True

	def isorthogonal(self):
		deg90 = np.pi/2
		EPS_deg = 1e-5
		if not FT.equal(self.alpha , deg90, EPS_deg):
			return False
		if not FT.equal(self.beta , deg90, EPS_deg):
			return False
		if not FT.equal(self.gamma , deg90, EPS_deg):
			return False
		return True

	def iscubic(self):
		EPS_a = 1e-4
		if not self.isorthogonal():
			return False

		if not FT.equal(self.a , self.b, EPS_a):
			return False

		if not FT.equal(self.b , self.c, EPS_a):
			return False

		return True

	def ishcp(self):
		deg90 = np.pi/2
		deg120 = np.deg2rad(120)
		EPS_deg = 1e-5
		EPS_a = 1e-4
		if not FT.equal(self.alpha , deg90, EPS_deg):
			return False
		if not FT.equal(self.beta , deg90, EPS_deg):
			return False
		if not FT.equal(self.gamma , deg120, EPS_deg):
			return False

		if not FT.equal(self.a , self.b, EPS_a):
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
	def rcp_matrix(self):
		if hasattr(self, '_rcp_matrix'):
			return self._rcp_matrix

		def Calc_rcp_vectors(v):
			return np.vstack((np.cross(v[1],v[2]),np.cross(v[2],v[0]),np.cross(v[0],v[1])))/linalg.det(v)
		m = Calc_rcp_vectors(self.direct_matrix)

		logging.debug('Reciprocal matrix is\n %s'%(str(m)))
		self._rcp_matrix = m
		return m

	def show(self, fout = None):
		if fout != None:
			printfile = functools.partial(print, file = fout)
		else:
			printfile = print

		printfile('\ta = ', self.a, '\n\tb = ', self.b, '\n\tc = ',self.c)
		printfile('\talpha = ', np.rad2deg(self.alpha), '\n\tbeta  = ', np.rad2deg(self.beta), '\n\tgamma = ', np.rad2deg(self.gamma))





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
	def __new__(cls, coor):
		if len(coor) != 3:
			raise ValueError('Coordinates must be 3-dimensional!')

		if any((i < 0 or i > 1 for i in coor)):
			raise ValueError('Fractional coordinates must be 0~1!')

		return np.asarray(coor).view(cls)

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
	def __init__(self, material = None, structure = None, args = None):
		super(Lattice, self).__init__()

		self.LP = None
		self.atoms = []

		if not material is None or (not structure is None and not args is None):
			if not material is None:
				# material is input
				structure, args = LUT.dict_material[material]
				latp, atomsymbols = Replace_placeholder(LUT.dict_structure[structure], args, diction = {})
			else:
				# structure and args is input
				latp, atomsymbols = Replace_placeholder(LUT.dict_structure[structure], args, diction = {})

			self.Add_latticeparameters(LatticeParameter(latp))
			for atomsymbol in atomsymbols[1:]:
				self.Add_atom(Atom(*atomsymbol))

			logging.info('init of lattice is done!')
		else:
			logging.info('No args input and empty lattice is initialized!')

	def Add_latticeparameters(self, latticeparam):
		if not type(latticeparam) is LatticeParameter:
			raise TypeError('argument should be type \'LatticeParameter\'')

		if not self.LP is None:
			logging.warning('the existing Lattice Parameters will be replaced!!')

		self.LP = latticeparam
		logging.info('Appending of lattice parameters is done!')

	def Add_atom(self, atom):

		if isinstance(atom, Atom):
			self.atoms.append(atom)
			logging.debug('Appending of atom %s at %s is done!'% (atom.element.name, str(atom.coor)))
		elif isinstance(atom, Iterable):
			for a in atom:
				if not isinstance(a, Atom):
					raise TypeError('argument should be type \'Atom\'')

				self.atoms.append(a)
				logging.debug('Appending of atom %s at %s is done!'% (a.element.name, str(a.coor)))

	def Atoms_from_structure(self, structure, element):

		atomsymbols = Replace_placeholder(LUT.dict_structure[structure][1], (element), diction = {})
		atoms = []
		for atomsymbol in atomsymbols[1:]:
			atoms.append(Atom(*atomsymbol))
		return atoms

	def show(self):
		self.LP.report()
		for atom in self.atoms:
			print("atom %s at %s"%(atom.element.name, str(atom.coor)))

	def Gen_strct_factor(self, hkls):
		'''
		 	Generator
		 	Structure factor for a group of hkl
		 	Note: scattering factors of all elements is considered as 1
		 		if real scattering factors should be considered, please refer to function sc_factor() 
		'''
		if isinstance(hkls, index):
			hkls = hkls[None,:]

		for hkl in hkls:
			yield self.strct_factor(hkl)

	def strct_factor(self, hkl):
		'''
			Structure factor for one hkl
			If hkl is Familyindex(), multiplicity is considered
			Note: scattering factors of all elements is considered as 1
		 		if real scattering factors should be considered, please refer to function sc_factor()
		'''
		if not isinstance(hkl, index):
			raise ValueError('Not index')
		if len(self.atoms) == 0:
			logging.error('No atoms in Lattice!')
			raise ValueError

		sf = 0
		for atom in self.atoms:
			kr = float(np.sum(atom.coor * hkl))
			sf0 = atom.element.sc_factor_coeff[2] * np.exp(2j*np.pi*kr)
			sf += sf0

		if np.abs(sf) < EPS_sf:
			logging.debug('The Structure Factor of %s is %r, too small\n'%(str(hkl), sf))
			return None
		
		if isinstance(hkl, Familyindex):
			sf = hkl.multiplicity*sf

		logging.debug('The Structure Factor of %s is %r\n'%(str(hkl), sf))
		return sf	
	
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
	def rcp_matrix(self):
		return self.LP.rcp_matrix

	def isextinct(self, hkl):
		'''
			Return True when this hkl is structurally extinct
		'''
		sta = self.strct_factor(hkl)
		return True if sta is None else False

	def Gen_d_spacing(self, hkls):
		'''
			Generator
			return the d-spacing of one or a group of given hkls
			if hkl is extinct, None will be returned 
		'''

		for vec in self.Gen_vec_in_rcp(hkls):
			yield 1/vec.length

	def d_spacing(self, hkl):
		'''
			only for one index
		'''
		return 1/(self.vec_in_rcp(hkl).length)

	def Gen_vec_in_rcp(self , hkls, rcp_matrix = None):
		'''
			transform index to Vector in lattice orthogonal frame 
		'''
		if isinstance(hkls, index):
			hkls = hkls[None,:]

		for hkl in hkls:
			yield self.vec_in_rcp(hkl, rcp_matrix = rcp_matrix)

	def vec_in_rcp(self, hkl, rcp_matrix = None):
		if hkl is None:
			return None

		if rcp_matrix is None:
			rcp_matrix = self.rcp_matrix

		return Vector(np.dot(hkl, rcp_matrix))

	def isperpendicular(self, hkl1, hkl2):
		'''
			if these two hkls are perpendicular, return True
			else return False 
		'''

		vec1, vec2 = self.vec_in_rcp(hkl1), self.vec_in_rcp(hkl2)
		return vec1.isperpendicular(vec2)

	def isparallel(self, hkl1, hkl2):
		'''
			if these two hkls are parallel and have same direction, return 1;
			if reverse-parallel, return 2;
			else return 0
		'''
		vec1, vec2 = self.vec_in_rcp(hkl1), self.vec_in_rcp(hkl2)
		return vec1.isparallel(vec2)


class index(Vector):
	'''
		index of lattice plane in a Xtal
	'''
	def __new__(cls, input_array):

		if len(input_array) != 3:
			raise ValueError('Must be 3-dimension')
		input_array = [int(i) for i in input_array]

		return np.array(input_array).view(cls)

	@property
	def h(self):
		return self[0]

	@property
	def k(self):
		return self[1]

	@property
	def l(self):
		return self[2]

	@property
	def str(self):
		return str(self[0]) + ' ' + str(self[1]) + ' ' + str(self[2])


class Familyindex(index):
	'''
		index of a family of lattice plane in a Xtal
	'''

	def __new__(cls, input_array, symmetry = None):
		obj = super(Familyindex, cls).__new__(cls, input_array)
		obj.symmetry = symmetry
		return obj

	def __array_finalize__(self, obj):
		if obj is None: return
		self.symmetry = getattr(obj, 'symmetry', None)

	def add_sons(self, hkl):
		if not hasattr(self, '_sons'):
			self._sons = [hkl]
		else:
			self._sons.append(hkl)

	def sons(self):
		'''
			give the sons of this family of lattice plane
		'''
		if hasattr(self, '_sons'):
			return self._sons

		if self.symmetry == 'cubic':
			store = set()
			perms = list(itertools.permutations(self, 3))
			signs = list(itertools.product((1,-1),(1,-1),(1,-1)))
			for perm in perms:
				for sign in signs:
					i = [int(p*s) for p,s in zip(perm, sign)]
					store.add(tuple(i))
			self._sons = [index(s) for s in store]
			return self._sons

		if self.symmetry == 'hcp':
			store = set()
			perms_elem = (self.h, self.k, -(self.h + self.k))
			perms = list(itertools.permutations(perms_elem, 3))
			for perm in perms:
				i = (perm[0], perm[1], self.l)
				store.add(i)
				i = (perm[0], perm[1], -self.l)
				store.add(i)
			self._sons = [index(s) for s in store]
			return self._sons

		return ()

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


def Gen_hklfamilies(hklrange = (10,10,10), lattice = None):
	'''
		give an array of hklfamilies from given range of index
		NOTE:
		This is a GENERATOR
	'''
	if lattice is None:
		hkls = Gen_hkls(hklrange = hklrange)
		for hkl in hkls:
			hklfamily = Familyindex(hkl)
			hklfamily.add_sons(hkl)
			yield hklfamily
	else:

		if lattice.LP.iscubic():
			for h in range(hklrange[0] + 1):
				for k in range(min(h, hklrange[1]) + 1):
					for l in range(min(k, hklrange[2]) + 1):
						if h == 0 and k == 0 and  l == 0:
							continue
						hklfamily = Familyindex((h,k,l), symmetry = 'cubic')
						if not lattice.isextinct(hklfamily):
							yield hklfamily

		elif lattice.LP.ishcp():
			for h in range(hklrange[0] + 1):
				for k in range(min(h, hklrange[1]) + 1):
					for l in range(hklrange[2] + 1):
						if h == 0 and k == 0 and l == 0:
							continue
						hklfamily = Familyindex((h,k,l), symmetry = 'hcp')
						if not lattice.isextinct(hklfamily):
							yield hklfamily

		else:
			hkls = FT.tolist(Gen_hkls(hklrange = hklrange, lattice = None))
			d_spacing = lattice.D_spacing(hkls)
			d_spacing_0 = 0
			EPS_d = 1e-5
			for hkl,d_spacing in sorted(zip(hkls, d_spacing), key = lambda x: x[1], reverse = True):
				if FT.equal(d_spacing, d_spacing_0, EPS_d):
					hklfamily.add_sons(hkl)
				else:
					yield hklfamily
					hklfamily = Familyindex(hkl)
					hklfamily.add_sons(hkl)
			yield hklfamily



def Gen_hkls(hklrange = (10,10,10), lattice = None):
	'''
		give an array of hkl from given hklfamilies
		NOTE:
		This is a GENERATOR
	'''
	def sorted_range(end):
		return sorted(range(-end, end+1), key = abs)

	for h in sorted_range(hklrange[0]):
		for k in sorted_range(hklrange[1]):
			for l in sorted_range(hklrange[2]):
				if h == 0 and k == 0 and l == 0:
					continue
				hkl = index((h,k,l))
				if (lattice is None) or (not lattice.isextinct(hkl)):
					yield hkl

'''
	TEST
'''

if __name__ == '__main__':
	# l = Lattice(material = 'Mg')
	# l.LP.show()
	# hkls = FT.tolist(Gen_hklfamilies((3,3,3), lattice = l))
	# vs = l.vec_in_lattice(hkls)
	# ds = l.d_spacing(hkls)
	# for hkl,v, d in zip(hkls,vs, ds):
	# 	print(hkl,hkl.multiplicity,v, d)
	hkls = FT.tolist(Gen_hkls((3,3,3)))
	# print(l.Vec_in_lattice((index((1,0,0)), index((0,1,0)), index((0,0,1)))))
	# ele = Element('Cu')
	# print(ele.sc_factor(np.deg2rad(30), 0.5))
	# print(index((0,0,-1)).str)