# -*- coding: utf-8 -*-
# @Date    : 2017-12-26 10:18:07
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging


import lattice as LTTC
import myfunctools as FT
from Vec import Vector

class PolyXtal(object):
	'''
		Class for poly crystal
		Attributes:
			lattice: lattice information of this poly crystal
			rcp_matrix: a copy of reciprocal_matrix attribute in Lattice

		Method:
			d_spacing: d_spacing for specific hkl
	'''
	def __init__(self, lttc):
		super(PolyXtal, self).__init__()

		self.lattice = lttc
		self.rcp_matrix = self.lattice.rcp_matrix

	def Gen_d_spacing(self, hklfamilies):

		return self.lattice.Gen_d_spacing(hklfamilies)

	def d_spacing(self, hklfamily):
		return self.lattice.d_spacing(hklfamily)

	def Gen_intn(self, hklfamilies):
		for strct_factor in self.lattice.Gen_strct_factor(hklfamilies):
			if strct_factor is None:
				yield 0
			else:
				yield np.abs(strct_factor)**2

	def Calc_rcp_space(self, hklrange = (10,10,10), sorting = False, *, hkls = None):
		if hkls:
			hklfamilies = []
			for hkl in hkls:
				if self.lattice.isextinct(hkl):
					logging.info("%s is eliminated because of extinction"%(hkl.str))
				else:
					hklfamilies.append(hkl)				
		else:
			hklfamilies = list(LTTC.Gen_hklfamilies(hklrange = hklrange, lattice = self.lattice))
		
		if sorting:
			hklfamilies = sorted(hklfamilies, key = self.d_spacing, reverse = True)

		self.rcp_space_num = len(hklfamilies)
		self.rcp_space_hkls = hklfamilies
		self.rcp_space_d = [1/d for d in self.Gen_d_spacing(hklfamilies)]
		self.rcp_space_intn = list(self.Gen_intn(hklfamilies))

	def Save_simple_rcp_space(self, filename):
		if (not hasattr(self, 'rcp_space_num')) or self.rcp_space_num == 0:
			raise ValueError('No rcp_space')

		with open(filename, 'w') as f:
			f.writelines('%d\n'%(self.rcp_space_num))
			f.writelines('h k l radius intn\n')
			for hklfamily, radius, intn in zip(self.rcp_space_hkls, self.rcp_space_d, self.rcp_space_intn):
				f.writelines('%s %f %f\n'%(hklfamily.str, radius, intn))

	def Save_proj_rcp_space(self, filename, res_gamma = 0.1):
		if (not hasattr(self, 'rcp_space_num')) or self.rcp_space_num == 0:
			raise ValueError('No rcp_space')

		gammas = np.linspace(0, 2*np.pi, 360// res_gamma, endpoint = False)
		num = self.rcp_space_num*len(gammas)
		with open(filename, 'w') as f:
			f.writelines('%d\n'%(num))
			f.writelines('h k l kx ky intn\n')
			for hklfamily, radius, intn in zip(self.rcp_space_hkls, self.rcp_space_d, self.rcp_space_intn):
				for gamma in gammas:
					kx = radius * np.cos(gamma)
					ky = radius * np.sin(gamma)
					f.writelines('%s %f %f %f\n'%(hklfamily.str, kx, ky, intn))

	def Save_rcp_space(self, filename, res_tth = 0.1, res_gamma = 0.1):
		if (not hasattr(self, 'rcp_space_num')) or self.rcp_space_num == 0:
			raise ValueError('No rcp_space')

		def mesh_tthgam(tths, gammas):
			for tth in tths:
				if tth in (0,np.pi):
					yield (tth, 0)
				else:
					for gamma in gammas:
						yield (tth, gamma)

		tths = np.linspace(0, np.pi, 180//res_tth, endpoint = True)
		gammas = np.linspace(0, 2*np.pi, 360//res_gamma, endpoint = False)
		num = self.rcp_space_num*len(list(mesh_tthgam(tths,gammas)))
		with open(filename, 'w') as f:
			f.writelines('%d\n'%(num))
			f.writelines('h k l kx ky kz intn\n')
			for hklfamily, radius, intn in zip(self.rcp_space_hkls, self.rcp_space_d, self.rcp_space_intn):
				for tth, gamma in mesh_tthgam(tths, gammas):
					vec = Vector(tth = tth, gamma = gamma)
					# print(tth, gamma, vec, radius)
					f.writelines('%s %f %f %f %f\n'%(hklfamily.str, radius*vec[0], radius*vec[1], radius*vec[2], intn))


if __name__ == '__main__':
	Cu = LTTC.Lattice(material = 'Cu')
	Mg = LTTC.Lattice(material = 'Mg')
	X = PolyXtal(Mg)
	hkls = (LTTC.Familyindex((1,0,0)), LTTC.Familyindex((0,0,2)), LTTC.Familyindex((1,0,1)))
	# print(list(X.Gen_d_spacing(hkl)))
	# hkls = LTTC.Familyindex(())
	# X.Calc_rcp_space((3,3,3))
	# X.Save_proj_rcp_space('rcp3.dat', res_gamma = 1)
	X.Calc_rcp_space(hkls = hkls)
	X.Save_proj_rcp_space('rcp.dat', res_gamma = 1)
