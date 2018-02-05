# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 10:55:55
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import myfunctools as FT

class Xray(object):
	'''
		Class for x-ray
		Attributes:
			wavelength: wavelength or wavelengthes of X-ray
			energy: energy or energies of X-ray
			intensity: corresponding intensities of X-ray
			spectrum: spectrum of X-ray, zip of wavelength and intensity
	'''
	def __init__(self,*, wavelength = None, energy = None, intensity = None, filename = None, islambda = True, EPS = None):
		super().__init__()

		if wavelength is None and energy is None:
			self.Gen_from_file(filename, islambda)
			if not (EPS is None):
				filt = (self.intensity > EPS)
				self.intensity = self.intensity[filt]
				self.wavelength = self.wavelength[filt]
				self.energy = self.energy[filt]
			self.Gen_tag()
			return

		if wavelength is None:
			self.energy = np.array(energy)
			if self.energy.shape == ():
				self.energy = self.energy[None]
			self.wavelength = self.lambda_from_energy(self.energy)

		if energy is None:
			self.wavelength = np.array(wavelength)
			if self.wavelength.shape == ():
				self.wavelength = self.wavelength[None]
			self.energy = self.energy_from_lambda(self.wavelength)

		if intensity is None:
			self.intensity = np.ones(self.wavelength.shape)
		else:
			self.intensity = np.array(intensity)
			if self.intensity.shape != self.wavelength.shape:
				logging.error('Shape of intensity and wavelength do not match!')
				raise ValueError
		self.Gen_tag()

	def lambda_from_energy(self, energy):
		return 12.4144/energy

	def energy_from_lambda(self, wavelength):
		return 12.4144/wavelength

	def Gen_from_file(self, filename, islambda = True):
		data = np.genfromtxt(filename)
		if islambda:
			self.wavelength = data[:,0]
			self.energy = self.energy_from_lambda(self.wavelength)
		else:
			self.energy = data[:,0]
			self.wavelength = self.lambda_from_energy(self.energy)
		self.intensity = data[:,1]

	def Add_energy(self, energy, intn = 1, tag = True):
		self.energy = np.concatenate((self.energy, np.array(energy)[None]))
		self.wavelength = np.concatenate((self.wavelength, np.array(self.lambda_from_energy(energy))[None]))
		self.intensity = np.concatenate((self.intensity, np.array(intn)[None]))
		if tag:
			self.Add_tag(energy = energy)
		if hasattr(self, 'inter_spectrum'):
			del self.inter_spectrum

	def Add_wavelength(self, wl, intn = 1, tag = True):
		self.wavelength = np.concatenate((self.wavelength, np.array(wl)[None]))
		self.energy = np.concatenate((self.energy, np.array(self.energy_from_lambda(wl))[None]))
		self.intensity = np.concatenate((self.intensity, np.array(intn)[None]))
		if tag:
			self.Add_tag(wl = wl)

		if hasattr(self, 'inter_spectrum'):
			del self.inter_spectrum

	def Add_tag(self, wl = None, energy = None):
		def find_nearest_idx(arr, value):
			return np.abs(arr - value).argmin()

		if wl is None and energy is None:
			raise TypeError("Empty Parameter")
		if not wl is None:
			w = wl
		elif not energy is None:
			w = self.lambda_from_energy(energy)

		w_idx = find_nearest_idx(self.wavelength, w)

		if not hasattr(self, 'tag_wl'):
			self.tag_wl = [self.wavelength[w_idx]]
			self.tag_wl_intn = [self.intensity[w_idx]]
		else:
			self.tag_wl.append(self.wavelength[w_idx])
			self.tag_wl_intn.append(self.intensity[w_idx])

		logging.debug('wl %f and intn %f is added tag'%(self.wavelength[w_idx], self.intensity[w_idx]))
	
	def Gen_tag(self, range_wl = None, range_energy = None):
		if not range_energy is None and range_wl is None:
			range_wl = self.lambda_from_energy(range_energy)

		if range_wl is None:
			wl = self.wavelength
			intn = self.intensity
		else:
			range_wl = np.sort(range_wl)
			idx = np.bitwise_and(range_wl[0] < self.wavelength, self.wavelength < range_wl[1])
			wl = self.wavelength[idx]
			intn = self.intensity[idx]

		idx_max = intn.argmax()
		self.Add_tag(wl[idx_max])

	def intn_wl(self, wl, EPS = 0):

		idx = FT.where(self.wavelength, wl)
		if idx:
			return self.intensity[idx - 1]

		if hasattr(self, 'inter_spectrum'):
			f_intn = self.inter_spectrum
		else:
			f_intn = interp1d(self.wavelength, self.intensity, kind  = 'cubic')
			self.inter_spectrum = f_intn

		try:
			intn = f_intn(wl)
		except ValueError as e:
			return 0

		return 0 if intn < EPS else intn
		
	@property
	def spectrum(self, islambda = True):
		xx = self.wavelength if islambda else self.energy
		if xx.shape == ():
			return np.array((xx, self.intensity))[None,:]
		return zip(xx, self.intensity)

	@property
	def num_tag(self):
		if not hasattr(self, 'tag_wl'):
			return 0
		return len(self.tag_wl)

	def show(self, islambda = True):
		xx = self.wavelength if islambda else self.energy
		xlabel = 'wavelength (A)' if islambda else 'energy (keV)'
		y = self.intensity

		plt.figure()
		plt.plot(xx, y,'.')
		if hasattr(self, 'tag_wl'):
			xx_tag = self.tag_wl if islambda else self.energy_from_lambda(self.tag_wl)
			yy_tag = self.tag_wl_intn
			plt.plot(xx_tag, yy_tag, 'o', mfc = 'none')
			for x_tag, y_tag in zip(xx_tag, yy_tag):
				plt.annotate('%.3f'%(x_tag), xy = (x_tag, y_tag), ha = 'center', va = 'bottom')
		plt.title('Spectrum')
		plt.xlabel(xlabel)
		plt.ylabel('Intensity (a. u.)')
		plt.show()

if __name__ == '__main__':
	# x = Xray(wavelength = np.arange(0.4, 0.5, 0.001))
	x = Xray(filename = 'u18_gap12mm.txt', islambda = False, EPS = 1e12)
	# x.Gen_tag(range_wl = (0.3,0.1))
	# x = Xray(wavelength = 0.5)
	# x.Gen_tag()
	print(x.intn_wl(0.53, EPS = 1e10))
	x.show()
	