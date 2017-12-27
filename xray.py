# -*- coding: utf-8 -*-
# @Date    : 2017-12-27 10:55:55
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import logging
import matplotlib.pyplot as plt

class Xray(object):
	'''
		Class for x-ray
		Attributes:
			wavelength: wavelength or wavelengthes of X-ray
			energy: energy or energies of X-ray
			intensity: corresponding intensities of X-ray
			spectrum: spectrum of X-ray, zip of wavelength and intensity
	'''
	def __init__(self,*, wavelength = None, energy = None, intensity = None, filename = None, islambda = True):
		super().__init__()

		if wavelength is None and energy is None:
			self.Gen_from_file(filename, islambda)
			return

		if wavelength is None:
			self.energy = np.array(energy)
			self.wavelength = self.lambda_from_energy(self.energy)

		if energy is None:
			self.wavelength = np.array(wavelength)
			self.energy = self.energy_from_lambda(self.wavelength)

		if intensity is None:
			self.intensity = np.ones(self.wavelength.shape)
		else:
			self.intensity = np.array(intensity)
			if self.intensity.shape != self.wavelength.shape:
				logging.error('Shape of intensity and wavelength do not match!')
				raise ValueError

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

	@property
	def spectrum(self, islambda = True):
		xx = self.wavelength if islambda else self.energy
		if xx.shape == ():
			return np.array((xx, self.intensity))[None,:]
		return zip(xx, self.intensity)

	def show(self, islambda = True):
		xx = self.wavelength if islambda else self.energy
		xlabel = 'wavelength (A)' if islambda else 'energy (keV)'
		y = self.intensity
		plt.plot(xx, y)
		plt.title('Spectrum')
		plt.xlabel(xlabel)
		plt.ylabel('Intensity (a. u.)')
		plt.show()

if __name__ == '__main__':
	x = Xray(filename = 'u18_gap12mm.txt', islambda = False)
	x.show(islambda = False)
	