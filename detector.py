# -*- coding: utf-8 -*-
# @Date    : 2018-01-30 16:24:01
# @Author  : J. W. Huang (huangjasper@126.com)

import numpy as np
import sys
import matplotlib.pyplot as plt
import logging; logging.basicConfig(level = logging.INFO)

import singleXtal as SX
import polyXtal as PX
import lattice as LTTC
import xray as XR
import peak as PK
from Vec import Vector
import LUT
import simu2d as S2D
import simu1d as S1D

class Detector(object):
	def __init__(self, normal = None, size = None, poni = (0,0), ps = None, x = None, y = None, dist = None, *, det = None):
		super(Detector, self).__init__()
		if not det is None:
			self.__init__(*LUT.detector[det])

		self.normal = normal
		self.dist = dist
		self.size = size
		self.poni = self.str2poni(poni)
		self.ps = ps
		self.x = x
		self.y = y
		self.Calc_geometry()

	def set(self, normal = None, dist = None, size = None, poni = None, ps = None, x = None, y = None, *, det = None):
		if not det is None:
			self.set(*LUT.detector[det])

		if not normal is None:
			self.normal = normal
		if not dist is None:
			self.dist = dist
		if not size is None:
			self.size = size
		if not x is None:
			self.x = x
		if not y is None:
			self.y = y
		if not poni is None:
			self.poni = self.str2poni(poni)
		if not ps is None:
			self.ps = ps
		self.Calc_geometry()

	def str2poni(self, ss):
		if ss is None:
			return None

		position = {
				'c' : 0.5,
				'l' : 0,
				'r' : 1,
				'u' : 1,
				'd' : 0,
			}

		p = [0,0]
		for i,s in enumerate(ss):
			if not isinstance(s, str):
				p[i] = s
				continue

			if not s in position:
				raise ValueError('Unknown Position!')

			if not hasattr(self, 'size'):
				raise ValueError('Unknown Size of Detector!')

			p[i] = self.size[i] * position[s]

		logging.debug('poni is %s'%(str(p)))
		return p

	def Calc_geometry(self):
		num_of_None = np.sum([1 if vec is None else 0 for vec in (self.normal, self.x, self.y)])
		logging.debug('num of None: %d'%(num_of_None))

		if num_of_None >= 2:
			if not self.normal is None:
				logging.info('Only normal is known, x & y will be guessed')
				normal = self.normal.norm
				if normal.isparallel(Vector(0,0,1)):
					x = Vector((-1,0,0))
					y = Vector((0,1,0))
				elif normal.isparallel(Vector(0,0,-1)):
					x = Vector((1,0,0))
					y = Vector((0,1,0))
				else:
					x = Vector(normal[1],-normal[0],0)
					y = normal.cross(x)
				self.normal = normal
				self.x = x
				self.y = y
				logging.info('x and y are guessed as %s, %s, respectively'%(str(x), str(y)))
			else:
				raise ValueError('Geometry is not decided!')

		if num_of_None == 0:
			normal = Vector(self.normal).norm
			x = Vector(self.x).norm
			y = Vector(self.y).norm
			if not normal.isperpendicular(x) or not normal.isperpendicular(y) or not x.isperpendicular(y):
				raise ValueError('x,y,normal of detector is not perpendicular!')
			self.normal = normal
			self.x = x
			self.y = y
			return

		if num_of_None == 1:
			if self.normal is None:
				logging.info('Normal is empty')
				x = Vector(self.x).norm
				y = Vector(self.y).norm
				if not x.isperpendicular(y):
					logging.debug('x is %r, y is %r'%(x, y))
					raise ValueError('x,y of detector is not perpendicular!')
				self.x = x
				self.y = y
				self.normal = x.cross(y)
				logging.debug('normal is %r'%(self.normal))
				return

			if self.x is None:
				logging.info('x is empty')
				normal = Vector(self.normal).norm
				y = Vector(self.y).norm
				if not normal.isperpendicular(y):
					logging.debug('normal is %r, y is %r'%(normal, y))
					raise ValueError('normal,y of detector is not perpendicular!')
				self.normal = normal
				self.y = y
				self.x = y.cross(normal)
				logging.debug('x is %r'%(self.x))
				return

			if self.y is None:
				logging.info('y is empty')
				normal = Vector(self.normal).norm
				x = Vector(self.x).norm
				if not normal.isperpendicular(x):
					logging.debug('normal is %r, x is %r'%(normal, x))
					raise ValueError('normal,x of detector is not perpendicular!')
				self.normal = normal
				self.x = x
				self.y = normal.cross(x)
				logging.debug('y is %r'%(self.y))
				return

	def rotate_by(self, q = None, axis = None, degree = None):
		self.normal = self.normal.rotate_by(q = q, axis = axis, degree = degree)
		self.x = self.x.rotate_by(q = q, axis = axis, degree = degree)
		self.y = self.y.rotate_by(q = q, axis = axis, degree = degree)

	def on(self, v):

		vx, vy = v
		if 0 <= vx <= self.size[0] and 0 <= vy <= self.size[1]:
			return True
		return False

	def Calc_tthgam_map(self, inc = None, x = None):
		if inc is None:
			inc = self.normal

		if x is None:
			x = Vector(inc[1], -inc[0], 0)
			if inc[1] < 0:
				x = - x

		n = self.normal

		logging.debug('inc = %r, vx = %r, normal = %r'%(inc, x, n))

		def tthgam_map(tth, gamma, on = True):
			v = Vector(tth = tth, gamma = gamma, z = inc, x = x)
			p = n.dot(v)
			logging.debug('v: %r, p: %r'%(v, p))
			if p <= 0:
				return None
			v_in_det = self.dist * (v/p - n)
			vx = v_in_det.dot(self.x) + self.poni[0]
			vy = v_in_det.dot(self.y) + self.poni[1]
			logging.debug('tth: %3f, gamma: %3f, vx: %.3f, vy: %.3f'%(np.rad2deg(tth), np.rad2deg(gamma), vx,vy))
			if on and not self.on((vx,vy)):
				return None
			return (vx,vy)

		self.map = tthgam_map
		self.direct_beam = self.map(0, 0, on = False)
		if self.direct_beam is None:
			logging.warn('Direct beam is not on plane of detector!')
		else:
			logging.debug('direct beam: (%.3f,%.3f)'%(self.direct_beam[0], self.direct_beam[1]))

	def Calc_peaks(self, p):

		self.Calc_tthgam_map(inc = p.inc, x = p.vx)

		if isinstance(p, S1D.Profile1D):
			self.patterndim = '1D'
			self.Calc_peaks_1D(p)

		elif isinstance(p, S2D.Pattern2d):
			self.patterndim = '2D'
			self.Calc_peaks_2D(p)

		else:
			raise ValueError('Unknown Pattern or Profile!')


	def Calc_peaks_2D(self, pattern):

		peaks = []
		for peak in pattern.peaks:
			coor = self.map(peak.tth, peak.gamma)
			if coor is None:
				continue
			peak.x_on_det, peak.y_on_det = coor
			logging.debug('peak: %s,grain: %d,tth: %f,gamma: %f, x: %f, y:%f'%(str(peak.index), peak.grain.number, peak.tth, peak.gamma, peak.x_on_det, peak.y_on_det))
			peaks.append(peak)

		peaks = np.array(peaks)
		self.peaks = peaks[None] if peaks.ndim == 0 else peaks

	def Calc_peaks_1D(self, profile):
		peaks = []
		for peak in profile.peaks:
			tth = peak.tth
			line = []
			isformer = False
			for gamma in np.linspace(0, 2*np.pi, 1000):
				v = self.map(tth, gamma)
				if v is None:
					if isformer: # end of line
						line.append(np.array(line_tmp))
						del line_tmp
					else: # among non-line
						pass
					isformer = False
				else: 
					if isformer: # among line
						line_tmp.append(v)
					else: # start of line
						line_tmp = [v]
					isformer = True
			if 'line_tmp' in dir(): # end of line
				line.append(np.array(line_tmp))
			if len(line) > 0:
				peak.line = list(line)
				peaks.append(peak)

		peaks = np.array(peaks)
		self.peaks = peaks[None] if peaks.ndim == 0 else peaks

	def show(self, *, islabel = True):
		if not hasattr(self, 'patterndim'):
			raise ValueError('No pattern!')

		if self.patterndim == '1D':
			self.show_1D(islabel = islabel)
		elif self.patterndim == '2D':
			self.show_2D(islabel = islabel)
		else:
			raise ValueError('No pattern!')

	def show_2D(self, *, islabel = True):
		def onpick(event):
			ind = event.ind
			peak = self.peaks[ind[0]]
			fig_txt_pop = plt.figure('Properties',figsize = (2,2))
			ax = fig_txt_pop.add_axes([0,0,1,1])
			peak_attr = 'index: %d %d %d\n'%(peak.index[0], peak.index[1], peak.index[2])\
						+ 'tth , gam' + ': %.2f, %.2f\n'%(np.rad2deg(peak.tth), np.rad2deg(peak.gamma))\
						+ 'x, y: %.3f, %.3f\n'%(peak.x_on_det, peak.y_on_det)\
						+ 'grain: %d\n'%(peak.grain.number)\
						+ 'd_spacing: %.4f A\n'%(peak.d_spacing)\
						+ 'wavelength: %.5f A\n'%(peak.wl)\
						+ 'energy: %.5f keV\n'%(peak.energy)\
						+ 'intensity: %.2e\n'%(peak.intn)
			ax.text(0, 1, peak_attr, weight = 'bold', style = 'italic', ha = 'left', va = 'top', transform =  ax.transAxes)
			ax.set_axis_off()
			plt.show()

		intn = np.array([peak.intn for peak in self.peaks])
		if len(intn) == 0:
			logging.error('No peaks on detector!')
			return
		intn_min, intn_max = intn.min(), intn.max()
		grey = plt.get_cmap('Greys')
		color = grey(intn / intn_max)
		fig = plt.figure(figsize = (self.size[0]/50, self.size[1]/50))
		ax = fig.gca()
		ax.set(xlabel = r'x', ylabel = r'y', title = 'Pattern on detector', xlim = (0, self.size[0]), ylim = (0, self.size[1]))
		if not self.poni is None:
			ax.scatter((self.poni[0],),(self.poni[1]), s = 100, marker = 's', c = 'tab:orange', label = 'poni')
		if not self.direct_beam is None: 
			ax.scatter((self.direct_beam[0]), (self.direct_beam[1]), s = 100, marker = '*', c = 'r', label = 'direct beam')
		s = np.array([(peak.x_on_det, peak.y_on_det) for peak in self.peaks])
		ax.scatter(s[:,0],s[:,1], c = color, s = 10, picker = True)
		for peak in self.peaks:
			logging.debug('peak: index: %s, x: %f, y:%f'%(peak.index.str, peak.x_on_det, peak.y_on_det))
			# print('peak: index: %s,tth: %f, gamma: %f,x: %f, y:%f'%(peak.index.str, np.rad2deg(peak.tth), np.rad2deg(peak.gamma), peak.x_on_det, peak.y_on_det))
			# ax.scatter(peak.x_on_det, peak.y_on_det, c = color_f(peak.intn), s = 10)
			if islabel:
				ax.annotate(peak.index.str, xy = (peak.x_on_det, peak.y_on_det), ha = 'center', va = 'bottom', size = 6)
		fig.canvas.mpl_connect('pick_event', onpick)

		plt.legend()
		plt.show()

	def show_1D(self,*, islabel = True):
		intn = np.array([peak.intn for peak in self.peaks])
		if len(intn) == 0:
			logging.error('No peaks on detector!')
			return
		intn_min, intn_max = intn.min(), intn.max()
		size_f = lambda x: x/intn_max
		fig = plt.figure(figsize = (self.size[0]/50, self.size[1]/50))
		ax = fig.gca()
		ax.set(xlabel = r'x', ylabel = r'y', title = 'Pattern on detector', xlim = (0, self.size[0]), ylim = (0, self.size[1]))
		if not self.poni is None:
			ax.scatter((self.poni[0],),(self.poni[1]), s = 100, marker = 's', c = 'tab:orange', label = 'poni')
		if not self.direct_beam is None: 
			ax.scatter((self.direct_beam[0]), (self.direct_beam[1]), s = 100, marker = '*', c = 'r', label = 'direct beam')
		for i,peak in enumerate(self.peaks):
			logging.info('1d peak: index %s, tth %f, intn %f'%(peak.index.str, np.rad2deg(peak.tth), peak.intn))
			# print(peak.index.str)
			for line in peak.line:
				# print(line.shape)
				ax.plot(line[:,0], line[:,1], linewidth = 2*size_f(peak.intn), color = LUT.colors[i+1], label = peak.index.str)
			# ax.annotate(peak.index.str, xy = (peak.line[0][0,0], peak.line[0][0,1]), ha = 'left', va = 'bottom', size = 8)

		plt.legend()
		plt.show()

if __name__ == '__main__':
	Ta = LTTC.Lattice(material = 'Ta')
	Cu = LTTC.Lattice(material = 'Cu')
	Mg = LTTC.Lattice(material = 'Mg')
	Al = LTTC.Lattice(material = 'Al')
	# sx = SX.SingleXtal(Al, z = (1,1,1), x = (-1,1,0))
	# sx.rotate_by(axis = (1,0,0), degrees = -6)
	# sx.strain1d(axis = (1,0,0), ratio = -0.03)
	px = PX.PolyXtal(Al)
	# xr = XR.Xray(filename = 'u18_gap12mm.txt', islambda = False, EPS = 1e13)
	# xr.show()
	xr = XR.Xray(wavelength = 0.53)
	inc = Vector(0,0,-1)
	# inc = inc.rotate_by(axis = (1,0,0), degree = 10)
	# p = S2D.Pattern2d(sx = sx, xray = xr, inc = Vector(0,0,-1))
	# p.Calc(hklrange = (5,5,5))
	p = S1D.Profile1D(xray = xr, px = px)
	p.Calc(range_2th = (10, 50), precision = 0.001)
	# p.show()
	p.Add_geometry(inc = inc, x = Vector(1,0,0))
	# det1 = Detector(normal = -inc , size = (400, 400), poni = ('c','c'), x = (1,0,0), dist = 250)
	det1 = Detector(normal = inc , size = (400, 400), poni = ('c','c'), dist = 250)
	# det1.rotate_by(axis = (1,1,0), degree = 60)
	# det1.set(poni = (-50,'c'))
	det1.Calc_peaks(p)
	det1.show()
