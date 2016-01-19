import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d as gausfilt
import dill as pickle
import math
import random


class Track():
	def __init__(self,param,resolution,energy,randz=False):
		self.param=param
		self.l=2*self.param.EARTHRADIUS
		self.theta=0	
		self.resolution=resolution #steps/km
		self.step=1.0/resolution
		self.E=energy

		self.intersections=[None,None]
		self.shellradii=[0.1917,0.5462]
		self.shellangles=np.zeros(2)
		self.delta=0.009456


		for rad in enumerate(self.shellradii):
			self.shellangles[rad[0]]=2*math.asin(rad[1]/self.param.EARTHRADIUS)
		if randz:
			self.randomize_zenith()

	def calc_l(self,theta):
		R=self.param.EARTHRADIUS
		self.l = 2*R*math.cos(theta/2.0) 
	
	def r(self,x):
		r=math.sqrt(1.0+4.0*math.cos(self.theta/2.0)**2*(x**2-x))
		#return as fraction of total earth radius
		return r	

	def randomize_zenith(self):
		'''Generate a theta distribution corresponding
		to a uniform distribution over the surface of a sphere.
		Because the density profile and chord distance are both 
		phi-invariant, we do not bother with that variable here'''

		rand = random.random()		
		self.theta=math.acos(2*rand-1)
		
		self.calc_l(self.theta)

	def shell_intersection():
		self.calc_l(self.theta)

		for angle in enumerate(self.shellangles):
			if self.theta<=angle[1]:
				self.intersections[angle[0]]=True
			else:
				self.intersections[angle[0]]=False

		R=self.param.EARTHRADIUS
		s=np.zeros(2)
		for j in range(0,2):
			s[j]=R*math.cos(self.theta/2.0)-math.sqrt(r[j]**2-R**2*math.sin(self.theta/2.0)**2)
		
		return s
				






