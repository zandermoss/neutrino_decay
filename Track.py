import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d as gausfilt
import dill as pickle
import math
import random

## A class used to handle track geometry
# Provides functions to calculate track length from angle and to convert
# from the distance along a given track to earth radius (necessary for matter
# potential calculations in HamGen::Update() ). Also provides a function to 
# randomize the zenith angle, drawing from a distribution which is flat over the 
# surface of the earth (not flat in [0,pi]).

class Track():


## The constructor.
# Initializes parameters
# @par param a PhysConst object
# @par resolution the number of steps/km to be used by the numerical integrator.
# @par energy the neutrino energy
# @par zenith the neutrino zenith angle
# @par randz a boolean switch allowing the randomization of zenith angle.

	def __init__(self,param,resolution,energy,theta,randz=False):
		self.param=param
		self.l=2*self.param.EARTHRADIUS
		self.theta=theta	
		self.resolution=resolution #steps/km
		self.step=1.0/resolution
		self.E=energy

		self.intersections=[None,None]
		self.shellradii=[0.1917,0.5462]
		self.shellangles=np.zeros(2)
		self.delta=0.009456


		for rad in enumerate(self.shellradii):
			self.shellangles[rad[0]]=2*math.asin(rad[1])
		if randz:
			self.randomize_zenith()

##Calculates track length
# Uses theta parameter to determine the distance between the south pole and
# the point of neutrino production in km.

	def calc_l(self):
		R=self.param.EARTHRADIUS
		self.l = 2*R*math.cos(self.param.PI-self.theta) 
	

##Calculates radial position of neutrino based on position along track.
# @par x the fractional distance traveled along the neutrino track.
# @return the radial component of neutrino position as a fraction of earth radius.
	def r(self,x):
		r=math.sqrt(1.0+4.0*math.cos(self.param.PI-self.theta)**2*(x**2-x))
		#return as fraction of total earth radius
		return r	

##Samples a zenith angle from a distribution flat over the surface of the earth
# Samples the angle and sets the member variable self.theta equal to it.

	def randomize_zenith(self):
		'''Generate a theta distribution corresponding
		to a uniform distribution over the surface of a sphere.
		Because the density profile and chord distance are both 
		phi-invariant, we do not bother with that variable here'''

		rand = random.random()		
		self.theta=self.param.PI-math.acos(2*rand-1)/2.0
		
		self.calc_l()


##Experimental: was used in an attempt to speed propagation.
# The idea was to do wide steps away from density transitions and concentrate
# the integration resolution at the transition point. As it happens, the algorithm
# (dopri5) has an adaptive step-size which is designed to handle regions of great
# variablity.

	def shell_intersection(self):
		self.calc_l()

		for angle in enumerate(self.shellangles):
			if 2*(self.param.PI-self.theta)<=angle[1]:
				self.intersections[angle[0]]=True
			else:
				self.intersections[angle[0]]=False

		R=self.param.EARTHRADIUS
		s=np.zeros(2)
		for j in range(0,2):
			if self.intersections[j]==True: 
				s[j]=R*math.cos(self.param.PI-self.theta)-math.sqrt(R**2*self.shellradii[j]**2-R**2*math.sin(self.param.PI-self.theta)**2)
				s[j]/=self.l
		
		return s
				






