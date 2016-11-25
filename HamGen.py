import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from numpy import linalg as LA

import PhysConst as PC
import SUGen as SU
import random
import cmath
import math


## Generates and updates the propagation hamiltonian. 
# HamGen generates the propagation hamiltonian according to the model parameters provided through param (PhysConst), Um (PMNSGen), Ug (SUGen). Also updates this hamiltonian throughout neutrino propagation (according to varying earth density) using geometry information from a track (Track) and density splines.


class HamGen(object):

	## The constructor
	# The utility of the constructor, other than providing local variables to the 
	# multiple input objects, is to provide a method for switching between different
	# propagation modes. There are two physical decisions to make before propagation.
	# The first concerns neutrino decay. If myeig_dcy is not provided as an argument,
	# The hamiltonian will be constructed with no decay term, leading to standard 
	# neutrino oscillation physics. If decay is desired, a SUGen object, myeig_dcy
	# is provided, with the details of the decay model.
	# The second decision concerns the matter potential. If no splines object
	# is provided, the matter potential term will be left off, 
	# and the propagation will proceed in vacuum. If splines is provided as a float,
	# then propagation will proceed in a constant potential corresponding to the 
	# value of that float (though it is easier just to alter the func() function 
	# in DeSolve to the desired constant. If a Splines object (class Splines) is 
	# provided, then the potential term will vary with the earth density, though 
	# the hamiltonian must be explicitly updated at each propegation step using the 
	# update() function.

	def __init__(self,param,Um,tau,nu_mass,phi_mass,track,decay,splines=None):
		self.param=param
		self.ugen=SU.SUGen(param)
		self.splines=splines
		self.doupdate=False



		#Randomized self.parameters to generate conjugation matrices
		#We will work in the flavor basis, so Um and Ug map from 
		#the mass basis to the flavor basis and the decay basis 
		#to the flavor basis, respectively
		#Ugen=PC.PhysicsConstants()
		self.track=track
		erange = self.track.erange
		self.H0=np.zeros([len(erange),self.param.numneu,self.param.numneu],complex)	
		self.Gamma=np.zeros([len(erange),self.param.numneu,self.param.numneu],complex)	
		self.Int=np.zeros([self.param.numneu,self.param.numneu],complex)
		
		#Check if matter effects or decay is present

		if self.splines==None:
			matter=False
			vmatter=False
			self.doupdate=False
		elif type(self.splines)==float:
			matter=True
			vmatter=False
			self.doupdate=False
		else:
			self.doupdate=True
			matter=True
			vmatter=True
			self.yespline=splines.GetYe()
			#self.dspline=splines.GetEarthLine()
			self.dspline=splines.GetEarth()


		#Generate conj matrices
		#use known mixing self.parameters for M
		#Um=MT.calcU(self.param)
	
		#if decay:	
			#self.ugen.sample_params()
			#Ug=self.ugen.matrix_gen()	
		
		#Ugen.randomize_trig()
		#Ug=MT.calcU(Ugen)
		
		#Fill in mass and decay eigenvalues
		Md=np.zeros([self.param.numneu,self.param.numneu],complex)

		#Add interaction term
		if (matter)&(vmatter==False):
			potential=self.splines 
			#assume electron density is neutron density. This is roughly true.
			for flv in range(0,3):
				if flv==0:
					self.Int[flv,flv]=potential/2
				else:
					self.Int[flv,flv]=-potential/2


		for ei in range(0,len(erange)):
			for i in range(0,self.param.numneu):
				Md[i,i]= self.param.dm2[1,i+1]/(2*erange[ei])
				self.H0[ei,:,:]=M

		if decay:	
			for ei in range(0,len(erange)):
				for i in range(0,self.param.numneu):
					for j in range(i+1,self.param.numneu):
						self.Gamma[ei,i,j] = (erange[ei]**(self.param.decay_power))(nu_mass[j]/tau[i,j])

		if matter:			
			for ei in range(0,len(erange)):
				self.H0[ei] += np.dot(Um.conj().T,np.dot(self.Int,Um))	


	## This function updates the propagation hamiltonian as the earth radius changes.
	# Uses the track object to compute earth radius from fractional progress along 
	# the neutrino track, and then uses the earth density spline and electron 
	# fraction splines from the splines object to calculate the matter potential.
	# Returns a hamiltonian matrix with properly updated matter potential term.
	# checks the param object to determine whether the particle propagating is a 
	# neutrino or an antineutrino. If antineutrino, the sign of the potential term
	# is reversed.
	# @par x is the fractional distance traversed along the track (distance/track length).
	# @return an updated hamiltonian matrix (2d np array)


	def H0_Update(self,x):
		if self.doupdate==True:
			r=self.track.r(x)
			
			ye=self.yespline(r)
			density=self.dspline(r) #g/cm^3
			nd=density*self.param.gr/(self.param.GeV*self.param.proton_mass) #convert to #proton/cm^3
			npotential=math.sqrt(2)*self.param.GF*(self.param.cm)**(-3)*nd*(1-ye) #convert cm to 1/eV

			ppotential=math.sqrt(2)*self.param.GF*(self.param.cm)**(-3)*nd*ye #convert cm to 1/eV
			#print x*self.track.l*self.param.km, ",",npotential
	
			for flv in range(0,3):
				#assume electron density is neutron density. This is roughly true.
				if flv==0:
					self.Int[flv,flv]=ppotential-0.5*npotential
				else:
					self.Int[flv,flv]=-0.5*npotential

			ret_array=np.zeros([len(self.track.erange),self.param.numneu,self.param.numneu],complex)	
		
			MassInt = np.dot(Um.conj().T,np.dot(self.Int,Um))

			if (self.param.neutype=='antineutrino'):
				for ei in range(0,len(self.track.erange)):
					ret_array[ei,:,:] = self.H0+MassInt*-1.0

			else:
				for ei in range(0,len(self.track.erange)):
					ret_array[ei,:,:] = self.H0+MassInt

			return ret_array
		else:
			return self.H0
	
	def GetGamma(self):
		return self.Gamma
