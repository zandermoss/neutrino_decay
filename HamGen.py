import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from numpy import linalg as LA

import MatrixOps as MO
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


	def __init__(self,param,Um,tau,nu_mass,phi_mass,track,decay,dcy_channels,splines=None,regen=False):
		self.param=param
		self.ugen=SU.SUGen(param)
		self.splines=splines
		self.doupdate=False
		self.Um = Um
		self.tau=tau
		self.regen = regen
	
		#Sanity check on particle masses and pstar matrix precalculation	
		self.pstar = np.zeros([self.param.numneu,self.param.numneu],np.float64)
		for i in range(0,self.param.numneu):
			for j in range(i+1,self.param.numneu):
				if (dcy_channels[i,j] == True):
					if nu_mass[i] + phi_mass > nu_mass[j]:
						print "BAD: masses don't match up!"
					else:
						self.pstar[i,j] =(1.0/(2.0*nu_mass[j]))*((nu_mass[j]**2 - (nu_mass[i] + phi_mass)**2)*(nu_mass[j]**2 - (nu_mass[i] - phi_mass)**2))
		self.m_nu = nu_mass
		self.m_phi = phi_mass	
		self.dcy_channels = dcy_channels

		#Randomized self.parameters to generate conjugation matrices
		#We will work in the flavor basis, so Um and Ug map from 
		#the mass basis to the flavor basis and the decay basis 
		#to the flavor basis, respectively
		#Ugen=PC.PhysicsConstants()
		self.track=track
		erange = self.track.erange
		self.Hstart=np.zeros([len(erange),self.param.numneu,self.param.numneu],np.complex128)	
		self.Gamma=np.zeros([len(erange),self.param.numneu,self.param.numneu],np.complex128)	
		self.Int=np.zeros([self.param.numneu,self.param.numneu],np.complex128)
		
		#Check if matter effects or decay is present

		if self.splines==None:
			self.matter=False
			vmatter=False
			self.doupdate=False
		elif type(self.splines)==float:
			self.matter=True
			vmatter=False
			self.doupdate=False
		else:
			self.doupdate=True
			self.matter=True
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
		Md=np.zeros([self.param.numneu,self.param.numneu],np.complex128)

		#Add interaction term
		if (self.matter)&(vmatter==False):
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
			self.Hstart[ei,:,:]=Md


		if decay:	
			for ei in range(0,len(erange)):
				for i in range(0,self.param.numneu):
					for j in range(i+1,self.param.numneu):
						self.Gamma[ei,i,j] = (erange[ei]**(self.param.decay_power))*(nu_mass[j]/tau[i,j])

#		if matter:			
#			for ei in range(0,len(erange)):
#				self.H0[ei] += np.dot(Um.conj().T,np.dot(self.Int,Um))	

				
			

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


