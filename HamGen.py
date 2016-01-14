import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import SUGen as SU
import random
import cmath
import math


class HamGen(object):

	def __init__(self,param,Ug,track,myeig_dcy=None,splines=None):
		self.param=param
		self.ugen=SU.SUGen(param)
		self.eig_dcy=myeig_dcy
		self.splines=splines

		#Randomized self.parameters to generate conjugation matrices
		#We will work in the flavor basis, so Um and Ug map from 
		#the mass basis to the flavor basis and the decay basis 
		#to the flavor basis, respectively
		#Ugen=PC.PhysicsConstants()
		self.track=track
		self.H=np.zeros([self.param.numneu,self.param.numneu],complex)	
		self.Int=np.zeros([self.param.numneu,self.param.numneu],complex)
		
		#Check if matter effects or decay is present
		if self.eig_dcy==None:
			decay=False
		else:
			decay=True

		if type(self.splines)==float:
			vmatter=False
		else:
			vmatter=True
			self.yespline=splines.GetYe()
			self.dspline=splines.GetEarth()


		if self.splines==None:
			matter=False
		else:
			matter=True


		#Generate conj matrices
		#use known mixing self.parameters for M
		Um=MT.calcU(self.param)
	
		#if decay:	
			#self.ugen.sample_params()
			#Ug=self.ugen.matrix_gen()	
		
		#Ugen.randomize_trig()
		#Ug=MT.calcU(Ugen)
		
		#Fill in mass and decay eigenvalues
		Md=np.zeros([self.param.numneu,self.param.numneu],complex)
		if decay:	
			Gd=np.zeros([self.param.numneu,self.param.numneu],complex)

		#Add interaction term
		if (matter)&(vmatter==False):
			potential=self.splines 
			#assume electron density is neutron density. This is roughly true.
			for flv in range(0,3):
				if flv==0:
					self.Int[flv,flv]=potential/2
				else:
					self.Int[flv,flv]=-potential/2

		for i in range(0,self.param.numneu):
		
			Md[i,i]= self.param.dm2[1,i+1]/(2*self.track.E)
		
			#print "Md:", Md	
			if decay:	
				Gd[i,i]= self.eig_dcy[i]*(self.track.E**(self.param.decay_power)) #Power law energy dependence 
		
	
		M= np.dot(Um,np.dot(Md,Um.conj().T))
		if decay:	
			G= np.dot(Ug,np.dot(Gd,Ug.conj().T))

	
		#Assemble Hamiltonian
		self.H=M

		if decay:
			self.H += -1j*G 

		if matter:			
			self.H += self.Int	

	
	def update(self,x):
		r=self.track.r(x)
		ye=self.yespline(r)
		density=self.dspline(r) #g/cm^3
		nd=density*self.param.gr/(self.param.GeV*self.param.proton_mass) #convert to #proton/cm^3
		npotential=math.sqrt(2)*self.param.GF*(self.param.cm)**(-3)*nd*(1-ye) #convert cm to 1/eV
		ppotential=math.sqrt(2)*self.param.GF*(self.param.cm)**(-3)*nd*ye #convert cm to 1/eV


		for flv in range(0,3):
			#assume electron density is neutron density. This is roughly true.
			if flv==0:
				self.Int[flv,flv]=ppotential-0.5*npotential
			else:
				self.Int[flv,flv]=-0.5*npotential

		return self.H+self.Int
