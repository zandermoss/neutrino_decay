import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import random
import cmath
import math
import SUGen

## A daughter of SUGen specialized to the generation of PMNS matrices.

class PMNSGen(SUGen.SUGen):

	## The constructior
	# Initiaizes as a daughter of SUGen, and loads in the 3x3 mixing angles
	# (best global fit) from PhysConst
	def __init__(self,param):

		SUGen.SUGen.__init__(self,param)


		self.thetas=[param.th12,param.th13,param.th23]


	##Generates a matrix corresponding to the lambda parameters.
	# slight difference from matrix_gen() in SUGen in that no phases 
	# are rolled in. Ultimately, they will be, but they are currently 
	# neglected for the sake of simplicity.

	def matrix_gen(self):
		d=self.param.numneu
		pmns_pairs=[[0,1],[0,2],[1,2]]
		U=np.identity(d,complex)
		m=0
		while (m<=d-2):
			n=m+1
			while (n<=d-1):
				R=self.rotation_gen(m,n)
				#P=self.phase_gen(m,n) no cp
				#U=np.dot(U,P)
				U=np.dot(U,R)
				n+=1
			m+=1	
	
		if (self.param.neutype=='antineutrino'):
			return U.T

		else:
			
			return U.conj().T


	##Samples lambda parameters.
	# this is all the same as the SUGen sampling function, except that the 3x3
	# mixing angles are drawn from the PhysConst class (best global fit angles)
	# instead of being sampled. Ultimately, if there is a sterile neutrino, and 
	# it does decay, the global fits will change to some extent, so these angles 
	# will have to be fit cotemporally with the decay structure and extra mixing/
	# phasing angles. For the time being, they will be fixed for simplicity of 
	# analysis.

	def sample_params(self):
		d=self.param.numneu
		pmns_pairs=[[0,1],[0,2],[1,2]]

		m=1
		while (m<=d-1):
			n=0
			while (n<=m-1):
				if ([n,m] in pmns_pairs):
					self.lamb[n,m]=-1.0*self.thetas[pmns_pairs.index([n,m])]
				else:
					self.lamb[n,m]=-1.0*self.sample_angle(n,m)
				self.lamb[m,n]=self.sample_phase()

				n+=1
			m+=1

		for i in range(0,4):
			self.lamb[i,3]=0.0;



