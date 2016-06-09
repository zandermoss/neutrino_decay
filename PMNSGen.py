import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import random
import cmath
import math
import SUGen


class PMNSGen(SUGen.SUGen):

	def __init__(self,param):

		SUGen.SUGen.__init__(self,param)


		self.thetas=[param.th12,param.th13,param.th23]

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



