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
				print "M,N:",m,n
				if ([m,n] in pmns_pairs):
					R=self.pmnsr_gen(m,n,self.thetas[pmns_pairs.index([m,n])])
				else:
					R=self.rotation_gen(m,n)
				#P=self.phase_gen(m,n) no cp
				#U=np.dot(U,P)
				U=np.dot(U,R)
				n+=1
			m+=1	

		return U.conj().T

	def pmnsr_gen(self,m,n,angle):
		"""
		Here, we generate the rotation-type matrices used
		to construct the matrix element of SU(N). This is 
		simply modifying the identity with sin and cos in 
		the appropriate positions.
		"""
		print m,n,"ANGLE:",angle	

		R=np.identity(self.param.numneu,complex)
		lcos=cmath.cos(angle)
		lsin=cmath.sin(-1.0*angle)

		R[m,m]=lcos
		R[n,n]=lcos
		R[n,m]=(-1.0)*lsin
		R[m,n]=lsin
		
		return R

