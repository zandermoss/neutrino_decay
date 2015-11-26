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


class SUGen(object):

	def __init__(self,param):
		self.param=param

		self.lamb=np.zeros([self.param.numneu,self.param.numneu],complex)
		self.lamb[0,1]=0.5
		self.lamb[1,0]=0.25

	def sample_params(self):
		#do stuff
		H=np.zeros([self.param.numneu,self.param.numneu],complex)
	


	def matrix_gen(self):
		#do stuff
		d=self.param.numneu
		U=np.identity(d,complex)
		for m in range(0,d-1):
			for n in range(m+1,d):
				print "M,N:",m,",",n
				R=self.rotation_gen(m,n)
				P=self.phase_gen(m,n)
				U=np.dot(U,P)
				U=np.dot(U,R)
			

		return U
	
	def rotation_gen(self,m,n):
		"""
		Here, we generate the rotation-type matrices used
		to construct the matrix element of SU(N). This is 
		simply modifying the identity with sin and cos in 
		the appropriate positions.
		"""

		R=np.identity(self.param.numneu,complex)
		lcos=cmath.cos(self.lamb[m,n])
		lsin=cmath.sin(self.lamb[m,n])

		R[m,m]=lcos
		R[n,n]=lcos
		R[m,n]=(-1.0)*lsin
		R[n,m]=lsin

		return R

	def phase_gen(self,m,n):
		"""
		Now, we generate the phasing-type matrices used in 
		the construction. We again modify the identity, 
		but this time with complex phases.
		"""
	
		P=np.identity(self.param.numneu,complex)
		lpphase=cmath.exp(1.0j*self.lamb[n,m])
		lmphase=cmath.exp(-1.0j*self.lamb[n,m])

		P[m,m]=lpphase
		P[n,n]=lmphase

		return P	
