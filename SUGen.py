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

	def sample_params(self):
		d=self.param.numneu

		m=1
		while (m<=d-1):
			n=0
			while (n<=m-1):
				self.lamb[m,n]=self.sample_phase()
				self.lamb[n,m]=self.sample_angle(n,m)
				n+=1
			m+=1

	def sample_phase(self):	
		u=random.random()
		phi=u*math.pi*2
		return phi

	def sample_angle(self,m,n):	
		u=random.random()
		k=2*(n-m)-1
		theta=math.acos((1-u)**(1.0/float(k+1)))	
		return theta

	def matrix_gen(self):
		d=self.param.numneu
		U=np.identity(d,complex)
		m=0
		while (m<=d-2):
			n=m+1
			while (n<=d-1):
				R=self.rotation_gen(m,n)
				P=self.phase_gen(m,n)
				U=np.dot(U,P)
				U=np.dot(U,R)
				n+=1
			m+=1	

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
		R[n,m]=(-1.0)*lsin
		R[m,n]=lsin

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
