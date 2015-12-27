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


class NumSolve(object):

	def __init__(self,H,param):
		self.H = H
		self.param=param
		w,N = np.linalg.eig(H)
		self.w=w
		self.N=N
		self.Ht=np.zeros([param.numneu,param.numneu],complex)
		self.Ni=np.linalg.inv(self.N)


	
	def prop(self,l):

		for i in range(0,self.param.numneu):
			self.Ht[i,i] = cmath.exp(-1j*self.w[i]*l)
		
		T= np.dot(self.N,np.dot(self.Ht,self.Ni))
		return T


	def scalar_prop(self,l,i,j):
	
		p=self.prop(l)
		P=p*(p.conj())
		amp=P[j,i].real
		return amp	
