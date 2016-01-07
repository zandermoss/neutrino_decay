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


class DeSolve(object):

	def __init__(self,myhamgen,param):
		self.hamgen = myhamgen
		self.param=param

		self.falf=0.0
		self.fstep=0.0


		#Set up the solver

		self.norm=1.0
		#Generate basis vectors
		self.b=[]
		for x in range(0,param.numneu):
			self.b.append(np.zeros(param.numneu))
			self.b[x][x]=1.0
		#---------------------

	#Time independent hamiltonian,
	#define as matrix product for DE solver:
	def func(self,t,y):
		H=self.hamgen.update(t/self.norm)
		#H=self.hamgen.update(0.5)
		return -1j*np.dot(H,y)

		
	

	def prop(self,x,i,j):	
		#Initial value: pure neutrino-0 state
		self.r=ode(self.func).set_integrator('zvode', method='adams',rtol=1e-10)

		self.norm=x[-1]

		y0=self.b[i]
		x0=0
	
		xf=x[-1]
		step=x[1]-x[0]

	
		self.r.set_initial_value(y0, x0)
		
		#Solve!
		dist=[]
		output=[]
		
		output.append(y0)

		self.fstep=step
		while self.r.successful() and self.r.t < xf:
			output.append(self.r.integrate(self.r.t+step))
			dist.append(self.r.t)
			self.falf=self.r.t
	
		amp=np.zeros(len(output))


		for k in range(0,len(output)):
			ip=np.dot(self.b[j],output[k])	
			amp[k]=np.absolute(ip)**2
	
		return amp
