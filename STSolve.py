import numpy as np

import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import random
import cmath
import math

import DeSolve

class STSolve(DeSolve.DeSolve):

	def __init__(self,myhamgen,param):

		DeSolve.DeSolve.__init__(self,myhamgen,param)

		self.H0=self.hamgen.H0


	#Time dependent hamiltonian,
	#define as matrix product for DE solver:
	def func(self,t,y):
		self.hamgen.update(t/self.norm)
		Int=self.hamgen.Int
		expm=LA.expm(1.0j*self.H0*t)
		iexpm=LA.expm(-1.0j*self.H0*t)
		H_I1=np.dot(expm,np.dot(Int,iexpm))
		
		#H=self.hamgen.update(0.5)
		return -1j*np.dot(H_I1,y)

	def printfunc(self,t):
		print self.H0
		self.hamgen.update(t/self.norm)
		Int=self.hamgen.Int
		expm=LA.expm(1.0j*self.H0*t)
		iexpm=LA.expm(-1.0j*self.H0*t)
		H_I1=np.dot(expm,np.dot(Int,iexpm))
		
		#H=self.hamgen.update(0.5)
		return H_I1 
		
	

	def prop(self,track,i,j):	
		#Initial value: pure neutrino-0 state


		y0=self.b[i]
		x0=0.0
	
		xf=self.param.km*track.l
		self.norm=xf
		step=self.param.km*track.step

	
		self.r.set_initial_value(y0, x0)
		
		#Solve!
		dist=[]
		output=[]
		
		output.append(y0)



		while self.r.successful() and self.r.t <= xf:
			output.append(self.r.integrate(self.r.t+step))
			dist.append(self.r.t)
		
		amp=np.zeros(len(dist))


		for k in range(0,len(dist)):
			print "K:",k
			print "dist:",len(dist)
			print "output:",len(output)
			sch_out=np.dot(LA.expm(-1.0j*dist[k]*self.H0),output[k])
			ip=np.dot(self.b[j],sch_out)	
			amp[k]=np.absolute(ip)**2
	
		return amp


