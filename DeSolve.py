import numpy as np

import scipy as sp
from scipy.integrate import ode
from scipy.integrate import complex_ode
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numpy import linalg as LA



import PhysConst as PC
import MinimalTools as MT
import random
import cmath
import math
import Splines

class DeSolve(object):

	def __init__(self,myhamgen,param):
		self.hamgen = myhamgen
		self.param=param

		#Set up the solver
		self.norm=0
		#self.r=ode(self.func).set_integrator('zvode', method='bdf',rtol=1e-6)
		#self.r=complex_ode(self.func).set_integrator('dopri5',rtol=1e-6)
		self.r=complex_ode(self.func).set_integrator('dopri5',rtol=1e-6,nsteps=100000)


		#Generate basis vectors
		self.b=[]
		for x in range(0,param.numneu):
			self.b.append(np.zeros(param.numneu))
			self.b[x][x]=1.0

		self.splines=Splines.Spline()

		#---------------------

	#Time independent hamiltonian,
	#define as matrix product for DE solver:
	def func(self,t,y):
		H=self.hamgen.update(t/self.norm)
		#H=self.hamgen.update(0.5)
		#H=self.hamgen.H

		return -1j*np.dot(H,y)

	def printfunc(self,t):
		H=self.hamgen.update(t/self.norm)
		#H=self.hamgen.update(0.5)
		print H
		return H 
		
	

	def prop(self,track,i,j):	
		#Initial value: pure neutrino-0 state


		y0=self.b[i]
		x0=0.0
	

		xf=self.param.km*track.l
		#print "BASELINE:", self.param.km*track.l
		self.norm=xf
		step=self.param.km*track.step

		self.r.set_initial_value(y0, x0)
		output=self.r.integrate(xf)

		ip=np.dot(self.b[j],output)	
		amp=np.absolute(ip)**2
	
		return amp

		"""
		outarr = []
		distarr=[]
		while self.r.successful() and self.r.t < xf:
			distarr.append(self.r.t/self.param.km)
			outarr.append(self.r.integrate(self.r.t+step))


		np_out=np.asarray(outarr)
		np_dist=np.asarray(distarr)
		return np_dist,np_out
		"""

