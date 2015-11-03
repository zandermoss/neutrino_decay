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

	def __init__(self,H,param):
		self.H = H
		self.param=param

		#Set up the solver
		self.r=ode(func).set_integrator('zvode', method='bdf')

		#Generate basis vectors
		self.b=[]
		for x in range(0,param.numneu):
			self.b.append(np.zeros(param.numneu))
			self.b[x][x]=1.0
		#---------------------
		
	

	def prop(dist,i,j):	
		#Initial value: pure neutrino-0 state
		y0=self.b[i]
		x0=0
	
		xf=dist[-1]
		step=dist[1]-dist[0]

	
		self.r.set_initial_value(y0, x0)
		
		#Solve!
		dist=[]
		output=[]
		
		while r.successful() and r.t < xf:
			output.append(r.integrate(r.t+step))
			dist.append(r.t)
	
		amp=np.zeros(len(output))


		for x in range(0,len(output)):
			ip=np.dot(self.b[j],output[x])	
			amp[x]=np.absolute(ip)**2
	
		return amp
