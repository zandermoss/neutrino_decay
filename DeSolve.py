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

		#Set up the solver
		self.norm=0
		self.r=ode(self.func).set_integrator('zvode', method='adams',rtol=1e-8)

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
		self.norm=xf
		step=self.param.km*track.step

	
		self.r.set_initial_value(y0, x0)
		

		dists=track.shell_intersection()

		if (track.intersections[0]==False)&(track.intersections[1]==False):
			output=self.r.integrate(xf)
		elif (track.intersections[0]==False)&(track.intersections[1]==True):
			output=self.r.integrate(dists[1]-track.delta)	
			self.r.set_initial_value(output,dists[1]+track.delta)
			output=self.r.integrate(track.l-dists[1]-track.delta)
			self.r.set_initial_value(output,dists[1]+track.delta)
			output=self.r.integrate(xf)
		elif (track.intersections[0]==True)&(track.intersections[1]==True):
			output=self.r.integrate(dists[1]-track.delta)	
			self.r.set_initial_value(output,dists[1]+track.delta)
			output=self.r.integrate(dists[0]-track.delta)	
			self.r.set_initial_value(output,dists[0]+track.delta)
			output=self.r.integrate(track.l-dists[0]-track.delta)
			self.r.set_initial_value(output,track.l-dists[0]+track.delta)
			output=self.r.integrate(track.l-dists[1]-track.delta)
			self.r.set_initial_value(output,dists[1]+track.delta)
			output=self.r.integrate(xf)



#		while self.r.successful() and self.r.t <= xf:
#			output.append(self.r.integrate(self.r.t+step))
#			dist.append(self.r.t)
		
	


		ip=np.dot(self.b[j],output)	
		amp=np.absolute(ip)**2
	
		return amp
