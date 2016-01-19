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
import Splines

class DeSolve(object):

	def __init__(self,myhamgen,param):
		self.hamgen = myhamgen
		self.param=param

		#Set up the solver
		self.norm=0
		self.r=ode(self.func).set_integrator('zvode', method='adams',rtol=1e-10)

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
		
		"""
		dists=track.shell_intersection()
		dists=dists*self.param.EARTHRADIUS*self.param.km
		tracklen=track.l*self.param.km
		trackdel=4*track.delta*self.param.EARTHRADIUS*self.param.km
		if (track.intersections[0]==False)&(track.intersections[1]==False):
			output=self.r.integrate(xf)
		elif (track.intersections[0]==False)&(track.intersections[1]==True):
			output=self.r.integrate(dists[1]-trackdel)	
			self.r.set_initial_value(output,dists[1]+trackdel)
			output=self.r.integrate(tracklen-dists[1]-trackdel)
			self.r.set_initial_value(output,tracklen-dists[1]+trackdel)
			output=self.r.integrate(xf)
			print "THETA:",track.theta
			print "Sprime:",tracklen/self.param.km-dists[1]/self.param.km
			print "LEN:",tracklen/self.param.km
			print "DIST1:",dists[1]/self.param.km
			dspline=self.splines.GetEarth()
			dspace=np.linspace((dists[1]-trackdel)/(self.param.EARTHRADIUS*self.param.km),(dists[1]+trackdel)/(self.param.EARTHRADIUS*self.param.km),100)
			density=np.zeros(len(dspace))
			for j in range(0,len(density)):
				density[j]=dspline(track.r(dspace[j]))
		
			plt.plot(dspace,density)
			plt.show()
	
		elif (track.intersections[0]==True)&(track.intersections[1]==True):
			output=self.r.integrate(dists[1]-trackdel)	
			self.r.set_initial_value(output,dists[1]+trackdel)
			output=self.r.integrate(dists[0]-trackdel)	
			self.r.set_initial_value(output,dists[0]+trackdel)
			output=self.r.integrate(tracklen-dists[0]-trackdel)
			self.r.set_initial_value(output,tracklen-dists[0]+trackdel)
			output=self.r.integrate(tracklen-dists[1]-trackdel)
			self.r.set_initial_value(output,dists[1]+trackdel)
			output=self.r.integrate(xf)
<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes

		else:
			print "IMPOSSIBLE INTERSECTION PATTERN"

		"""

#		while self.r.successful() and self.r.t <= xf:
#			output.append(self.r.integrate(self.r.t+step))
#			dist.append(self.r.t)
		
	
		dist=[]
		output=[]
		while self.r.successful() and self.r.t <= xf:
			output.append(self.r.integrate(self.r.t+step))
			dist.append(self.r.t)
	
		amp=np.zeros(len(output))


		for k in range(0,len(output)):
			ip=np.dot(self.b[j],output[k])	
			amp[k]=np.absolute(ip)**2

		return amp

		#ip=np.dot(self.b[j],output)	
		#amp=np.absolute(ip)**2
	
		#return amp
