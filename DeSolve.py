import numpy as np

import scipy as sp
from scipy.integrate import ode
from scipy.integrate import complex_ode
from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from numpy import linalg as LA
import MatrixOps as MO
import sys


import PhysConst as PC
import random
import cmath
import math
import Splines



## DeSolve implements a complex ode solver algorithm from scipy.integrate
# This class is called once all neutrino model parameters and the neutrino track 
# are set up. The constructor function instantiates a complex_ode solver object, 
# and sets integration parameters. The prop() function evolves an initial pure
# neutrino flavor through the earth, making calls to the func() function, which
# fetches an updated hamiltonian at each distance step, according to the geometry 
# information from the track object passed to prop().




class DeSolve(object):


	## The Constructor
	# accepts a PhysConst object, as well as a HamGen object, and assigns 
	# local variables to them. The integrator object is initialized with a specified
	# relative error tolerance and maximum number of steps. Finally, flavor basis 
	# vectors are defined for calculating oscillation amplitudes in prop().
	# @par myhamgen a HamGen object to be used for hamiltonian updating over the trajectory of the neutrino.
	# @par param a PhysConst object to be used by the solver.

	def __init__(self,myhamgen,param):
		self.hamgen = myhamgen
		self.param=param
		self.rhoshapes=[0,0,0]
		self.shapeprod=0.0
		self.track = self.hamgen.track 
		self.theta = self.track.theta

		self.length = self.hamgen.track.l
		self.res = 1000
		self.xsteps = np.linspace(0,1.0,self.res)
		self.yevals  = np.zeros(len(self.xsteps))
		self.dvals  = np.zeros(len(self.xsteps))
		for i in range(0,len(self.xsteps)):
			self.yevals[i] = self.hamgen.yespline(self.xsteps[i])
			self.dvals[i] = self.hamgen.dspline(self.xsteps[i])
	

		#Set up the solver
		self.norm=0
		#self.r=ode(self.func).set_integrator('zvode', method='bdf',rtol=1e-6)
		#self.r=complex_ode(self.func).set_integrator('dopri5',rtol=1e-6)
		self.r=complex_ode(self.func).set_integrator('dopri5',rtol=1e-6,nsteps=100000)

		if (param.neutype=='neutrino'):
			self.anti=False	
		elif (param.neutype=='antineutrino'):
			self.anti=True
		else:
			print "BAD NEUTRINO TYPE"
		"""
		self.projmats=[]
		for i in range(0,4):
			mat=np.zeros((4,4))
			mat[i,i]=1
			print mat
			self.projmats.append(mat)
		"""
		#---------------------


	## Updating function for the hamiltonian as the neutrino propagates through the earth.
	# Calls hamgen.update() with the progress along the total track length and current neutrino state. Returns the product of the updated hamiltonian and the neutrino state times a phase. This constitutes the RHS of the schrodinger equation (divided by i. hbar is 1 so the LHS is simply a time derivative, for the convenience of the solver).
	#par t the distance (time) travelled along the track.
	#par y deprecated!
	#return the RHS of the schrodinger equation divided by i.
	def func(self,t,rho):
		#H=self.hamgen.update(0.5)
		#H=self.hamgen.H
		rho.shape = self.rhoshapes

		r = MO.r(t, self.norm, self.theta)
		#r_arg = MO.find_nearest(self.xsteps,r)


		ye =self.hamgen.yespline(r)
		#ye = 0.5
		#ye =self.yevals[r_arg]
		density=self.hamgen.dspline(r) #g/cm^3
		#density=3.0 #g/cm^3
		#density =self.dvals[r_arg]


		H0 = self.hamgen.Hstart		
		#H0 = np.zeros(self.rhoshapes, np.complex128)
		if (self.hamgen.matter==True):
			H0=MO.H0_Update(ye, density, self.param.numneu, self.hamgen.Um, self.anti, self.track.erange, self.hamgen.Hstart)

		R = np.zeros(self.rhoshapes, np.complex128)
		if (self.hamgen.regen==True):
			R = MO.Regen(rho, self.track.erange, self.param.numneu, self.hamgen.dcy_channels, self.hamgen.pstar, self.hamgen.m_nu, self.hamgen.m_phi, self.hamgen.tau)
		"""
		prstr=""
		for i in range(0,4):
			prstr+=(str(MO.PyTrace(np.dot(rho[5,:,:],self.projmats[i])))+" ")
		print prstr
		"""
		drho_dt = MO.DrhoDt(rho, H0, R, self.hamgen.Gamma)	
		drho_dt.shape = self.shapeprod

		return drho_dt


	## Propagates a neutrino of flavor i along a track defined by the track argument. Calculates amplitude for oscillation into the j flavor.
	# Calculates the total track length from the track object (as well as the step resolution). Steps the initial neutrino state forward along the track, calling func() for hamiltonian updates at each stage. At the end of the evolution, the complex amplitude squared of the inner product of the final state with the jth flavor state is calculated and returned.
	# @par track a Track object with the track geometry.
	# @par i the initial neutrino flavor state.
	# @par j the target neutrino flavor state.
	# @return the amplitude for transition from the ith to jth flavor state.

	def prop(self,track,rho0):	
		self.shapeprod = (rho0.shape[0])*(rho0.shape[1])*(rho0.shape[2])
		self.rhoshapes[0] = rho0.shape[0]
		self.rhoshapes[1] = rho0.shape[1]
		self.rhoshapes[2] = rho0.shape[2]

		rho0.shape = self.shapeprod

		x0=0.0
	
		xf=self.param.km*track.l
		#print "BASELINE:", self.param.km*track.l
		self.norm=xf
		step=self.param.km*track.step

		self.r.set_initial_value(rho0, x0)
		rhof=self.r.integrate(xf)

		if (self.r.successful()!=True):
			sys.exit("Integrator Failed. Check for stiffness UserWarning above.")

		rhof.shape = self.rhoshapes
		return rhof

	## Propagates a neutrino of flavor i along a track defined by the track argument. Returns the entire history of propagation (distances, transition amplitudes to flavor j).
	# Calculates the total track length from the track object (as well as the step resolution). Steps the initial neutrino state forward along the track, calling func() for hamiltonian updates at each stage. At the end of the evolution, a vector of the complex amplitude squared of the inner product of each state in the iterative integration process with the jth flavor state is calculated and returned, along with a vector of relevant distances (distance at each propagation step).
	# @par track a Track object with the track geometry.
	# @par i the initial neutrino flavor state.
	# @par j the target neutrino flavor state.
	# @return (distances, amplitudes) a list of lists. distances is a list of the distances corresponding to each integration step, and amplitude is the corresponding amp-squared inner product of that intermediate state with the jth flavor state.
#
#	def prop_hist(self,track,i,j):	
#		#Initial value: pure neutrino-0 state
#
#
#		y0=self.b[i]
#		x0=0.0
#	
#
#		xf=self.param.km*track.l
#		#print "BASELINE:", self.param.km*track.l
#		self.norm=xf
#		step=self.param.km*track.step
#
#		self.r.set_initial_value(y0, x0)
#
#	
#		outarr = []
#		distarr = []
#		while self.r.successful() and self.r.t < xf:
#			distarr.append(self.r.t/self.param.km)
#			complex_state_vector=self.r.integrate(self.r.t+step)
#			ip=np.dot(self.b[j],complex_state_vector)	
#			amp=np.absolute(ip)**2
#			outarr.append(amp)
#
#		np_out=np.asarray(outarr)
#		np_dist=np.asarray(distarr)
#		return np_dist,np_out
#
