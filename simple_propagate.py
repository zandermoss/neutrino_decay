import pickle
import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA
import MatrixOps as MO

import PhysConst as PC
import Splines  
import SUGen as SU
import Track
import random
import cmath
import math

import DeSolve
import HamGen
import PMNSGen
import SUGen


##A collection of functions useful for performing simulations which tie 
# together all the moving parts of NuSHEEP.
# the function currently used is AtmosphericNeutrinoOscillationProbability().


#Real parameters

#Earth Density and Electron Fraction Splines
splines = Splines.Spline()

resolution=10.0**4/6371.0



## A function used to calculate the amplitude for transition between two flavors.
# this function combines the different pieces of NuSHEEP, and performs a calculation
# of transition amplitude from initial_flavor to final_flavor.
# @param initial_flavor the initial neutrino flavor at the point of production.
# @param final_flavor the target flavor for detection at the end of the trajectory.
# @param energy the neutrino energy.
# @param theta the zenith angle.
# @param myparam a PhysConst object.
# @param pmnsgen a PMNSGen object specifying the flavor-mass mixing matrix.
# @param ugen a SUGen object specifying the flavor-decay mixing matrix.
# @param eig_dcy a vector of eigenvalues of the decay matrix.
# @return the transition amplitude.

def AtmosphericNeutrinoOscillationProbability(initial_flavor,final_flavor,
                           erange,theta,dcy_channels,tau,myparam,pmnsgen,nu_mass,phi_mass,regen):


	""" 
	Here, we have the option of sampling a random decay matrix, but currently
	we're specifying the decay parameters and feeding them in through the ugen 
	object. When it is time to search the parameter space by sampling, the first two
	lines would implement that sampling.
	"""

	#ugen=SU.SUGen(myparam)
	#ugen.sample_params()
	#Ug=ugen.matrix_gen()

	"""
	The Um matrix is generated from the PMNSGen object.
	""" 

	Um=pmnsgen.matrix_gen()


	"""
	The track object is instantiated, and the track length is calculated.
	"""

	track=Track.Track(myparam,resolution,erange,theta,False)
	track.calc_l()


	"""
	The hamgen object is instantiated. Here are examples of different propagation
	modes. The first line has both eig_dcy and splines arguments, so the propagation
	will take place in earth matter, and decay effects will be included. The second
	line has no decay effects, but matter is included. The third has decay but no 
	matter, and the 4th is a simple vacuum propagation.
	"""

	vhamgen=HamGen.HamGen(myparam,Um,tau,nu_mass,phi_mass,track,True,dcy_channels,splines,regen)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,splines)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,None)
	#vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,None)

	#prem solution

	p0 = MO.ProjMat(initial_flavor,myparam.numneu)
	pf = MO.ProjMat(final_flavor,myparam.numneu)

	rho0_mass=np.dot(Um.conj().T,np.dot(p0,Um))
	rhof_mass=np.dot(Um.conj().T,np.dot(pf,Um))



	rho0=np.zeros([len(erange),myparam.numneu,myparam.numneu],complex)
	for ei in range(0,len(erange)):
		rho0[ei,:,:] = rho0_mass	

	
	"""
	A DeSolve object is instantiated, with the freshly minted hamiltonian generator
	and PhysConst parameters as arguments.
	"""

	desolve= DeSolve.DeSolve(vhamgen,myparam)

	"""
	The transition ampliude from initial_flavor to final_flavor is calculated using
	the DeSolve prop function. This amplitude is then returned to the calling script.
	"""

	rhof=desolve.prop(track,rho0)

	amps= np.zeros([len(erange)])
	 
	for i in range(0,len(erange)):
		amps[i] = np.real(MO.Trace(np.dot(rhof[i],rhof_mass)))

	return amps

#def PropagationHistory(initial_flavor,final_flavor,
#                           energy,theta,myparam,pmnsgen,ugen,eig_dcy):
#
#
#	""" 
#	Here, we have the option of sampling a random decay matrix, but currently
#	we're specifying the decay parameters and feeding them in through the ugen 
#	object. When it is time to search the parameter space by sampling, the first two
#	lines would implement that sampling.
#	"""
#
#	#ugen=SU.SUGen(myparam)
#	#ugen.sample_params()
#	Ug=ugen.matrix_gen()
#
#	"""
#	The Um matrix is generated from the PMNSGen object.
#	""" 
#
#	Um=pmnsgen.matrix_gen()
#
#
#	"""
#	The track object is instantiated, and the track length is calculated.
#	"""
#
#	track=Track.Track(myparam,resolution,energy,theta,False)
#	track.calc_l()
#
#
#	"""
#	The hamgen object is instantiated. Here are examples of different propagation
#	modes. The first line has both eig_dcy and splines arguments, so the propagation
#	will take place in earth matter, and decay effects will be included. The second
#	line has no decay effects, but matter is included. The third has decay but no 
#	matter, and the 4th is a simple vacuum propagation.
#	"""
#
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,splines)
##	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,splines)
##	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,None)
#	#vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,None)
#
#	#prem solution
#	
#	"""
#	A DeSolve object is instantiated, with the freshly minted hamiltonian generator
#	and PhysConst parameters as arguments.
#	"""
#
#	desolve= DeSolve.DeSolve(vhamgen,myparam)
#
#	"""
#	The transition ampliude from initial_flavor to final_flavor is calculated using
#	the DeSolve prop function. This amplitude is then returned to the calling script.
#	"""
#
#	distances,amplitudes=desolve.prop_hist(track,initial_flavor,final_flavor)
#
#	return distances,amplitudes 
#
###Deprecated: loops over energy and zenith angle are now done from the calling script.
#
#def gridrun():
#
#	param=PC.PhysicsConstants()
#	print "SHEEP"
#	print param.GeV
#	print param.TeV
#	theta = np.linspace(param.PI/2.0, param.PI, 10)
#	energy = np.logspace(9, 12, 10)
#	print theta
#	print energy
#	prob=np.zeros((10,10))
#	for e in enumerate(energy):
#		for t in enumerate(theta):
#			print e,t
#			print "OWL"
#			prob[e[0],t[0]]=AtmosphericNeutrinoOscillationProbability(1,1,e[1],t[1],param)
#			
#
#	print prob
#	E, T = np.meshgrid(energy, theta)
#
#
#	plt.figure()
#	CS = plt.contourf(E,T,prob,cmap=plt.cm.jet)
#	plt.clabel(CS, inline=1, fontsize=10)
#	cbar = plt.colorbar(CS)
#	plt.show()	
#
#
#
##gridrun()
#
#
