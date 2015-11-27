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

import ApproxSolve
import NumSolve
import DeSolve
import HamGen 

#Real parameters
param = PC.PhysicsConstants()

hamgen=HamGen.HamGen(param)

eig_dcy=np.zeros(param.numneu)
for i in range(0,len(eig_dcy)):
	eig_dcy[i]=random.random()*1e-5
	#eig_dcy[i]=0

for x in range(1,10):
	H=hamgen.gen(eig_dcy,param.MeV)
	
	print "H:", H
	
	asolve = ApproxSolve.ApproxSolve(H,param)
	nsolve = NumSolve.NumSolve(H,param)
	desolve= DeSolve.DeSolve(H,param)
	
	xdist=np.arange(0,1000)
	
	dist=xdist
	dist=dist*1000 #km to m
	
	#we can just convert the distance to MKS, as it is the only
	#parameter
	
	dist=dist*1/(param.hbar*param.sol*param.Joule)
	
	#converting into time*MKS conversion factors, so we'll end up with:
	#J*s/M*M*1/hbar in the exponent. Perfect!
	
	#we should now be looking at effective propagation in meters
	#print xdist
	
	a_amp=np.zeros(len(dist))
	n_amp=np.zeros(len(dist))
	for i in range(0,len(dist)):
	#	a_amp[i] = asolve.P_ee(dist[i])
		n_amp[i]= nsolve.scalar_prop(dist[i],0,1)
	
	
	d_amp=desolve.prop(dist,0,0)
	
	
	
	#Plot oscillation amplitudes
	fig, ax = plt.subplots()
	ax.plot(xdist,n_amp,'r-',label='P(e->e): Diagonalized')
	#ax.plot(xdist,a_amp,'b-',label='P(e->e): Approximate')
	#ax.plot(xdist,d_amp,'g-',label='P(e->e): Numerical')
	
	ax.set_xlabel("Distance (Meters)")
	ax.set_ylabel("Oscillation Amplitude")
	#ax.set_title("Evolution of 3 Flavors with a Random Hamiltonian")
	#ax.set_title("Comparison of Numerical Evolution to AdG Approximation")
	
	plt.xlim([0,300])
	
	legend = ax.legend(loc='upper right', shadow=False)
	
	plt.show()
