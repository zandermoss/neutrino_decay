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

param.numneu=3

print "NUMNEU:",param.numneu
hamgen=HamGen.HamGen(param)

eig_dcy=np.zeros(param.numneu)

osc_test=True
#osc_test=False
#matter=False
matter=True

eig_dcy[0]=3.1989727405321533e-07
eig_dcy[1]=9.283607843296119e-07
eig_dcy[2]=6.567215979512332e-07

#for i in range(0,len(eig_dcy)):
#	eig_dcy[i]=random.random()*1e-6
#	#eig_dcy[i]=0



nruns=100

for x in range(0,nruns):

	#Progress display
	if (x+1)%((nruns)/10)==0:
		print "Done: ",x+1,"/",nruns
	

	H=hamgen.gen(eig_dcy,param.GeV,matter)
	
	
	#asolve = ApproxSolve.ApproxSolve(H,param)
	nsolve = NumSolve.NumSolve(H,param)
	desolve= DeSolve.DeSolve(H,param)
	
	#xdist=np.arange(0,param.EARTHRADIUS)
	xdist=np.linspace(0,param.EARTHRADIUS,res)
	#xdist=np.divide(xdist,10.0)	
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
		n_amp[i]= nsolve.scalar_prop(dist[i],0,0)
	
	d_amp=desolve.prop(dist,0,0)
	#print "RAW", len(d_amp)
	#d_amp=d_amp[0:len(xdist)]

	if osc_test:
		fig, ax = plt.subplots()

		ax.plot(xdist,n_amp,'r-',label='P(e->e): Diagonalized')
		#ax.plot(xdist,a_amp,'b-',label='P(e->e): Approximate')
		ax.plot(xdist,d_amp,'g--',label='P(e->e): Numerical')
		
		ax.set_xlabel("Distance (km)")
		ax.set_ylabel("Oscillation Amplitude")
		#ax.set_title("Evolution of 3 Flavors with a Random Hamiltonian")
		#ax.set_title("Comparison of Numerical Evolution to AdG Approximation")
		
		plt.xlim([0,param.EARTHRADIUS])
		plt.show()
		break		

	if x==0:
		xrun=xdist
		yrun=n_amp
	else:
		xrun=np.concatenate((xrun,xdist))
		yrun=np.concatenate((yrun,n_amp))


if (osc_test==False):	
	nbinsx = 1000
	nbinsy = 100
	H, xedges, yedges = np.histogram2d(xrun,yrun,bins=(nbinsx,nbinsy))
	
	# H needs to be rotated and flipped
	H = np.rot90(H)
	H = np.flipud(H)
	
	 
	# Mask zeros
	Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
	# Plot 2D histogram using pcolor
	fig2 = plt.figure()
	fig2.patch.set_facecolor('white')
	plt.pcolormesh(xedges,yedges,Hmasked)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.xlabel("Distance (km)")
	plt.ylabel("Oscillation Amplitude")
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Counts') 
	plt.xlim([0,param.EARTHRADIUS])
	plt.show()
	
	
