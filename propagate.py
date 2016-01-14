import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import Splines  
import SUGen as SU
import Track
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

#Earth Density and Electron Fraction Splines
splines = Splines.Spline()

eig_dcy=np.zeros(param.numneu)
dcy_ord=-18.0

eig_dcy[0]=3.1989727405321533*(10.0)**dcy_ord
eig_dcy[1]=9.283607843296119*(10.0)**dcy_ord
eig_dcy[2]=6.567215979512332*(10.0)**dcy_ord

#for i in range(0,len(eig_dcy)):
#	eig_dcy[i]=random.random()*1e-6
#	#eig_dcy[i]=0


print "NUMNEU:",param.numneu


#osc_test=True
osc_test=False
#matter=False
#matter=True



resolution=10.0**4/param.EARTHRADIUS

#Oscillation Channel
channel=[1,1]

nruns=100

maxdiffs=np.zeros(nruns)
argdiffs=np.zeros(nruns)
thetas=np.zeros(nruns)

loopth=np.linspace(0,3.141592,nruns)
fname="data/run"

for x in range(0,nruns):


	ugen=SU.SUGen(param)
	ugen.sample_params()
	Ug=ugen.matrix_gen()

	energy=0.01+random.random()*0.99

	#track=Track.Track(param,resolution,energy*param.TeV,True)
	track=Track.Track(param,resolution,param.TeV,False)
	track.theta=loopth[x]
	track.calc_l(track.theta)
	thetas[x]=track.theta

	#shamgen=HamGen.HamGen(param,track,eig_dcy,3.90863690777e-13)
	shamgen=HamGen.HamGen(param,Ug,track,eig_dcy)
	vhamgen=HamGen.HamGen(param,Ug,track,eig_dcy,splines)

	#Progress display
	if (x+1)%((nruns)/10)==0:
		print "Done: ",x+1,"/",nruns
	

	Hn=shamgen.H
	Hs=vhamgen.update(0.5)
	
	#asolve = ApproxSolve.ApproxSolve(H,param)
	nsolve = NumSolve.NumSolve(Hs,param)
	nnsolve = NumSolve.NumSolve(Hn,param)
	desolve= DeSolve.DeSolve(vhamgen,param)
	

	

	d_amp=desolve.prop(track,channel[0],channel[1])
	#print "NSTEPS: ", len(d_amp)


	#print "RAW", len(d_amp)
	#print "SHAM",Hs
	xplot=np.linspace(0,track.l,len(d_amp),endpoint=True)
	#print "LEN", track.l
	dist=param.km*np.linspace(0,track.l,len(d_amp),endpoint=True)
	n_amp=np.zeros(len(dist))
	nn_amp=np.zeros(len(dist))
	for i in range(0,len(dist)):
		n_amp[i]= nsolve.scalar_prop(dist[i],channel[0],channel[1])
		nn_amp[i]= nnsolve.scalar_prop(dist[i],channel[0],channel[1])

	diff=np.absolute(n_amp-d_amp)
	maxdiffs[x]=np.max(diff)
	argdiffs[x]=np.argmax(diff)

	#np.savez(fname+str(x),resolution=resolution,energy=energy,theta=track.theta,dist=xplot,numerical=d_amp,static=n_amp,vacuum=nn_amp)

	if osc_test:
		fig, ax = plt.subplots()
		plotx=np.linspace(0,track.l,len(d_amp),endpoint=True)

		ax.plot(plotx,n_amp,'r-',label='P(mu->mu): Diagonalized')
		ax.plot(plotx,nn_amp,'b--',label='P(mu->mu): Diagonalized, nomatter')
		#ax.plot(xdist,a_amp,'b-',label='P(e->e): Approximate')
		print "THETA", track.theta/3.1415
		print "Length", track.l
		ax.plot(plotx,d_amp,'g--',label='P(mu->mu): Numerical')
		
		ax.set_xlabel("Distance (km)")
		ax.set_ylabel("Oscillation Amplitude")
		ax.set_title("Diagonal Method with Static Matter vs. Numerical with Full Matter.")
		#ax.set_title("Comparison of Numerical Evolution to AdG Approximation")
		

		legend = plt.legend(loc='upper right', shadow=False, fontsize='x-large')

		# Put a nicer background color on the legend.
		#legend.get_frame().set_facecolor('#00FFCC')
		'''
		#Look for max difference		
		diff=np.absolute(n_amp-d_amp)
		#plt.cla()
		#ax.plot(plotx,diff,'r-',label='Diff')
		print
		print "DIFF", np.max(diff)
		print "ARGDIFF", np.argmax(diff)
		print "LENDIFF:", len(diff)
		print
		'''

		plt.xlim([0,track.l])
		plt.show()
		break		


np.save("linear_thetaout.npy",thetas)
np.save("linear_diffsout.npy",maxdiffs)
fig, ax = plt.subplots()
ax.plot(thetas,maxdiffs,'o')
ax.set_xlabel("theta (radians)")
ax.set_ylabel("Maximum Probability Split")
ax.set_title("Maximum Difference between Numerical and Diagonal Probabilities")
plt.show()

"""
if x==0:
	xrun=xdist
	yrun=n_amp
else:
	xrun=np.concatenate((xrun,xdist))
	yrun=np.concatenate((yrun,n_amp))

"""

"""
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
	
"""	
