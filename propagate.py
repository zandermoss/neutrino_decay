import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import Splines  
import SUGen as SU
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
spman = Splines.Spline()
dspline=spman.GetEarth()
yespline=spman.GetYe()

eig_dcy=np.zeros(param.numneu)
dcy_ord=-18.0

eig_dcy[0]=3.1989727405321533*(10.0)**dcy_ord
eig_dcy[1]=9.283607843296119*(10.0)**dcy_ord
eig_dcy[2]=6.567215979512332*(10.0)**dcy_ord

#for i in range(0,len(eig_dcy)):
#	eig_dcy[i]=random.random()*1e-6
#	#eig_dcy[i]=0


print "NUMNEU:",param.numneu
shamgen=HamGen.HamGen(param,eig_dcy,3.90863690777e-13)
vhamgen=HamGen.HamGen(param,eig_dcy,dspline,yespline)


osc_test=True
#osc_test=False
#matter=False
#matter=True





#Oscillation Channel
channel=[1,1]

nruns=100

for x in range(0,nruns):

	ugen=SU.SUGen(param)
	ugen.sample_params()
	Ug=ugen.matrix_gen()



	res=10**4
	#Progress display
	if (x+1)%((nruns)/10)==0:
		print "Done: ",x+1,"/",nruns
	

	Hv=vhamgen.gen(param.TeV,Ug)
	Hs=vhamgen.update(0.5)
	
	
	#asolve = ApproxSolve.ApproxSolve(H,param)
	nsolve = NumSolve.NumSolve(Hs,param)
	desolve= DeSolve.DeSolve(vhamgen,param)
	
	#xdist=np.arange(0,param.EARTHRADIUS)
	xdist=np.linspace(0,param.EARTHRADIUS,res)
	dist = xdist*param.km #km to eV	

	
	a_amp=np.zeros(len(dist))
	n_amp=np.zeros(len(dist))
	for i in range(0,len(dist)):
	#	a_amp[i] = asolve.P_ee(dist[i])
		n_amp[i]= nsolve.scalar_prop(dist[i],channel[0],channel[1])
	
	d_amp=desolve.prop(dist,channel[0],channel[1])
	print "RAW", len(d_amp)
	d_amp=d_amp[0:len(xdist)]
	print "SHAM",Hs

	if osc_test:
		fig, ax = plt.subplots()

		ax.plot(xdist,n_amp,'r-',label='P(mu->mu): Diagonalized')
		#ax.plot(xdist,a_amp,'b-',label='P(e->e): Approximate')
		ax.plot(xdist,d_amp,'g--',label='P(mu->mu): Numerical')
		
		ax.set_xlabel("Distance (km)")
		ax.set_ylabel("Oscillation Amplitude")
		ax.set_title("Diagonal Method with Static Matter vs. Numerical with Full Matter.")
		#ax.set_title("Comparison of Numerical Evolution to AdG Approximation")
		

		legend = plt.legend(loc='upper right', shadow=False, fontsize='x-large')

		# Put a nicer background color on the legend.
		#legend.get_frame().set_facecolor('#00FFCC')

		#Look for max difference		
		diff=np.absolute(n_amp-d_amp)
		print
		print "DIFF", np.max(diff)
		print


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
	
	
