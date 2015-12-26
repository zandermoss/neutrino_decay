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


#Oscillation Channel:
channel=[0,0]

#ODE solve methods
ode_methods=['BDF','Adams']


H=hamgen.gen(eig_dcy,param.GeV,matter)


#asolve = ApproxSolve.ApproxSolve(H,param)

dists=[]
amps=[]
deamps=[]

for pwr in range(4,7):
	res=10**pwr
	print res
	
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

	nsolve = NumSolve.NumSolve(H,param)
	n_amp=np.zeros(len(dist))
	for i in range(0,len(dist)):
	#	a_amp[i] = asolve.P_ee(dist[i])
		n_amp[i]= nsolve.scalar_prop(dist[i],channel[0],channel[1])
	dists.append(xdist)
	amps.append(n_amp)

	for method in ode_methods:
		
		desolve= DeSolve.DeSolve(H,param,method)
			
		d_amp=desolve.prop(dist,channel[0],channel[1])
		print "RAW", len(d_amp)
		d_amp=d_amp[0:len(xdist)]
		
		deamps.append(d_amp)
	

plt.style.use('ggplot')

fig, ax = plt.subplots()

colors=['r','g','b','m']

linestyles=['-','--',':']
index=0

for pwr in range(0,3):
	dlab=('Diagonalized: '+str(10**(pwr+4)))
	ax.plot(dists[pwr],amps[pwr],colors[pwr]+linestyles[0],label=dlab)
	for x in range(0,2):
		delab=('Numerical Method: '+ode_methods[x])
		ax.plot(dists[pwr],deamps[index],colors[pwr]+linestyles[x+1],label=delab)
		index+=1
	
ax.set_xlabel("Distance (km)")
ax.set_ylabel("Oscillation Amplitude")
#ax.set_title("Evolution of 3 Flavors with a Random Hamiltonian")
#ax.set_title("Comparison of Numerical Evolution to AdG Approximation")

# Now add the legend with some customizations.
legend = ax.legend(loc='upper right', shadow=True)

plt.xlim([0,param.EARTHRADIUS])
plt.show()

	
