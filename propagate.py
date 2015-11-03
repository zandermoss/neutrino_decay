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
import HamGen 

#Real parameters
param = PC.PhysicsConstants()

hamgen=HamGen.HamGen(param)

eig_dcy=np.zeros(param.numneu)
for i in range(0,len(eig_dcy)):
	eig_dcy[i]=random.random()*1e-5

H=hamgen.gen(eig_dcy,1)

print "H:", H

asolve = ApproxSolve.ApproxSolve(H,param)
nsolve = NumSolve.NumSolve(H,param)


dist=np.arange(0,100000)

dist*=10

a_amp=np.zeros(len(dist))
n_amp=np.zeros(len(dist))
for i in range(0,len(dist)):
	a_amp[i] = asolve.P_ee(dist[i])
	n_amp[i]= nsolve.scalar_prop(dist[i],0,0)


#Plot oscillation amplitudes
fig, ax = plt.subplots()
ax.plot(dist,n_amp,'r-',label='P(e->e): Numerical')
ax.plot(dist,a_amp,'b-',label='P(e->e): Approximate')

ax.set_xlabel("Distance (A.U.)")
ax.set_ylabel("Oscillation Amplitude")
#ax.set_title("Evolution of 3 Flavors with a Random Hamiltonian")
ax.set_title("Comparison of Numerical Evolution to AdG Approximation")

#plt.xlim([0,5])

legend = ax.legend(loc='upper right', shadow=False)

plt.show()
