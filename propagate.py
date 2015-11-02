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

H=hamgen.gen(eig_dcy)

asolve = ApproxSolve.ApproxSolve(H,param)
nsolve = NumSolve.NumSolve(H,param)


dist=np.arange(0,100000)

n_amp=nsolve.scalar_prop(dist,0,0)
a_amp=asolve.P_ee(dist)




"""
dist*=10
for x in range(0,len(dist)):
	 
	p=prop(w,N,Ni,dist[x],Ht)

plt.plot(dist,amp,'b-')
#plt.show()


dist=np.arange(0,100000)
amp=np.zeros(100000)

dist*=10
for x in range(0,len(dist)):
	amp[x]=P_ee(dist[x],approx_param) 

plt.plot(dist,amp,'r-')
plt.show()
"""

"""
#Plot oscillation amplitudes
fig, ax = plt.subplots()
ax.plot(time,probs[0],'r-',label='P(1->1)')
ax.plot(time,probs[1],'g-',label='P(1->2)')
ax.plot(time,probs[2],'b-',label='P(1->3)')

ax.set_xlabel("Distance (A.U.)")
ax.set_ylabel("Oscillation Amplitude")
ax.set_title("Evolution of 3 Flavors with a Random Hamiltonian")

plt.xlim([0,5])

legend = ax.legend(loc='upper right', shadow=False)

#plt.show()
"""
