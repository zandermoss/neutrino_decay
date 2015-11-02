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

#Real parameters
param = PC.PhysicsConstants()

#---------------------
#Generate basis vectors
f=[]
for x in range(0,param.numneu):
	f.append(np.zeros(param.numneu))
	f[x][x]=1.0
#---------------------

#Generate ham!

#Randomized parameters to generate conjugation matrices
#We will work in the flavor basis, so Um and Ug map from 
#the mass basis to the flavor basis and the decay basis 
#to the flavor basis, respectively
Ugen=PC.PhysicsConstants()

#Generate conj matrices
#use known mixing parameters for M
Um=MT.calcU(param)

Ugen.randomize_trig()
Ug=MT.calcU(Ugen)

#Fill in mass and decay eigenvalues
Md=np.zeros([param.numneu,param.numneu],complex)
Gd=np.zeros([param.numneu,param.numneu],complex)

for i in range(0,param.numneu+0):

	print "m2:", param.dm2[1,i+1] #FIXME: mass constants-> Add in energy dependence!
	Md[i,i]= param.dm2[1,i+1] #FIXME: mass constants-> Add in energy dependence!
	#Md[i,i]= 1 #FIXME: mass constants-> Add in energy dependence!
	Gd[i,i]= random.random()*1e-5 #FIXME


print Md

M= np.dot(Um,np.dot(Md,Um.conj().T))
G= np.dot(Ug,np.dot(Gd,Ug.conj().T))

#Assemble Hamiltonian
H=np.zeros([param.numneu,param.numneu],complex)

H= M -1j*G


print H



dist=np.arange(0,100000)
amp=np.zeros(100000)

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


m2d=np.dot(Ug.conj().T,Um)

Gm = np.dot(m2d.conj().T,np.dot(Gd,m2d))

#Let's just look at the survival probability for electrons
"""
def P_ee(Gm,x):
	
	B=np.zeros(2)
	B[0]=Gm[0,0]
	B[1]=Gm[1,1]
	
	b=Gm[1,2]
	
	delta=Md[1,1]-Md[0,0]
	
	z=-2j*b/delta
	

#	P_ee = math.exp(-2*B[0]*x)*math.cos(




#Generate hamiltonian
#note: uniform distribution over [(+/-)1,(+/-)1]
mr=2*np.random.rand(3,3)-1
mi=2*np.random.rand(3,3)-1
m=np.asarray(mr+1j*mi)

print "EIGENVALUES: "
w,v =  LA.eig(m)
print w

#Time independent hamiltonian,
#define as matrix product for DE solver:
def func(t,y):
	return np.dot(m,y) 


#Initial value: pure neutrino-0 state
y0=f[0]
t0=0

#Set up the solver
r=ode(func).set_integrator('zvode', method='bdf')
r.set_initial_value(y0, t0)

#Determine integration range
t1=5
dt=0.01


#Solve!
time=[]
output=[]

while r.successful() and r.t < t1:
	output.append(r.integrate(r.t+dt))
	time.append(r.t)


#Calculate oscillation amplitudes
probs=np.empty([3,len(output)],dtype=complex)

for i in range(0,probs.shape[0]):
	for j in range(0,probs.shape[1]):
		probs[i,j] = np.absolute(np.dot(output[j],f[i]))**2

#print time
#print probs

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
