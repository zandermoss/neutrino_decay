import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt



#---------------------
#Generate basis vectors
f=[]
for x in range(0,3):
	f.append(np.zeros(3))
	f[x][x]=1.0
#---------------------

#Generate hamiltonian
#note: uniform distribution over [(+/-)1,(+/-)1]
mr=2*np.random.rand(3,3)-1
mi=2*np.random.rand(3,3)-1
m=np.asarray(mr+1j*mi)

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

print time
print probs

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

plt.show()

