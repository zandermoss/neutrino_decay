import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.interpolate import interp1d
import numpy.fft as FFT

from scipy.signal import savgol_filter
from scipy.ndimage.filters import gaussian_filter1d as gausfilt
import PhysConst as PC
import MinimalTools as MT
import random
import cmath
import math

import EarthRadial as ER

er=ER.EarthRadial()
param=PC.PhysicsConstants()

nspl=10000

dist=np.linspace(0,1,nspl)
newdist=np.linspace(0,1,100000)
kmdist=np.multiply(dist,param.EARTHRADIUS)
newkmdist=np.multiply(newdist,param.EARTHRADIUS)

dense=np.zeros(len(dist))

for x in range(0,len(dist)):
	dense[x]=er.rdensity(dist[x])

diff=np.diff(dense)
diffpad=np.zeros(len(dense))
diffpad[:len(diff)]=diff


pre=int(float(nspl)*0.1)
print pre

diffmean=np.mean(diffpad[:pre])
print diffmean

t=-0.025

pairs=[]
starts=[]
stops=[]
thresh=False
for x in range(0,len(diffpad)):
	if ((diffpad[x]<=t)&(thresh==False)):
		thresh=True
		start=x-6	
		starts.append(x-6)
		print "start"
	elif ((diffpad[x]>=t)&(thresh==True)):
		if x+5>=len(dist)-1:
			jump=(len(dist)-1)-x
		else:
			jump=5
		thresh=False
		stop=x+jump
		stops.append(x+jump)
		pairs.append([start,stop])	
		print "stop"
	
#Now assemble the splining ranger
splinecount=1000

local_nbins=int(float(splinecount)*0.25)
global_nbins=int(splinecount-local_nbins)

dynadist=np.linspace(0,1,global_nbins)

peak_nbins=int(float(local_nbins)/float(len(pairs)))

print
print #---------------#
print "Spline Count:",splinecount
print "Global NBins:",global_nbins
print "Peak NBins:",peak_nbins
print "NPeaks:", len(pairs)
print #---------------#
print

for pair in pairs:
	local_dist=np.linspace(dist[pair[0]],dist[pair[1]],peak_nbins,endpoint=True)
	dynadist=np.append(dynadist,local_dist)
	

dynadist=np.unique(dynadist)


print dynadist

smooth=gausfilt(dense,10)

smoothfunc= interp1d(dist,smooth,kind='linear')

dynadense=np.zeros(len(dynadist))
for x in range(0,len(dynadist)):
	dynadense[x]=smoothfunc(dynadist[x])
	

smalldense=np.zeros(splinecount)
smalldist=np.linspace(0,1,splinecount)

for x in range(0,len(smalldist)):
	smalldense[x]=smoothfunc(smalldist[x])

dynafunc = interp1d(dynadist, dynadense, kind='cubic')
func = interp1d(smalldist, smalldense, kind='cubic')

#linedense=func(dynadist)

#func3 = interp1d(dynadist, linedense, kind='cubic')
	
fig, ax = plt.subplots()


avgdense=np.zeros(len(dense))
'''
for x in range(0,6):
	avgdense[x]=dense[x]
for x in range(0,len(dist)-10):
	avgdense[x]=np.mean(dense[x-5:x+5])

for x in range(len(dist)-5,len(dist)):
	avgdense[x]=dense[x]
'''
ax.plot(dist,smooth,'b-')


#ax.plot(dist,diffpad,'--')
#ax.plot(dist,dense,'-')
#ax.plot(dist,avgdense,'-')
#ax.plot(dynadist,np.zeros(len(dynadist)),'mo')
ax.plot(newdist,dynafunc(newdist),'r--')
ax.plot(newdist,func(newdist),'g--')
#ax.plot(dynadist,func(dynadist),'go')
#ax.plot(dist[starts],diffpad[starts],'go')
#ax.plot(dist[stops],diffpad[stops],'ro')

print len(dense)
print len(np.diff(dense))
#ax.plot(newkmdist,f(newdist),'-')




"""
fourier=FFT.rfft(dense)
nfreq=np.arange(0,len(fourier))
ffreq=FFT.rfftfreq(dense.size,1)


rdense=FFT.irfft(fourier)
#f = interp1d(dist, dense, kind='quartic')

fig, ax = plt.subplots()

#ax.plot(kmdist,dense,'-')
#ax.plot(newkmdist,f(newdist),':')

#ax.plot(nfreq,fourier,"b-")
ax.plot(dist,rdense,"b-")
fourier[:100]=0
#ax.plot(nfreq,fourier,"r-")
nrdense=FFT.irfft(fourier)
ax.plot(dist,nrdense,"r-")
"""
#ax.set_xlabel("Distance (km)")
#ax.set_ylabel("Electron Density")
#ax.set_title("Evolution of 3 Flavors with a Random Hamiltonian")
#ax.set_title("Comparison of Numerical Evolution to AdG Approximation")

#plt.xlim([0,param.EARTHRADIUS])
plt.show()
