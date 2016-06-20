#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import random
import sys
#energy=np.logspace(10,12,10)
from pylab import rcParams
rcParams['figure.figsize'] = 14, 10
plt.rc('font', size=28)


import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.set_cmap(cmaps.inferno)


fig,ax=plt.subplots()		

f=np.load("1_1ev_2d_osc_dcy.npz")


e_edges=f['e']
t_edges=f['t']
probs=f['p']

print t_edges
#plt.rc('font', size=14)
plt.pcolor(t_edges,e_edges,probs,cmap="viridis",vmin=0,vmax=1)
plt.colorbar(label="Muon Antineutrino Survival Probability")

print t_edges.shape
print probs.shape

ax.set_ylim(e_edges[0],e_edges[-1])
ax.set_xlim(t_edges[0],t_edges[-1])

ax.set_yscale('log')
ax.set_xlabel(r'Cos($\theta_{Z}$)')
ax.set_ylabel('Energy (GeV)')
ax.set_title('Antineutrino Oscillation. Mixed Decay')
#ax.set_xlim([1e2,1e3])
#ax.set_yscale('log')


plt.show()



