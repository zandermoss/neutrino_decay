#! /usr/bin/python
import Verosim as VS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
from numpy import sqrt
from pylab import rcParams

rcParams['figure.figsize'] = 12, 10

plt.rc('font', size=28)


cuts=np.zeros(2)
cuts[0]=6
cuts[1]=50-16


f = np.load("exp_and_edge.npz")
print f.files

eprox=f['ep']
cosz=f['cos']


eprox_cut=eprox[cuts[0]:cuts[1]]

null=np.load("null_oscillation_simp2.npz")


cut_exp=null['exp'][cuts[0]:cuts[1],:]
cut_dat=null['dat'][cuts[0]:cuts[1],:]

exp_cosz=np.sum(cut_exp,axis=0)
exp_eprox=np.sum(cut_exp,axis=1)

dat_cosz=np.sum(cut_dat,axis=0)
err_cosz=np.sqrt(dat_cosz)
dat_eprox=np.sum(cut_dat,axis=1)
err_eprox=np.sqrt(dat_eprox)


bins_center_cosz = np.asarray([ (cosz[i]+cosz[i+1])/2.0 for i in range(len(cosz)-1) ])
bins_center_eprox = np.asarray([ (eprox[i]+eprox[i+1])/2.0 for i in range(len(eprox)-1) ])

bins_center_eprox_cut=bins_center_eprox[cuts[0]:cuts[1]]

fix,ax=plt.subplots()

#n,b,p=ax.hist(bins_center_cosz, bins=cosz, weights=exp_cosz, color='Green',histtype="step",label="Expectation",lw=1)

#ax.errorbar(bins_center_cosz,dat_cosz,yerr=err_cosz,color="black",fmt="o",label="Data")

n,b,p=ax.hist(bins_center_eprox_cut, bins=eprox_cut, weights=exp_eprox, color='Green',histtype="step",label="Expectation",lw=2)

ax.errorbar(bins_center_eprox_cut,dat_eprox,yerr=err_eprox,color="black",fmt="o",label="Data")


# Now add the legend with some customizations.
legend = ax.legend(loc='upper right', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(20)

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width


#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', serif='helvetica')


ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_ylim([0,8000])
ax.set_xlabel('Energy Proxy')
ax.set_ylabel('Events')
#ax.set_title("Event Distribution over Cos (Zenith Angle)")
#plt.semilogy()
plt.show()
