import numpy as np
import matplotlib.pyplot as plt

f=np.load("scatterdata.npz")

fig = plt.figure()

zenith=f['thetas']/2.0
e_gev=f['energies']*1000.0


cutindex=np.where(e_gev>500)
#cutindex=np.where(e_gev>-1)


#plt.scatter(zenith,e_gev,edgecolors='none',c=np.log(f['nprobs']))

fig=plt.figure()
ax=fig.add_subplot(111,axisbg='black')
mappable=ax.scatter(zenith[cutindex],e_gev[cutindex],edgecolors='none',c=f['dnprobs'][cutindex])
cbar=plt.colorbar(mappable,pad=0.05)
cbar.set_label('Nu Mu Survival Probability')
ax.set_title('Oscillation Probabilities (mu->mu) (Numerical Dynamic Matter)')

ax.set_xlabel("Zenith Angle (Radians)")
ax.set_ylabel("Energy (GeV)")

plt.show()
#fig.canvas.draw()

