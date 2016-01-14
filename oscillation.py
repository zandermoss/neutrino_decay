from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np


nprobs=np.zeros(10000)
dsprobs=np.zeros(10000)
dnprobs=np.zeros(10000)
energies=np.zeros(10000)
thetas=np.zeros(10000)

nruns=10000
fname="data/run"
for run in range(0,nruns):


	f=np.load(fname+str(run)+".npz")

	energies[run]=f['energy']
	thetas[run]=f['theta']
	nprobs[run]=f['numerical'][-1]
	dsprobs[run]=f['static'][-1]
	dnprobs[run]=f['vacuum'][-1]

	#print "E:",energies[run]
	#print "Theta:",thetas[run]
	#print "Prob:",nprobs[run]
	#Progress display
	if (run+1)%((nruns)/10)==0:
		print "Done: ",run+1,"/",nruns

np.savez("scatterdata",energies=energies,thetas=thetas,nprobs=nprobs,dsprobs=dsprobs,dnprobs=dnprobs)

