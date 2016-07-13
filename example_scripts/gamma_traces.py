#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import PMNSGen
import SUGen
#from tqdm import tqdm

"""
Hi Marjon!

This script is used to generate the vectors for the spaghetti plots. The basic idea is to iterate over energies, calculate
survival amplitudes at each energy and plot energy against survival amplitude. This is done using antineutrinos so we see
the MSW effect (matter resonance around ~4TeV). You will notice a loop over "t" outside the energy loop. This is iterating
over different values of phi, or the "decay mixing angle", which determines rotation in the muon-sterile neutrino subspace.
We are looking at 4 neutrinos: 3 active, 1 sterile. The PMNS matrix is as given by the "nufit" results, including best-
guess sterile-muon mixing angle. The sterile-electron mass-squared splitting is 1.0, which is close to the short-baseline
anomaly fits. All decay eigenvalues are zero except for the 4th eigenvalue. When phi is zero, this means that the only decaying
flavor vector will be the sterile. However, vacuum oscillation and matter resonance will push muon flavor into the sterile flavor,
which is a death-trap (because of that non-zero eigenvalue). That is to say, even though we start out with a pure muon flavor,
the sum of all probabilities will decay to zero over a baseline tending to infinity. As phi becomes zero, both muon and sterile 
flavors begin to decay (the eigenvalues of decay are mixtures, or superpositions, in the muon-sterile subspace). So, what we see
when we loop over t (in this case, from phi=0 to phi=3*pi/10, though you should play with this!), we are seeing the effects of 
the decay mixing as a function of energy at different phi values. It is the interplay between varying decay mixing parameters, 
vacuum oscillation, and matter resonance that makes these spaghetti plots so gnarly (and so exciting!).
Something to keep in mind: when I refer to "decay mixing", that is rotation of the gamma matrix, not rotation of the mass-
splitting matrix due to the PMNS matrix. The PMNS matrix is constant in all this! What's interesting to see is that,
with increasing phi, the (muon) survival probability can actually increase! This is unexpected, and we think it's due to some 
strange interplay between decay mixing and the MSW effect. Hopefully you can elucidate the origin of this effect!

Use gamma_plotter.py to plot these curves!
"""


"""
Loading Carlos's "mini-pdg" and setting the sterile mass-splitting to 1ev^2
"""

param=pc.PhysicsConstants(1.0)

"""
Generating a vector of energies with which to propagate. Evenly spaced on a log scale.
This goes from 100GeV to 1000 TeV.
The final argument, currently 100, is the number of samples. This is quick at 100, 
but the curves are rough. I would recommend changing it to 1000 before making pretty plots,
for smoothness.
"""

energy=np.logspace(np.log10(100),np.log10(100000),100)
#theta=np.linspace(3.141592/2,3.141592,100)

pi=3.141592

"""
1 is antineutrino, 0 neutrino. We only want antineutrinos because the matter resonance does not occur in neutrinos.
"""
ntype=1 


"""
Generating an extended PMNS matrix with the nufit results in the 3x3 submatrix defining active flavor mixing.
The 4th mixing angle, which is the sterile-muon angle, is set to a value roughly in agreement with short-baseline
fits to a sterile in 3+1.
"""

pg=PMNSGen.PMNSGen(param)
pg.sample_params()
pg.lamb[1,3]=-0.2318


"""
SUGen is initialized. There is currently no mixing, which means the eigenbasis for gamma will be the flavor basis.
Combined with the fact that the 4th eigenvalue is the only nonzero one, this statement implies that only the sterile
flavor is decaying. In the t loop, we will change the muon-sterile mixing angle "phi" to see what phenomena the concomitant
"decay mixing" induces!
"""

ug=SUGen.SUGen(param)
eig_dcy=np.zeros(4)
eig_dcy[3]=1e-14

"""
Converting a nice integer designation of particle/antiparicle into my stupid string formalism. I need to change this 
in NuSHEEP.... Mental note...
"""

if ntype==0:
	pg.param.neutype='neutrino'
	ug.param.neutype='neutrino'
	param.neutype='neutrino'
elif ntype==1:
	pg.param.neutype='antineutrino'
	ug.param.neutype='antineutrino'
	param.neutype='antineutrino'
else:
	print "BAD NEUTRINO TYPE"



"""
Begin loop over phi angles. This will iterate over [0,pi/10,2*pi/10,3*pi/10]. I know Janet wanted
you to extend this search. Manipulating this loop is a great place to start!
"""


for t in range(0,4):
	
	print "Running samples with phi=",str(t),"*pi/10"


	"""
	Set the ug mixing angle according to the desired phi.
	"""

	ug.lamb[1,3]=float(t)*pi/10.0
		
	amplitudes=np.zeros(len(energy))
	

	"""
	Begin loop over energies. You'll notice that, in the commented line, I've wrapped my iteration in something called 'tqdm'.
	This is a slick python package for displaying progress bars. If you up the resolution on the energy
	list to get smooth plots, these things will take a little while to run. It's nice to be able to keep time.
	If you want to use tqdm, uncomment the 'import tqdm' line at the top, and run 'pip install tqdm', if you have pip.
	If not, you can google 'python pip' and install it. It's just a package manager.
	This is totally non-essential, but pretty and convenient.
	"""

	#for e in tqdm(range(0,len(energy))):
	for e in range(0,len(energy)):
		amp=sp.AtmosphericNeutrinoOscillationProbability(1,1,energy[e]*param.GeV,param.PI,param,pg,ug,eig_dcy)
		amplitudes[e]=amp
	
		"""
		Write the energy list and amplitude list to a file. 4 files are written. One for each phi value.
		"""
	
		np.savez("mult_"+str(t)+"dcy_mix"+"1",e=energy,p=amplitudes)
	
