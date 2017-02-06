#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import random
from math import pi
import math
from tqdm import tqdm
import emcee
from getdist import plots, MCSamples
import getdist

import Verosim

#Set paths to data, fluxes, and responses.
data_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/data/observed_events.dat"
flux_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/atmospheric_flux/averaged/PolyGonato_QGSJET-II-04.h5"
effective_area_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/systematics_response_arrays/"


#Load oscillation probability arrays
osc_file = np.load("3flav_matter_osc.npz")
nu_prob = osc_file['nu']
nubar_prob = osc_file['nubar']
#nu_prob = np.ones(len(nu_prob))
#nubar_prob = np.ones(len(nubar_prob))

#Calling Verosimilitud
V=Verosim.Verosim(4,data_path, flux_path, effective_area_path, nu_prob, nubar_prob)

#Testing Chi2
"""
  There are 5 nuisance parameters at the moment.
  #  Name                        Mean    Sigma
  --------------------------------------------
  1) Normalization coefficient   1.0     0.4
  2) Spectral index              0.0     0.05
  3) Kaon-pion ratio             1.0     0.1
  4) Nu-nubar ratio              1.0     0.025
  5) DOM efficiency              1.0     0.3
"""

def llh(params):
	return V.GetLLH(params)

ndim,nwalkers = 5,100
ivar = np.array([1.0,0.0,1.0,1.0,1.0])
p0 = [ivar + 1e-2*np.random.rand(ndim) for i in range(nwalkers)]

prior_cov = np.diag((0.4,0.05,0.1,0.025,0.3))	
prior_cov = prior_cov**2
print prior_cov
prior_means = np.array([1.0,0.0,1.0,1.0,1.0])
prior_vecs = np.random.multivariate_normal(prior_means,prior_cov,100000)


sampler = emcee.EnsembleSampler(nwalkers,ndim,llh)
pos,prob,state = sampler.run_mcmc(p0,1000)
sampler.reset()
sampler.run_mcmc(pos,1000)



names = ["N","\gamma","R_{k\pi}","R_{\\nu \\bar\\nu}","\epsilon"]
labels = ["N","\gamma","R_{k\pi}","R_{\\nu \\bar\\nu}","\epsilon"]

samples = MCSamples(samples=sampler.flatchain, names=names, labels=labels)
prior_samples = MCSamples(samples=prior_vecs, names=names, labels=labels)
g = plots.getSubplotPlotter()
#g.triangle_plot(samples)
g.triangle_plot([prior_samples,samples],filled=True, legend_labels = ['Prior','Posterior'],title="Gaussian Initialization, 100 Step Burn-in.")
#g.triangle_plot(samples,filled=True)
g.export('mcmc.pdf')



"""
#Minimizing Chi2
nparams = 5
low_bound=np.array([0.0,-1.0,0.0,0.0,0.91])
high_bound=np.array([2.0,1.0,2.0,2.0,1.1978])
param_to_minimize=np.ones(nparams)
min_ret = V.MinLLH(param,low_bound,high_bound,param_to_minimize)
print min_ret
"""
