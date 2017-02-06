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

import Verosim


#Commented out command line args for the time being.
#I think these will need to be re-enabled to accept parameters
#from the batch submission script!


#Loading in energy and zenith edges for nusheep calculation
edge_file = np.load("sterile_edges.npz")
e_vec = edge_file['true_energy_edges']
costh_vec = np.arccos(edge_file['costh_edges'])

nu_prob = np.ones((e_vec.size,costh_vec.size))
nubar_prob = np.ones((e_vec.size,costh_vec.size))


shapeprod = (nu_prob.shape[0])*(nu_prob.shape[1])
nu_prob.shape = (shapeprod)
nubar_prob.shape = (shapeprod)


#Change these as needed to point to the right stuff!
data_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/data/observed_events.dat"
flux_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/atmospheric_flux/averaged/PolyGonato_QGSJET-II-04.h5"
effective_area_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/systematics_response_arrays/"


#Running the beast
#Check years to run!
V=Verosim.Verosim(4,data_path, flux_path, effective_area_path, nu_prob, nubar_prob)



"""
#Minimizing Chi2
nparams = 3
param=np.zeros(nparams)
low_bound=np.zeros(nparams)
high_bound=np.zeros(nparams)
param_to_minimize=np.zeros(nparams)

min_ret = V.MinLLH(param,low_bound,high_bound,param_to_minimize)

#np.savez("sterile_oscillation_simp2",exp=expectation,exp_nopert=expectation_nopert,dat=data,chi2=retvec[0],nuisance=nuisance)
"""
