#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import random
from math import pi
import math

import Verosim


#Commented out command line args for the time being.
#I think these will need to be re-enabled to accept parameters
#from the batch submission script!

"""
mattervar=sys.argv[1]
regenarg = sys.argv[2]
n_calc=int(sys.argv[3])

if(regenarg=="regen"):
    regen=True
elif(regenarg=="noregen"):
    regen=False

f(mattervar=="matter"):
    mtr_switch=True
elif(mattervar=="nomatter"):
    mtr_switch=False
else:
    print "BAD INPUT: NEED TO SELECT MATTER"
"""

#Choice of matter and regeneration
regen=False
mtr_switch=True

#Choice of initial and final flavors (mu in this case)
flv_0 = 1
flv_f = 1


#Initialize PhysConst object
param=pc.PhysicsConstants()
print "NUMNEU: ",param.numneu

#Initialize masses
nu_mass = np.zeros(param.numneu, dtype=np.float64)
nu_mass[0] = 0.0
nu_mass[1] = math.sqrt(param.dm21sq)
nu_mass[2] = math.sqrt(param.dm31sq)
nu_mass[3] = 1.0
phi_mass=0.0
param.dm2[1,4] = nu_mass[3]**2 - nu_mass[0]**2

#Initialize mixing angles (CHECK BY PRINTING pg.lamb that all are correct!!)
pg=PMNSGen.PMNSGen(param)
pg.sample_params()
pg.lamb[0,3] = math.pi/4.0
pg.lamb[1,3] = math.pi/4.0
pg.lamb[2,3] = math.pi/4.0
print "LAMB: ",pg.lamb

#Setting switches on which decay channels to use
dcy_channels = np.zeros((param.numneu,param.numneu),dtype=np.uint8)
dcy_channels[0,1]=True
dcy_channels[0,2]=True
dcy_channels[0,3]=True
dcy_channels[1,2]=True
dcy_channels[1,3]=True
dcy_channels[2,3]=True

#This is a group lifetime for the time being
lifetime = 1.0e+30

#Setting decay channel lifetimes (mass basis)
tau = np.zeros((param.numneu,param.numneu),dtype=np.float64)
tau[0,1]=lifetime
tau[0,2]=lifetime
tau[0,3]=lifetime
tau[1,2]=lifetime
tau[1,3]=lifetime
tau[2,3]=lifetime

#Loading in energy and zenith edges for nusheep calculation
edge_file = np.load("energy_zenith_edges.npz")
e_vec = edge_file['e_edges']
theta_vec = np.arccos(edge_file['cos_z_edges'])

#Initializing probability arrays
nubar_prob=np.zeros((len(theta_vec), len(e_vec)))
nu_prob=np.zeros((len(theta_vec), len(e_vec)))

#Calculating neutrino probability array
ntype=0

if ntype==0:
    pg.param.neutype='neutrino'
    param.neutype='neutrino'
elif ntype==1:
    pg.param.neutype='antineutrino'
    param.neutype='antineutrino'
else:
    print "BAD NEUTRINO TYPE"

for ti in range(0,len(theta_vec)):
	nu_prob[ti,:] = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,theta_vec[ti],dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)

#Calculating antineutrino probability array
ntype=1

if ntype==0:
    pg.param.neutype='neutrino'
    param.neutype='neutrino'
elif ntype==1:
    pg.param.neutype='antineutrino'
    param.neutype='antineutrino'
else:
    print "BAD NEUTRINO TYPE"

for ti in range(0,len(theta_vec)):
	nubar_prob[ti,:] = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,theta_vec[ti],dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)


#Flattening probability arrays for transmission to verosimilitud
shapeprod = (nu_prob.shape[0])*(nu_prob.shape[1])
nu_prob.shape = (shapeprod)
nubar_prob.shape = (shapeprod)


#Change these as needed to point to the right stuff!
data_path="/home/pinkpig/physics/neutrino_decay/nusheep_change/neutrino_decay/verosimilitud/data/2010.dat"
flux_path="/home/pinkpig/physics/neutrino_decay/nusheep_change/neutrino_decay/verosimilitud/data/Marjon_Int_HondaGaisser.h5"
effective_area_path="/home/pinkpig/physics/neutrino_decay/nusheep_change/neutrino_decay/verosimilitud/data/effective_area.h5"
detector_correction_path="/home/pinkpig/physics/neutrino_decay/nusheep_change/neutrino_decay/verosimilitud/data/conventional_flux.h5"


#Running the beast
#Check years to run!
V=Verosim.Verosim(param.numneu,0,0,data_path, flux_path, effective_area_path, detector_correction_path, nu_prob, nubar_prob)

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
