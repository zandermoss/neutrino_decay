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

regen=False
mtr_switch=True

#Choice of initial and final flavors (mu in this case)
flv_0 = 1
flv_f = 1

#Initialize PhysConst object
param=pc.PhysicsConstants()
param.numneu=4
print "NUMNEU: ",param.numneu

#Initialize masses
nu_mass = np.zeros(param.numneu, dtype=np.float64)
nu_mass[0] = 0.0
nu_mass[1] = math.sqrt(param.dm21sq)
nu_mass[2] = math.sqrt(param.dm31sq)
nu_mass[3] = 1.0
phi_mass=0.0
param.dm2[1,4] = nu_mass[3]**2 - nu_mass[0]**2 #1ev^2


#Initialize mixing angles (CHECK BY PRINTING pg.lamb that all are correct!!)
pg=PMNSGen.PMNSGen(param)
pg.sample_params()
pg.lamb[0,3] = 0.0
pg.lamb[1,3] = 0.0
pg.lamb[2,3] = 0.0
print "LAMB: ",pg.lamb

#Setting switches on which decay channels to use
#Decay is from second index to first, like tau!
dcy_channels = np.zeros((param.numneu,param.numneu),dtype=np.uint8)
dcy_channels[0,1]=False
dcy_channels[0,2]=False
dcy_channels[0,3]=False
dcy_channels[1,2]=False
dcy_channels[1,3]=False
dcy_channels[2,3]=False
print "DCYC: ",dcy_channels

#This is a group lifetime for the time being
lifetime = 1.0e+4

#Setting decay channel lifetimes (mass basis)
tau = np.zeros((param.numneu,param.numneu),dtype=np.float64)
tau[0,1]=1e+60
tau[0,2]=1e+60
tau[0,3]=lifetime
tau[1,2]=1e+60
tau[1,3]=lifetime
tau[2,3]=lifetime


#Loading in energy and zenith edges for nusheep calculation
edge_file = np.load("sterile_edges.npz")
e_vec = edge_file['true_energy_edges']
e_vec *= param.GeV
cosz_vec = edge_file['costh_edges']
theta_vec = np.arccos(cosz_vec)

#Initializing probability arrays
nubar_prob=np.zeros((len(theta_vec), len(e_vec)))
nu_prob=np.zeros((len(theta_vec), len(e_vec)))

#Calculating neutrino probability array
pg.param.neutype='neutrino'
param.neutype='neutrino'

print "Calculating Neutrino Survival Probabilities"

for ti in tqdm(range(0,len(theta_vec))):
	nu_prob[ti,:] = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,theta_vec[ti],dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)

#Calculating antineutrino probability array
pg.param.neutype='antineutrino'
param.neutype='antineutrino'

print "Calculating Antineutrino Survival Probabilities"

for ti in tqdm(range(0,len(theta_vec))):
	nubar_prob[ti,:] = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,theta_vec[ti],dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)

#Flattening probability arrays for transmission to verosimilitud
shapeprod = (nu_prob.shape[0])*(nu_prob.shape[1])
nu_prob.shape = (shapeprod)
nubar_prob.shape = (shapeprod)

np.savez("3flav_matter_osc",ev=e_vec,tv=theta_vec,nu=nu_prob,nubar=nubar_prob)

