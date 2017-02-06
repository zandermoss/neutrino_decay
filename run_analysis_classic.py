#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import random
from math import pi
import math
import argparse
import Verosim

parser = argparse.ArgumentParser()
parser.add_argument("nu3mass", type=float, help="nu3mass")
parser.add_argument("theta24", type=float, help="theta24")
parser.add_argument("lifetime", type=float, help="lifetime")
args = parser.parse_args()

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
#regen=True
regen=False
mtr_switch=True
#mtr_switch=False

#Choice of initial and final flavors (mu in this case)
flv_0 = 1
flv_f = 1

#Initialize PhysConst object
param=pc.PhysicsConstants()
param.numneu=4
#print "NUMNEU: ",param.numneu

#Initialize masses
nu_mass = np.zeros(param.numneu, dtype=np.float64)
nu_mass[0] = 0.0
nu_mass[1] = math.sqrt(param.dm21sq)
nu_mass[2] = math.sqrt(param.dm31sq)
nu_mass[3] = args.nu3mass #input parameter
phi_mass=0.0
param.dm2[1,4] = nu_mass[3]**2 - nu_mass[0]**2 #1ev^2


#Initialize mixing angles (CHECK BY PRINTING pg.lamb that all are correct!!)
pg=PMNSGen.PMNSGen(param)
pg.sample_params()
pg.lamb[0,3] = 0.0
pg.lamb[1,3] = -1.0*args.theta24 #input parameter
pg.lamb[2,3] = 0.0

#for i in range (0,3):
#    for j in range(0,3):
#        pg.lamb[i,j]= 0.0

#print "LAMB: "
#print pg.lamb

#Setting switches on which decay channels to use
#Decay is from second index to first, like tau!
#Mass basis, nu_4 decays
dcy_channels = np.zeros((param.numneu,param.numneu),dtype=np.uint8)
dcy_channels[0,1]=False
dcy_channels[0,2]=False
dcy_channels[0,3]=True
dcy_channels[1,2]=False
dcy_channels[1,3]=True
dcy_channels[2,3]=True
#print "DCYC: "
#print dcy_channels

#This is a group lifetime for the time being
lifetime = args.lifetime #input parameter

#Setting decay channel lifetimes (mass basis)
tau = np.zeros((param.numneu,param.numneu),dtype=np.float64)
tau[0,1]=1e+60
tau[0,2]=1e+60
tau[0,3]=lifetime
tau[1,2]=1e+60
tau[1,3]=lifetime
tau[2,3]=lifetime


ntype=1

if ntype==0:
    pg.param.neutype='neutrino'
    param.neutype='neutrino'
elif ntype==1:
    pg.param.neutype='antineutrino'
    param.neutype='antineutrino'
else:
    print "BAD NEUTRINO TYPE"

#e_vec=param.TeV*np.logspace(-1,2,5)
#print "ENERGY: ",e_vec
#ret = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,math.pi,dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)
#print ret

#Loading in energy and zenith edges for nusheep calculation
path="/home/carguelles/work/NeutrinoDecay/neutrino_decay/"
edge_file = np.load(path+"energy_zenith_edges.npz")
e_vec = edge_file['e_edges']*1.0e9
theta_vec = np.arccos(edge_file['cos_z_edges'])

#We should cut e_vec above index 150
#e_vec = edge_file['e_edges'][0:1]
#e_vec = np.logspace(11,13,100)
#theta_vec = np.arccos(edge_file['cos_z_edges'][-1:])

#print e_vec
#print theta_vec

#Initializing probability arrays
nu_prob=np.ones((len(theta_vec), len(e_vec)))
nubar_prob=np.ones((len(theta_vec), len(e_vec)))

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
        #print "starting ti: ", ti
	#nu_prob[ti,:]=ti
        nu_prob[ti,:] = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,theta_vec[ti],dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)
        #print "ti: ", ti

#print "NU PROB: "
#print nu_prob



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
        #print "starting ti: ", ti
	#nubar_prob[ti,:] = 2*ti
	nubar_prob[ti,:] = sp.AtmosphericNeutrinoOscillationProbability(flv_0,flv_f,e_vec,theta_vec[ti],dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)
	#print "ti: ", ti

#print "NUBAR PROB: "
#print nubar_prob


# TODO check that nu_prob and nubar_prob are correct. Plot them.


#Flattening probability arrays for transmission to verosimilitud
shapeprod = (nu_prob.shape[0])*(nu_prob.shape[1])
nu_prob.shape = (shapeprod)
nubar_prob.shape = (shapeprod)

#save oscillation probabilities for oscillogram
#np.save("/home/mmoulai/nu_oscillation",nu_prob)
#np.savez("/home/mmoulai/nu_oscillation", nu_prob=nu_prob, nubar_prob=nubar_prob, nu3mass=args.nu3mass, theta24=args.theta24, lifetime = args.lifetime)


#Change these as needed to point to the right stuff!
#data_path="/home/carguelles/work/NeutrinoDecay/neutrino_decay/"
#flux_path="/home/carguelles/work/NeutrinoDecay/verosimilitud/data/Marjon_Int_HondaGaisser.h5"
#effective_area_path="/home/carguelles/work/NeutrinoDecay/verosimilitud/data/effective_area.h5"
#detector_correction_path="/home/carguelles/work/NeutrinoDecay/verosimilitud/data/conventional_flux.h5"

# THIS FILES PATHS ARE FOR THE STERILE NEUTRINO MC; the new stuff
data_path="/home/carguelles/work/NeutrinoDecay/verosimilitud/data/observed_events.dat"
flux_path="/home/carguelles/work/NeutrinoDecay/verosimilitud/data/PolyGonato_QGSJET-II-04.h5"
effective_area_path="/home/carguelles/work/NeutrinoDecay/verosimilitud/data/"

# fix this paths as requiered
# note that datapath should point to the directory containtning 2010 and 2011.dat

#Running the beast
#Check years to run! -> those numbers are not 0 0 but something else. See header.

# this line is to be used with the old code
#V=Verosim.Verosim(param.numneu,0,2,data_path, flux_path, effective_area_path, detector_correction_path, nu_prob, nubar_prob)
V=Verosim.Verosim(param.numneu, data_path, flux_path, effective_area_path, nu_prob, nubar_prob)

nuis_param = np.array([1.0, 0.01, 1.0, 1.0])


# check that LLH outputs a number
#print V.LLH(nuis_param)

#Minimizing Chi2
#these bounds are +/- 3 sigma
#param_to_minimize = np.array([1, 1, 1, 1]) #1 is true, 0 is false
param_to_minimize = np.array([True, True, True, True]) #1 is true, 0 is false
low_bound = np.array([0.0001, -0.15, 0.7, 0.925])
high_bound = np.array([2.2, 0.15, 1.3, 1.075])


# TODO FIX ME
min_ret = V.MinLLH(nuis_param,low_bound,high_bound,param_to_minimize)
#print min_ret

#min_nuis = np.array([min_ret[0],min_ret[1],min_ret[2],min_ret[3]])
#print V.GetPertExpectationVec(min_nuis)

#print
print args.nu3mass, args.theta24, args.lifetime, min_ret[0], min_ret[1], min_ret[2], min_ret[3], min_ret[4]
