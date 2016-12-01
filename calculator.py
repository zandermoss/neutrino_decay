#! /usr/bin/python
import matplotlib as mpl

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import SUGen
import random
import sys
import math


from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=16)


#-----------------------------------------------------------------#
#Initialize vectors and physics parameters

param=pc.PhysicsConstants()

#regen=False

nu_mass = np.zeros(param.numneu, dtype=np.float64)
nu_mass[0] = 0.0
nu_mass[1] = math.sqrt(param.dm21sq)
nu_mass[2] = math.sqrt(param.dm31sq)
nu_mass[3] = 1.0
phi_mass=0.0
param.dm2[1,4] = nu_mass[3]**2 - nu_mass[0]**2

pg=PMNSGen.PMNSGen(param)
pg.lamb[0,3] = math.pi/4.0
pg.lamb[1,3] = math.pi/4.0
pg.lamb[2,3] = math.pi/4.0

dcy_channels = np.zeros((param.numneu,param.numneu),dtype=np.uint8)
dcy_channels[0,1]=False
dcy_channels[0,2]=False
dcy_channels[0,3]=True
dcy_channels[1,2]=False
dcy_channels[1,3]=True
dcy_channels[2,3]=True

lifetime = 1.0e+2

tau = np.zeros((param.numneu,param.numneu),dtype=np.float64)
tau[0,1]=lifetime
tau[0,2]=lifetime
tau[0,3]=lifetime
tau[1,2]=lifetime
tau[1,3]=lifetime
tau[2,3]=lifetime



ntype=1 #0-neutrino, 1-antineutrino

if ntype==0:
	pg.param.neutype='neutrino'
	param.neutype='neutrino'
elif ntype==1:
	pg.param.neutype='antineutrino'
	param.neutype='antineutrino'
else:
	print "BAD NEUTRINO TYPE"


pi=3.141592

calcvar=sys.argv[1]
mattervar=sys.argv[2]
regenarg = sys.argv[3]
n_calc=int(sys.argv[4])



if (calcvar=="energy"):
	valvec=np.logspace(-1,4,n_calc)

elif (calcvar=="theta"):
	valvec=np.linspace(0,2*pi,n_calc)
else:
	print "BAD INPUT: NEED TO SELECT MODE"

if(regenarg=="regen"):
	regen=True
elif(regenarg=="noregen"):
	regen=False


if(mattervar=="matter"):
	mtr_switch=True
	potential=2.65492195354e-13
elif(mattervar=="nomatter"):
	mtr_switch=False
	potential=0
else:
	print "BAD INPUT: NEED TO SELECT MATTER"


frac_vec=np.zeros(len(valvec))
num_vec=np.zeros(len(valvec))

#-----------------------------------------------------------------#
#Calculate Probabilities
#
#for val in range(0,len(valvec)):
#
#	if (calcvar=="energy"):
#		E=valvec[val]*param.TeV
#	else:
#		E=4.0*param.TeV
#	
#	if (calcvar=="theta"):
#		theta=valvec[val]
#	else:
#		theta=0.2318



	#diag_amp = sp.DiagProb(0,0,E,param.PI,param,pg,ug,eig_dcy,dcy_switch,mtr_switch,2*param.EARTHRADIUS*param.km)


E = valvec*param.TeV
print "E: ",E
amps= sp.AtmosphericNeutrinoOscillationProbability(3,3,E,param.PI,dcy_channels,tau,param,pg,nu_mass,phi_mass,regen,mtr_switch)

num_vec = amps

print "NV: ",num_vec
#print sumamps

np.savez(calcvar+"_"+regenarg+"_von_neumann_spline_"+mattervar,_E=E,_potential=potential,_valvec=valvec,_num_vec=num_vec,_ntype=ntype)

