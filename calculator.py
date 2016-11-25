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


from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=16)


#-----------------------------------------------------------------#
#Initialize vectors and physics parameters

param=pc.PhysicsConstants(1.0)


nu_mass = np.zeros(len(param.numneu))
nu_mass[3]=1.0
phi_mass=0.0
param.dm2[1,2] = nu_mass[1]**2 - nu_mass[0]**2
param.dm2[1,3] = nu_mass[2]**2 - nu_mass[0]**2
param.dm2[2,3] = nu_mass[2]**2 - nu_mass[1]**2
param.dm2[1,4] = nu_mass[3]**2 - nu_mass[0]**2
param.dm2[2,4] = nu_mass[3]**2 - nu_mass[1]**2
param.dm2[3,4] = nu_mass[3]**2 - nu_mass[2]**2

pg=PMNSGen.PMNSGen(param)

tau = np.zeros((param.numneu,param.numneu))
tau[0,1]=1.0e-16
tau[0,2]=1.0e-16
tau[0,3]=1.0e-16
tau[1,2]=1.0e-16
tau[1,3]=1.0e-16
tau[2,3]=1.0e-16



ntype=1 #0-neutrino, 1-antineutrino

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


pi=3.141592

calcvar=sys.argv[1]
mattervar=sys.argv[2]
n_calc=int(sys.argv[3])



if (calcvar=="energy"):
	valvec=np.logspace(0,6,n_calc)

elif (calcvar=="theta"):
	valvec=np.linspace(0,2*pi,n_calc)
else:
	print "BAD INPUT: NEED TO SELECT MODE"



dcy_switch=True

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



num_vec = sp.AtmosphericNeutrinoOscillationProbability(0,0,E,param.PI,tau,param,pg,nu_mass,phi_mass)


np.savez(calcvar+"_"+mattervar,_theta=theta,_E=E,_potential=potential,_valvec=valvec,_num_amp=num_amp,_ntype=ntype)

