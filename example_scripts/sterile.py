#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import random

import Verosim


V=Verosim.Verosim(4,"2010")
param=pc.PhysicsConstants(4)
ntype=1
testnum=1



pg=PMNSGen.PMNSGen(param)
pg.sample_params()
pg.lamb[1,3]=-0.2318238

print pg.lamb

#pg=sp.unpickle_PMNS(param)
	
ug=sp.unpickle_decay(param)
eig_dcy=sp.unpickle_dcyeig(param)



cuts=np.zeros(2)
cuts[0]=6
cuts[1]=50-16

print "SETTING CUTS"
V.SetEproxCuts(cuts)
print "SET CUTS"
V.SetSimpsNIntervals(2)

def wrapcall(arg):
	nparg=np.asarray(arg)

	if nparg[2]<0.4:
		pg.param.neutype='neutrino'
		ug.param.neutype='neutrino'
		param.neutype='neutrino'
		print "NEU"

	elif nparg[2]>0.6:
		pg.param.neutype='antineutrino'
		ug.param.neutype='antineutrino'
		param.neutype='antineutrino'
		print "ANTINEU"

	olv=sp.AtmosphericNeutrinoOscillationProbability(1,1,nparg[0]*param.GeV,nparg[1],param,pg,ug,eig_dcy)
	return olv[0]
V.SetDeSolver(wrapcall)

"""
ns=[]
vals=[]
diffs=[]
last=0
for x in range(1,5):
	val=V.SimpsAvg(1.57,3.14,5,10,0,2*x)

	ns.append(2*x)
	vals.append(val)
	diffs.append(abs(val-last))

#	print "N: ",2*x," AVG: ",val
#	print "N: ",2*x," DIFF: ",abs(val-last)
	last=val

for x in zip(ns,vals,diffs):
	print x
"""


V.CalculateExpectation()

init=np.zeros(2)
init[0]=1
init[1]=0

retvec=V.Chi2MinNuisance(init)
print retvec

nuisance=np.zeros(2)
nuisance[0]=retvec[1]
nuisance[1]=retvec[2]

eprox=np.asarray( V.GetEproxEdges())
cosz=np.asarray(V.GetCosZenithEdges())

raw_expectation= V.GetPertExpectationVec(nuisance)
raw_data= V.GetDataVec()
raw_expectation_nopert= V.GetExpectationVec()


baddims=V.GetExpDims()
dims=[0,0]
dims[0]=baddims[0]
dims[1]=baddims[1]
print "DIMS: ", dims

#deserialze
datasize=1
strides=[]
for i in range(0,2):
    strides.append(datasize)
    datasize*=dims[2-(i+1)]

index=0

expectation = np.zeros(dims)
expectation_nopert = np.zeros(dims)
data = np.zeros(dims)

def index(indices,sheep):
    myindex=0
    for i in range(0,2):
        myindex+=strides[i]*indices[2-(i+1)]
    #print raw_expectation[myindex] 
    return sheep[myindex]


for x in range(0,dims[0]):
    for y in range(0,dims[1]):
        myind=[x,y]
        expectation[x,y]=index(myind,raw_expectation)
        expectation_nopert[x,y]=index(myind,raw_expectation_nopert)
        data[x,y]=index(myind,raw_data)

print data 

np.savez("sterile_oscillation_simp2",exp=expectation,exp_nopert=expectation_nopert,dat=data,chi2=retvec[0],nuisance=nuisance)
