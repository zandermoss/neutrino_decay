#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import SUGen
import PMNSGen
import random
import sys

ntype=int(sys.argv[1])
param=pc.PhysicsConstants(1.0)
#energy=np.logspace(10,12,10)
#energy=np.linspace(10,1000,2)
energy=np.logspace(np.log10(1000),np.log10(100000),10)
#theta=np.linspace(3.141592/2,3.141592,10)
cosz=np.linspace(-1,0,10)
probs=np.zeros((len(energy),len(cosz)))
#pg=sp.unpickle_PMNS(param)
#ntype=int(raw_input("Ntype:"))
#testnum=int(raw_input("TestNum:"))

e_edges=np.zeros(len(energy)+1)
t_edges=np.zeros(len(cosz)+1)


dt=(cosz[1]-cosz[0])
for x in range(0,len(t_edges)):
	t_edges[x]=cosz[0]-dt/2.0+x*dt

de=(np.log10(energy[1])-np.log10(energy[0]))
for x in range(0,len(e_edges)):
	e_edges[x]=10.0**(np.log10(energy[0])-0.5*de+x*de)


theta=np.arccos(cosz)

print np.log10(energy)
print np.log10(e_edges)
print theta
print t_edges

testnum=2


#outfilename="dcypyescan"+str(testnum)+"nutype_"+str(ntype)

pg=PMNSGen.PMNSGen(param)
pg.sample_params()
#print param.dm41sq
pg.lamb[1,3]=-0.2318
#print pg.lamb
#print pg.matrix_gen()


#pg=sp.unpickle_PMNS(param)
#ug=sp.unpickle_decay(param)
ug=SUGen.SUGen(param)
ug.lamb[1,3]=3.141592/10
print ug.lamb
print ug.matrix_gen()
eig_dcy=np.zeros(4)
eig_dcy[3]=1e-14


#pg=sp.unpickle_PMNS(param)
#ug=sp.unpickle_decay(param)
#ig_dcy=sp.unpickle_dcyeig(param)
#eig_dcy=np.zeros(param.numneu)

#ntype=1
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

for t in enumerate(theta):
	for e in enumerate(energy):
		print t[0]," / ",len(theta),"          ",e[0]," / ",len(energy)
		ret=sp.AtmosphericNeutrinoOscillationProbability(1,1,e[1]*param.GeV,t[1],param,pg,ug,eig_dcy)
		probs[e[0],t[0]]=ret
		#print str(e[1])+","+str( probs[e[0]])

np.savez(str(ntype)+"_1ev_2d_osc_dcy",e=e_edges,t=t_edges,p=probs)
