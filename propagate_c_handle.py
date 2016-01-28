import pickle
import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import Splines  
import SUGen as SU
import Track
import random
import cmath
import math

import ApproxSolve
import NumSolve
import DeSolve
import HamGen
import ShellSkipSolve
import STSolve
import PMNSGen
import SUGen

#Real parameters


#Earth Density and Electron Fraction Splines
splines = Splines.Spline()
"""
eig_dcy=np.zeros(param.numneu)
dcy_ord=-18.0

eig_dcy[0]=3.1989727405321533*(10.0)**dcy_ord
eig_dcy[1]=9.283607843296119*(10.0)**dcy_ord
eig_dcy[2]=6.567215979512332*(10.0)**dcy_ord
"""
resolution=10.0**4/6371.0


def pickle_PMNS(param):

	pmnsgen=PMNSGen.PMNSGen(param)
	pmnsgen.sample_params()
	pickle.dump(pmnsgen,open("pmnsgen.p","wb"))
	return 0

def unpickle_PMNS(param):
	pmnsgen=pickle.load(open("pmnsgen.p","rb"))
	return pmnsgen

def pickle_decay(param):
	decaygen=SUGen.SUGen(param)
	decaygen.sample_params()
	pickle.dump(decaygen,open("decaygen.p","wb"))
	return 0

def unpickle_decay(param):
	decaygen=pickle.load(open("decaygen.p","rb"))
	return decaygen

def pickle_dcyeig(param,order):
	eig_dcy=np.zeros(param.numneu)
	for j in range(0,len(eig_dcy)):
		eig_dcy[j]=random.random()*order

	pickle.dump(eig_dcy,open("dcyeig.p","wb"))
	return 0

def unpickle_dcyeig(param):
	dcyeig=pickle.load(open("dcyeig.p","rb"))
	return dcyeig

def AtmosphericNeutrinoOscillationProbability(initial_flavor,final_flavor,
                           energy,theta,myparam,pmnsgen,ugen,eig_dcy):


	#ugen=SU.SUGen(myparam)
	#ugen.sample_params()
	Ug=ugen.matrix_gen()

	Um=pmnsgen.matrix_gen()
	track=Track.Track(myparam,resolution,energy,False)
	track.theta=theta
	track.calc_l(track.theta)


	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,splines)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,splines)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,None)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,None)

	# prem solution
	desolve= DeSolve.DeSolve(vhamgen,myparam)
	d_amp=desolve.prop(track,initial_flavor,final_flavor)

	#print track.l*myparam.km

	b=[]
	for x in range(0,myparam.numneu):
		b.append(np.zeros(myparam.numneu))
		b[x][x]=1.0

	amps=[]
	ptot=0
	for bv in b:
		ip=np.dot(bv,d_amp)	
		amp=np.absolute(ip)**2
		ptot+=amp
		amps.append(amp)

	print amps
	print ptot

	return amps[1],ptot
