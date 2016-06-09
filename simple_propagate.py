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


#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,splines)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,splines)
#	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,splines)
	vhamgen=HamGen.HamGen(myparam,Um,Ug,track,eig_dcy,None)
	#vhamgen=HamGen.HamGen(myparam,Um,Ug,track,None,None)

	# prem solution
	desolve= DeSolve.DeSolve(vhamgen,myparam)
	#stsolve= STSolve.STSolve(vhamgen,myparam)
	#sssolve= ShellSkipSolve.ShellSkipSolve(vhamgen,myparam)
	#d_amp=sssolve.prop(track,initial_flavor,final_flavor)
	dvec,pvec=desolve.prop(track,initial_flavor,final_flavor)
	#st_amp=stsolve.prop(track,initial_flavor,final_flavor)
	#print dvec,pvec

	##return [0,0]
	#print track.l*myparam.km
	b=[]
	for x in range(0,myparam.numneu):
		b.append(np.zeros(myparam.numneu))
		b[x][x]=1.0

	amps=np.zeros([len(dvec),4])
	for x in range(0,len(dvec)):
		for bv in enumerate(b):
			ip=np.dot(bv[1],pvec[x])	
			amps[x,bv[0]]=np.absolute(ip)**2


	return dvec,amps
def gridrun():

	param=PC.PhysicsConstants()
	print "SHEEP"
	print param.GeV
	print param.TeV
	theta = np.linspace(param.PI/2.0, param.PI, 10)
	energy = np.logspace(9, 12, 10)
	print theta
	print energy
	prob=np.zeros((10,10))
	for e in enumerate(energy):
		for t in enumerate(theta):
			print e,t
			print "OWL"
			prob[e[0],t[0]]=AtmosphericNeutrinoOscillationProbability(1,1,e[1],t[1],param)
			

	print prob
	E, T = np.meshgrid(energy, theta)


	plt.figure()
	CS = plt.contourf(E,T,prob,cmap=plt.cm.jet)
	plt.clabel(CS, inline=1, fontsize=10)
	cbar = plt.colorbar(CS)
	plt.show()	



#gridrun()


