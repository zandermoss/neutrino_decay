import pickle
import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import Splines  
import SUGen as SU
import Track
import random
import cmath
import math

import DeSolve
import HamGen
import PMNSGen
import SUGen


#Earth Density and Electron Fraction Splines
splines = Splines.Spline()

resolution=10.0**4/6371.0

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

	# prem solution
	desolve= DeSolve.DeSolve(vhamgen,myparam)
	d_amp=desolve.prop(track,initial_flavor,final_flavor)

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
