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

#Real parameters
param = PC.PhysicsConstants()

param.numneu=3

#Earth Density and Electron Fraction Splines
splines = Splines.Spline()

eig_dcy=np.zeros(param.numneu)
dcy_ord=-18.0

eig_dcy[0]=3.1989727405321533*(10.0)**dcy_ord
eig_dcy[1]=9.283607843296119*(10.0)**dcy_ord
eig_dcy[2]=6.567215979512332*(10.0)**dcy_ord

resolution=10.0**4/param.EARTHRADIUS

def AtmosphericNeutrinoOscillationProbability(initial_flavor,final_flavor,
                           energy,theta,param):

    ugen=SU.SUGen(param)
    ugen.sample_params()
    Ug=ugen.matrix_gen()

    track=Track.Track(param,resolution,energy,False)
    track.theta=theta
    track.calc_l(track.theta)

    vhamgen=HamGen.HamGen(param,Ug,track,eig_dcy,splines)

    # prem solution
    desolve= DeSolve.DeSolve(vhamgen,param)
    d_amp=desolve.prop(track,initial_flavor,final_flavor)

    return d_amp[-1]

