import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import random
import cmath
import math


class HamGen(object):

	def __init__(self,param):
		self.param=param


	def gen(self,eig_dcy):	
		#Randomized self.parameters to generate conjugation matrices
		#We will work in the flavor basis, so Um and Ug map from 
		#the mass basis to the flavor basis and the decay basis 
		#to the flavor basis, respectively
		Ugen=PC.PhysicsConstants()
		
		#Generate conj matrices
		#use known mixing self.parameters for M
		Um=MT.calcU(self.param)
		
		Ugen.randomize_trig()
		Ug=MT.calcU(Ugen)
		
		#Fill in mass and decay eigenvalues
		Md=np.zeros([self.param.numneu,self.param.numneu],complex)
		Gd=np.zeros([self.param.numneu,self.param.numneu],complex)
		
		for i in range(0,self.param.numneu+0):
		
		    print "m2:", self.param.dm2[1,i+1] #FIXME: mass constants-> Add in energy dependence!
		    Md[i,i]= self.param.dm2[1,i+1] #FIXME: mass constants-> Add in energy dependence!
		    #Md[i,i]= 1 #FIXME: mass constants-> Add in energy dependence!
		    Gd[i,i]= eig_dcy[i]*1e-5 #FIXME
		
		
		print Md
		
		M= np.dot(Um,np.dot(Md,Um.conj().T))
		G= np.dot(Ug,np.dot(Gd,Ug.conj().T))
		
		#Assemble Hamiltonian
		H=np.zeros([self.param.numneu,self.param.numneu],complex)
		
		H= M -1j*G
		
		print H
	
		return H
	


