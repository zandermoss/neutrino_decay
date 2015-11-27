import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import MinimalTools as MT
import SUGen as SU
import random
import cmath
import math


class HamGen(object):

	def __init__(self,param):
		self.param=param
		self.ugen=SU.SUGen(param)


	def gen(self,eig_dcy,E,matter):	
		#Randomized self.parameters to generate conjugation matrices
		#We will work in the flavor basis, so Um and Ug map from 
		#the mass basis to the flavor basis and the decay basis 
		#to the flavor basis, respectively
		#Ugen=PC.PhysicsConstants()
		
		#Generate conj matrices
		#use known mixing self.parameters for M
		Um=MT.calcU(self.param)
		
		self.ugen.sample_params()
		Ug=self.ugen.matrix_gen()	
		
		#Ugen.randomize_trig()
		#Ug=MT.calcU(Ugen)
		
		#Fill in mass and decay eigenvalues
		Md=np.zeros([self.param.numneu,self.param.numneu],complex)
		Gd=np.zeros([self.param.numneu,self.param.numneu],complex)

		#Add interaction term
		Int=np.zeros([self.param.numneu,self.param.numneu],complex)
		Int[0,0]=6.95e-11 #eV matter potential for a pure solid iron earth

	
		for i in range(0,self.param.numneu):
		
			Md[i,i]= self.param.dm2[1,i+1]/(2*E)
			Gd[i,i]= eig_dcy[i]/E #FIXME: precise form of energy dependence? 
		
	
		M= np.dot(Um,np.dot(Md,Um.conj().T))
		G= np.dot(Ug,np.dot(Gd,Ug.conj().T))

	
		#Assemble Hamiltonian
		H=np.zeros([self.param.numneu,self.param.numneu],complex)
	
		if matter:			
			H= M -1j*G + Int	
		else:
			H= M -1j*G 
		return H
	


