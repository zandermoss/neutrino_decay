import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from numpy import linalg as LA

import PhysConst as PC
import random
import cmath
import math

## A class for generating elements of SU(N)
# A class for generating elements of SU(N) either from given angles or phases,
# or from parameters sampled (by sample_phase() or sample_angle()) from a distribution
# flat with respect to the haar measure on SU(N). For details of the construction,
# see: "Composite parameterization and Haar measure for all unitary and
# special unitary groups", arXiv:1103.3408 [math-ph]
# This class, and it's daughter: PMNSGen, are used to generate the matrices which
# diagonalize the vacuum and decay terms in the hamiltonian. 
# To perform simulations with specific (determined) values in the lambda matrix,
# one can forego the sampling functions, and set the lambda matrix elements 
# explicity by accessing the lamb member function and then simply generate the desired# matrix with matrix_gen().


class SUGen(object):

	## The constructor
	# initializes the "lambda matrix", which stores all phases and mixing angles
	# for the construction of a SU(N) matrix.
	# @par param the PhysConst object. The number of neutrinos (the N in SU(N)) is specified here. 

	def __init__(self,param):
		self.param=param

		self.lamb=np.zeros([self.param.numneu,self.param.numneu],np.complex128)

	## Samples a random lambda matrix.
	# Constructs the lambda matrix by sampling phases and angles with
	# sample_phase() and sample_angle(). The construction and scalar sampling
	# is done so that the parent distribution is flat with respect to the haar
	# measure on SU(N) (see the paper cited in this class reference for details).

	def sample_params(self):
		d=self.param.numneu

		m=1
		while (m<=d-1):
			n=0
			while (n<=m-1):
				self.lamb[m,n]=self.sample_phase()
				self.lamb[n,m]=self.sample_angle(n,m)
				n+=1
			m+=1

	## Samples a random phase.
	# @return the sampled angle.

	def sample_phase(self):	
		u=random.random()
		phi=u*math.pi*2
		return phi

	##Samples a random angle. 
	# @par m the first of two indices corresponding to the two dimensions forming the plane of rotation.
	# @par n the second such index.
	# @return the sampled angle

	def sample_angle(self,m,n):	
		u=random.random()
		k=2*(n-m)-1
		theta=math.acos((1-u)**(1.0/float(k+1)))	
		return theta

	##Generates a SU(N) matrix from the lambda matrix.
	# Proceeds as in the cited paper. If the particle is an antineutrino,
	# the matrix is conjugate-transposed.
	# Calls rotation_gen() and phase_gen() to generate phasing and rotation 
	# matrices whose product will constute the final matrix. 
	
	def matrix_gen(self):
		d=self.param.numneu
		U=np.identity(d,np.complex128)
		m=0
		while (m<=d-2):
			n=m+1
			while (n<=d-1):
				R=self.rotation_gen(m,n)
				P=self.phase_gen(m,n)
				U=np.dot(U,P)
				U=np.dot(U,R)
				n+=1
			m+=1	

		if (self.param.neutype=='antineutrino'):
			return U.conj()
		else:
			return U

	## Generates a rotation matrix.
    # @par m the first of two indices corresponding to the two dimensions forming the plane of rotation.
    # @par n the second such index.
	# @return the rotation matrix as an np array


	def rotation_gen(self,m,n):
		"""
		Here, we generate the rotation-type matrices used
		to construct the matrix element of SU(N). This is 
		simply modifying the identity with sin and cos in 
		the appropriate positions.
		"""

		R=np.identity(self.param.numneu,np.complex128)
		lcos=cmath.cos(self.lamb[m,n])
		lsin=cmath.sin(self.lamb[m,n])

		R[m,m]=lcos
		R[n,n]=lcos
		R[n,m]=(-1.0)*lsin
		R[m,n]=lsin

		return R

	## Generates a phasing matrix.
    # @par m the first of two indices which are differentially phased.
    # @par n the second such index.
	# @return the rotation matrix as an np array

	def phase_gen(self,m,n):
		"""
		Now, we generate the phasing-type matrices used in 
		the construction. We again modify the identity, 
		but this time with np.complex128 phases.
		"""
	
		P=np.identity(self.param.numneu,np.complex128)
		lpphase=cmath.exp(1.0j*self.lamb[n,m])
		lmphase=cmath.exp(-1.0j*self.lamb[n,m])

		P[m,m]=lpphase
		P[n,n]=lmphase

		return P	
