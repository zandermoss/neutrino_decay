import numpy as np
cimport numpy as np

from libc.math cimport fabs
from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport pi

#DTYPE = np.complex_
DTYPE = np.complex128
BDTYPE = np.uint8
RDTYPE = np.float64
ctypedef np.complex128_t DTYPE_t
#ctypedef complex_t DTYPE_t
ctypedef np.float64_t RDTYPE_t
ctypedef np.uint8_t BDTYPE_t
ctypedef np.int_t IDTYPE_t



cdef np.ndarray[DTYPE_t, ndim=2] Comm(np.ndarray[DTYPE_t, ndim=2] A,np.ndarray[DTYPE_t, ndim=2] B):
	return np.subtract(np.dot(A,B),np.dot(B,A))
cdef np.ndarray[DTYPE_t, ndim=2] AntiComm(np.ndarray[DTYPE_t, ndim=2] A,np.ndarray[DTYPE_t, ndim=2] B):
	return np.add(np.dot(A,B),np.dot(B,A))

cdef np.ndarray[DTYPE_t, ndim=2] ProjMat(int i, int dim):
	cdef np.ndarray[DTYPE_t, ndim=2] p = np.zeros((dim,dim),dtype=DTYPE)
	p[i,i]=1.0
	return p

cdef DTYPE_t Trace(np.ndarray[DTYPE_t, ndim=2] A):
	cdef DTYPE_t tr=0.0
	for i in range(0,A.shape[0]):
		tr+= A[i,i]
	return tr

"""
cdef void runstuff():
	cdef np.ndarray[DTYPE_t, ndim=2] p = np.zeros((2,2), dtype=DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=2] q = np.zeros((2,2), dtype=DTYPE)
	p[0,0]=1.0
	q[1,0]=1.0
	print Comm(p,q)

runstuff()
"""

def Regen(np.ndarray[DTYPE_t, ndim=3] rho, np.ndarray[RDTYPE_t, ndim=1] erange, int numneu, np.ndarray[BDTYPE_t, ndim=2] dcy_channels, np.ndarray[RDTYPE_t, ndim=2] pstar,  np.ndarray[RDTYPE_t, ndim=1] m_nu, RDTYPE_t m_phi, np.ndarray[RDTYPE_t, ndim=2] tau):
	cdef np.ndarray[DTYPE_t, ndim=3] R = np.zeros((erange.shape[0],numneu,numneu), dtype=DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=2] p_i = np.zeros((numneu,numneu), dtype = DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=2] p_j = np.zeros((numneu,numneu), dtype = DTYPE)
	cdef RDTYPE_t Ef 
	cdef RDTYPE_t boost
	cdef RDTYPE_t E0
	cdef int E0_arg

	for ei in range(0,len(erange)):
		Ef = erange[ei]

		for i in range(0,numneu):
			p_i = ProjMat(i,numneu)
			for j in range(i+1,numneu):
				p_j = ProjMat(j,numneu)
				if (dcy_channels[i,j]==True):
					boost = Ef/pstar[i,j]
					E0 = boost*m_nu[j]
#				   print "(I,J): (",i,",",j,")"
#				   print "BOOST: ",boost
#				   print "MJ: ", self.m_nu[j]
#				   print "MI: ", self.m_nu[i]
#				   print "MPhi: ", self.m_phi
#				   print "PSTAR: ",self.pstar[i,j]
#				   print "EF: ",Ef,"  E0: ",E0
					E0_arg = np.argmin(np.abs(erange-E0))
#				   print "EOARG: ", E0_arg
					if(E0_arg==ei):
						pass
					else:
						R[ei,:,:] += Trace(np.dot(rho[E0_arg,:,:],p_j)) * (1.0/(boost*tau[i,j])) * p_i
#				   print R[ei,:,:]
	return R



def H0_Update(DTYPE_t ye, DTYPE_t density, int numneu, np.ndarray[DTYPE_t, ndim=2] Um, BDTYPE_t anti, np.ndarray[RDTYPE_t, ndim=1] erange, np.ndarray[DTYPE_t, ndim=3] Hstart):
	cdef DTYPE_t meter = 5.06773093741e6
	cdef DTYPE_t cm = 1.0e-2*meter
	cdef DTYPE_t kg = 5.62e35 
	cdef DTYPE_t gr = 1e-3*kg 
	cdef DTYPE_t sqrt2 = 1.4142135623730951
	cdef DTYPE_t GeV = 1.0e9	
	cdef DTYPE_t GF = 1.16639e-23
	cdef DTYPE_t proton_mass = 0.938272
	
	cdef np.ndarray[DTYPE_t, ndim=2] Int = np.zeros((numneu,numneu), dtype = DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=3] ret_array = np.zeros((erange.shape[0],numneu,numneu), dtype=DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=2] MassInt 
	cdef DTYPE_t nd
	cdef DTYPE_t npotential
	cdef DTYPE_t ppotential

	nd=density*gr/(GeV*proton_mass) #convert to #proton/cm^3
	npotential=sqrt2*GF*(cm)**(-3)*nd*(1-ye) #convert cm to 1/eV
	ppotential=sqrt2*GF*(cm)**(-3)*nd*ye #convert cm to 1/eV

	for flv in range(0,3):
		#assume electron density is neutron density. This is roughly true.
		if flv==0:
			Int[flv,flv]=ppotential-0.5*npotential
		else:
			Int[flv,flv]=-0.5*npotential

	MassInt = np.dot(Um.conj().T,np.dot(Int,Um))

	if (anti==True):
		for ei in range(0,erange.shape[0]):
			ret_array[ei,:,:] = Hstart[ei,:,:]+MassInt*-1.0

	else:
		for ei in range(0,erange.shape[0]):
			ret_array[ei,:,:] = Hstart[ei,:,:]+MassInt

	return ret_array


def DrhoDt(np.ndarray[DTYPE_t, ndim=3] rho, np.ndarray[DTYPE_t, ndim=3] H0, np.ndarray[DTYPE_t, ndim=3] R, np.ndarray[DTYPE_t, ndim=3] Gamma):
	cdef np.ndarray[DTYPE_t, ndim=3] drho_dt = np.zeros((rho.shape[0],rho.shape[1],rho.shape[2]), dtype=DTYPE)
	for ei in range(0,rho.shape[0]):
		drho_dt[ei,:,:] = -1.0j*Comm(H0[ei,:,:],rho[ei,:,:]) - 0.5*AntiComm(Gamma[ei,:,:],rho[ei,:,:]) + R[ei,:,:]
	return drho_dt



def find_nearest(np.ndarray[RDTYPE_t, ndim=1] array, RDTYPE_t value):
	cdef int idx = np.searchsorted(array, value, side="left")
	if idx > 0 and (idx == array.shape[0] or fabs(value - array[idx-1]) < fabs(value - array[idx])):
		return idx-1
	else:
		return idx

def r(RDTYPE_t t, RDTYPE_t norm, RDTYPE_t theta):
	cdef RDTYPE_t x = t/norm
	return sqrt(1.0+4.0*cos(pi-theta)**2*(x**2-x))
	


