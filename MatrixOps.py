import numpy as np

def Comm(A,B):
	return np.subtract(np.dot(A,B),np.dot(B,A))
def AntiComm(A,B):
	return np.add(np.dot(A,B),np.dot(B,A))
def ProjMat(i,dim):
	p = np.zeros((dim,dim))
	p[i,i]=1.0
	return p
def Trace(A):
	tr=0.0
	for i in range(0,A.shape[0]):
		tr+= A[i,i]
	return tr
