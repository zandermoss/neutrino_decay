import numpy as np

def Comm(A,B):
	return (np.dot(A,B)-np.dot(B,A))
def AntiComm(A,B):
	return (np.dot(A,B)+np.dot(B,A))
def ProjMat(i,dim):
	p = np.zeros((dim,dim))
	p[i,i]=1.0
	return ProjMat
