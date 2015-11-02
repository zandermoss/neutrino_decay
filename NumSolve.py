class num_solve(object):

	def __init__(self,H,param):
		self.H = H
		self.param=param
		self.w, self.N = np.linalg.eig(H)
		self.Ht=np.zeros([param.numneu,param.numneu],complex)
		self.Ni=np.linalg.inv(N)


	
	def prop(l):
	
		for i in range(0,self.param.numneu):
			self.Ht[i,i] = cmath.exp(-1j*self.w[i]*l)
		
		T= np.dot(self.N,np.dot(self.Ht,self.Ni))
		return T


	def scalar_prop(l,i,j):
	
		p=prop(l)
		P=p*(p.conj())
		amp=P[i,j].real
		return amp	
