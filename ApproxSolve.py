class approx_solve(object)

#Now, use AdG's small gamma approximation in 2 flavors to extract
#oscillation parameters from the mass basis!
#M= Md
#want to map mass->flavor->decay->G->flavor->mass
#U_g^T maps flavor->decay. U_m maps mass->flavor
#our op is (U_g^H)*(U_m)

	def __init__(self,H,param):

		if param.numeu!=2:
			print "Have only implemented AdG approximation for 2 neu case"

		self.H = H
		self.param=param
		self.w, self.N = np.linalg.eig(H)

		#Calculate H matrix (IP of eigvecs)
		hm=np.zeros([param.numneu,param.numneu],complex)
		
		for i in range(0,2):
			for j in range(0,2):
				print self.N[:,i].conj(),self.N[:,j]
				hm[i,j] = np.dot(self.N[:,i].conj(),self.N[:,j])
		
		delta_mat=hm-np.identity(2,complex)
		
		print "ID"
		print np.dot(N,np.dot((np.identity(2,complex)+delta_mat.T),N.conj().T))
		print
		
		V=np.dot(N,sp.linalg.sqrtm(hm))
		
		print np.dot(V,V.conj().T)
		
		theta = math.acos(V[0,0].real)
		delta = w[1].real-w[0].real
		b1=-w[0].imag
		b2=-w[1].imag
		z= delta_mat[0,1]
		
		z_p=cmath.polar(z)
		ep=z_p[0]
		zeta=z_p[1]
		
		self.par=dict(zip(['theta','delta','b1','b2','ep','zeta'],[theta,delta,b1,b2,ep,zeta]))

	def P_ee(x):
	
		return math.exp(-2*self.par['b1']*x)*math.cos(self.par['theta'])**4 + math.exp(-2*self.par['b2']*x)*math.sin(self.par['theta'])**4 + 0.5*math.exp(-1*(self.par['b1']+self.par['b2'])*x)*math.sin(2*self.par['theta'])*(math.sin(2*self.par['theta'])*math.cos(self.par['delta']*x) - 2*self.par['ep']*math.sin(self.par['zeta'])*math.sin(self.par['delta']*x))

