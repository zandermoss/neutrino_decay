import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d as gausfilt
import EarthRadial as ER


class DensitySpline():
	"""
	Generates a cubic spline of radial earth density profile.
	w/r/t spline range: x is a dimensional radius. x = 0 --> center. x = 1 : Earthradius
	"""

	def __init__(self):
		self.er=ER.EarthRadial()
		self.nspl=10000 #Number of samples to take from ER
		self.splinecount=1000 #Number of points in cubic spline. Restricted by time.

	def EarthSpline(self):
		"""
		Returns a cubic spline of gaussian-smoothed density profile.
		"""

		dist=np.linspace(0,1,self.nspl,endpoint=True) #Sampling range
	
		dense=np.zeros(len(dist)) #Sample from ER
		for x in range(0,len(dist)):
			dense[x]=self.er.rdensity(dist[x])
	
	
		smooth=gausfilt(dense,10) #Smooth the sample to eliminate discontinuities. Discontinuities cause the cubic spline fitter to ring at transition points.
		smoothfunc= interp1d(dist,smooth,kind='linear') #Fit a linear spline at full resolution (very quick) so we can sample from it and fit a cubic spline.
	
	
		smalldense=np.zeros(self.splinecount) #Sampling range (reduced resolution by 10 fold)
		smalldist=np.linspace(0,1,self.splinecount,endpoint=True) 
	
		for x in range(0,len(smalldist)): #Sampling
			smalldense[x]=smoothfunc(smalldist[x])
	
		func = interp1d(smalldist, smalldense, kind='cubic') #Finally, generate the spline and return.

		return func

