import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d as gausfilt
import dill as pickle
import matplotlib.pyplot as plt


def LineSpline():

	zero=13.08850

	dist=np.linspace(0,1.2,2)
	dense=np.linspace(zero,0,2)

	func = interp1d(dist,dense, kind='linear') #Finally, generate the spline and return.
	return func


line_spline=LineSpline()

pickle.dump(line_spline, open( "line_spline.p", "wb" ) )

"""
dist=np.linspace(0,1.2,1000)
dense=np.zeros(len(dist))

for x in range(0,len(dist)):
	dense[x]=line_spline(dist[x])

plt.plot(dist,dense,'ro')
plt.show()
"""



