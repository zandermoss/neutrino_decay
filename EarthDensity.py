import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d as gausfilt
import EarthRadial as ER


class EarthDensity():
    """
    Generates a cubic spline of radial earth density profile.
    w/r/t spline range: x is a dimensional radius. x = 0 --> center. x = 1 : Earthradius
    """

    self.EarthRadius = 6371.0    #[km]
    self.Radius = 6371.0    #[km]

    def __init__(self):
        self.nspl=10000 #Number of samples to take from ER
        self.splinecount=1000 #Number of points in cubic spline. Restricted by time.

    def rdensity(self,x):
        # Calcula la densidad de la Tierra segun el PREM
        # R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        # Arxiv : 9512364 pag. 23
        # x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        r = self.EarthRadius*x
        if r <= 1221.50 :
            dne = 13.08850-8.83810*x**2
        elif r>=1221.50 and r<3480 :
            dne=12.58150-1.26380*x-3.64260*x**2.-5.5280*x**3.
        elif r >=3480.0 and r < 5701.0 :
            dne=7.95650-6.47610*x+5.52830*x**2.-3.08070*x**3.
        elif r >= 5701.0 and r<5771.0 :
            dne=5.31970-1.48360*x
        elif r>=5771.0 and r<5971.0 :
            dne=11.24940-8.02980*x
        elif r>=5971.0 and r<6151.0 :
            dne=7.10890-3.80450*x
        elif r>=6151.0 and r<6346.60 :
            dne=2.6910+0.69240*x
        elif r >= 6346.60 and r < 6356.0 :
            dne = 2.9
        elif r >= 6356.0 and r < 6368 :
            dne = 2.6
        elif r<= self.EarthRadius :
            dne = 1.020
        elif r>=self.EarthRadius :
            dne=0.0
        return dne


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


