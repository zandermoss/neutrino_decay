from scipy.interpolate import interp1d
import pickle

class Spline():

	def __init__(self):
		self.ye_spline = pickle.load( open( "ye_spline.p", "rb" ) )
		self.earth_spline = pickle.load( open( "earth_spline.p", "rb" ) )
		self.line_spline = pickle.load( open( "line_spline.p", "rb" ) )

	def GetYe(self):
		return self.ye_spline

	def GetEarth(self):
		return self.earth_spline

	def GetEarthLine(self):
		return self.line_spline
