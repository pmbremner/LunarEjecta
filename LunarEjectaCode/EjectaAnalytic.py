import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

#  python .\EjectaAnalytic.py LatRunData_Meridian > fluxFactor_Meridian.txt
#  python .\EjectaAnalytic.py LatRunData_EastLimb > fluxFactor_EastLimb.txt
#  python .\EjectaAnalytic.py LatRunData_FarSide > fluxFactor_FarSide.txt
#  python .\EjectaAnalytic.py LatRunData_WestLimb > fluxFactor_WestLimb.txt


MEO_Directory = sys.argv[1] # LatRunData_Meridian


lat = np.linspace(-90,90,37)



for i in range(0,37):

	curDirectory = "../" + MEO_Directory + "/" + "lat" + str(int(lat[i]))

	FnLoDens = curDirectory + "/LoDensity/igloo_avg.txt"
	FnHiDens = curDirectory + "/HiDensity/igloo_avg.txt"

	#print(FnLoDens, FnHiDens)

	dataLoDens = np.loadtxt(FnLoDens, unpack=True, skiprows=8)
	dataHiDens = np.loadtxt(FnHiDens, unpack=True, skiprows=8)

	Nv    = int(np.shape(dataLoDens)[0] - 9)  # number of columns (speed)
	Nrows = int(np.shape(dataLoDens)[1])      # number of rows

	speedEdge = np.linspace(0., 80., Nv+1)
	speedCenter = (speedEdge[1:] + speedEdge[0:-1])/2 # km/s

	sumLo = 0.0
	sumHi = 0.0

	for j in range(0, Nrows):
		phi = dataLoDens[7][j] * np.pi / 180. # rads

		if phi > 0: # ignore below/at the horizon
			for k in range(0, Nv):
				speed = speedCenter[k]

				sumLo += dataLoDens[9+k][j] * (speed * np.sin(phi))**(3. * 0.4)
				sumHi += dataHiDens[9+k][j] * (speed * np.sin(phi))**(3. * 0.4)

	print(int(lat[i]), sumLo, sumHi)