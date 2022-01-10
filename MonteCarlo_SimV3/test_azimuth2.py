import sys
import numpy as np
import matplotlib.pyplot as plt


def Di(azmi, D, R, pm):
	# a = 1. - (np.sin(azmi) * np.sin(D+R))**2
	# b = -2. * np.cos(D+R) * np.cos(R)
	# c = (np.cos(D+R) * np.cos(R))**2 - (np.sin(D+R) * np.sin(R))**2 + (np.sin(D+R) * np.sin(azmi))**2

	a = (np.cos(R))**2 - (np.sin(D+R) * np.cos(azmi))**2
	b = -np.sin(2.*(D+R)) * np.cos(azmi)
	c = (np.cos(R))**2 - (np.cos(D+R))**2

	return np.arctan2((-b + pm*np.sqrt(b**2 - 4.*a*c)), (2.*a))

D = float(sys.argv[1])
R = float(sys.argv[2])

azm_max = np.arcsin((np.sin(R) / np.sin(D+R)))

azm = np.linspace(0., azm_max, 1000)

d0 = Di(azm, D, R, 1.)
d1 = Di(azm, D, R, -1.)


plt.plot(azm/np.pi*180., d0, label="d0")
plt.plot(azm/np.pi*180., d1, label="d1")
plt.legend()
plt.show()

