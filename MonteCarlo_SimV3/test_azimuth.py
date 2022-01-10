import sys
import numpy as np
import matplotlib.pyplot as plt


def cosai(azmi, D, R):
	return np.sin(azmi) * np.cos(D + R)

def ri(azmi, D, R):
	return np.arctan2(np.tan(azmi) * np.sin(D + R), 1.)

def c0i(azmi, D, R):
	return ri(azmi, D, R) * cosai(azmi, D, R) + np.sqrt(R**2 - (1. - cosai(azmi, D, R)**2) * ri(azmi, D, R)**2)

def c1i(azmi, D, R):
	return -ri(azmi, D, R) * cosai(azmi, D, R) + np.sqrt(R**2 - (1. - cosai(azmi, D, R)**2) * ri(azmi, D, R)**2)



D = float(sys.argv[1])
R = float(sys.argv[2])

#azm_max = np.arctan2(np.tan(R), np.sin(D + R))
azm_max = np.arccos((np.cos(R) - np.cos(D + R)**2) / (np.sin(D + R)**2))

azm = np.linspace(0., azm_max, 100)

c0 = c0i(azm, D, R)
c1 = c1i(azm, D, R)

print(np.sum((c0+c1)/(2.*R) * azm[1] / azm_max))

plt.plot(azm/np.pi*180., (c0+c1)/(2.*R), label="c0")
#plt.plot(azm, D + R - c1, label="c1")
plt.show()

