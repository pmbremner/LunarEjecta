import sys
import numpy as np
import matplotlib.pyplot as plt

def Ci(azmi, D, R):
	return 2.*np.sqrt((np.sin(R))**2 - (np.sin(D+R) * np.sin(azmi))**2)

def cosai(azmi, D, R):
	return -np.sqrt(1. - ((np.sin(D+R) * np.sin(azmi)) / (np.sin(R)))**2)

def Di(azmi, D, R):
	#c = Ci(azmi, D, R)

	#return np.arctan2(-2.*((np.sin(D+R))**2 - R**2), np.sin(2.*(D+R))*np.cos(azmi) - c)
	return np.arctan2((np.cos(R))**2 - (np.cos(D+R))**2, np.cos(D+R)*np.sin(D+R)*np.cos(azmi) - np.cos(R)*np.sin(R)*cosai(azmi, D, R))
	

D = float(sys.argv[1])
R = float(sys.argv[2])

azm_max = np.arcsin(np.sin(R) / np.sin(D+R))

azm = np.linspace(0., azm_max, 1000)

d0 = Di(azm, D, R)
d1 = d0 + Ci(azm, D, R)


plt.plot(azm/np.pi*180., d0, label="d0")
plt.plot(azm/np.pi*180., d1, label="d1")
plt.legend()
plt.show()

