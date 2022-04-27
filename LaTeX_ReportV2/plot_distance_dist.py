#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def v_speed(d, g, rs):
	return (((1./rs - np.cos(d))/(1. - np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.))**-(1/2)

def dM(d, g, h0, h1, a):
	mu = 0.4
	v0 = v_speed(d, g, h0) 
	v1 = v_speed(d, g, h1)
	v2 = v_speed(d + 2.*a, g, h1)
	v3 = v_speed(d + 2.*a, g, h0)

	varr = np.array([v0, v1, v2, v3, 1])

	vmin = np.min(varr)
	vmax = np.nanmax(varr)

	if np.isnan(vmin) or np.isnan(vmax):
		return 0.
	else:
		return vmin**(-3.*mu) - vmax**(-3.*mu)

dMvec = vectorize(dM)


def M(d, h0, h1, a):
	g = np.linspace(0., np.pi/2., 100)

	return np.sum(dMvec(d, g, h0, h1, a)) * (g[1] - g[0])	

Mvec = vectorize(M)


def Mtot(h0, h1, a):
	d = np.logspace(-4, np.log10(np.pi), 100)

	return np.sum(np.multiply(Mvec(d[:-1], h0, h1, a), d[1:] - d[:-1]))

Mtotvec = vectorize(Mtot)


a = np.logspace(-3, -1, 8)
h1 = np.logspace(-3, -1, 8)

# a = np.array([0.01, 0.02, 0.04])

# h1 = np.array([1.01, 1.02, 1.04])


# Mtot_array = np.zeros(h1.size)

# for i in range(h1.size):
# 	Mtot_array[i] = Mtot(1., h1[i], a[0])

for h1i in h1:
	plt.plot(a, Mtotvec(1., h1i, a))

plt.ylabel('Fraction of ejecta')
plt.xlabel('Asset radius (rm)')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.show()