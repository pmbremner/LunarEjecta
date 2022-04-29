#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def v_speed(d, g, rs):
	return (((1./rs - np.cos(d))/(1. - np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.))**-(1/2)

def gamma_pp(d, v, rs):
	return np.arctan(1./(v**2/np.tan(d/2) + np.sqrt(v**4 / np.tan(d/2)**2 + 2.*v**2 *((1./rs - np.cos(d)) / (1. - np.cos(d))) - 1.) ))

def gamma_pm(d, v, rs):
	return np.arctan(1./(v**2/np.tan(d/2) - np.sqrt(v**4 / np.tan(d/2)**2 + 2.*v**2 *((1./rs - np.cos(d)) / (1. - np.cos(d))) - 1.) ))

def F(v, g, rs):
	return -np.sqrt(1. + (rs - 1.)/(np.cos(g)**2) * (rs + 1. - rs/v**2))

def d_dist(v, g, rs):
	Fi = F(v,g,rs)
	return np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 )

rs = float(sys.argv[1])


g = np.linspace(0., np.pi/2, 500)
v = np.linspace(0, 1., 500)
G, V = np.meshgrid(g, v)


D = d_dist(V, G, rs)

#D0 = d_dist(V, G, 1.)

plt.pcolormesh(G, V, D, cmap='jet',shading='gouraud') #gist_ncar


plt.show()