#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def v_speed(d, g, rs):
	return (((1./rs - np.cos(d))/(1. - np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.))**-(1/2)
	#return ((((1. - rs)/rs + 2.*np.sin(d/2.))/(np.sin(d/2.)))*(np.sin(g)) + np.sin(2.*g)/np.tan(d/2.))**-(1/2)

def gamma_pp(d, v, rs):
	return np.arctan(1./(v**2/np.tan(d/2) + np.sqrt(v**4 / np.tan(d/2)**2 + 2.*v**2 *((1./rs - np.cos(d)) / (1. - np.cos(d))) - 1.) ))

def gamma_pm(d, v, rs):
	return np.arctan(1./(v**2/np.tan(d/2) - np.sqrt(v**4 / np.tan(d/2)**2 + 2.*v**2 *((1./rs - np.cos(d)) / (1. - np.cos(d))) - 1.) ))

def F(v, g, rs, d0):
	return -np.sqrt(1. + (rs - 1.)/(np.cos(g)**2) * (rs + 1. - rs/v**2))

def g_ap(d, rs):
	return np.arctan2(rs*np.tan(d/2.), rs - 1.)

def d_dist(v, g, rs, d0):
	Fi = F(v, g, rs, d0)
	if g > g_ap(d0, rs):
		Fi *= -1
		if v**2 > 1./(1. + 1./rs):
			return np.mod(np.pi + np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) , 2.*np.pi)
		else:
			return np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) 
	return np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) 

d_distvec = vectorize(d_dist)

def v_min(g, rs):
	return np.sqrt((2.*rs*(rs-1.)) / (np.cos(2.*g)-1. + 2.*rs**2))

rs = float(sys.argv[1])
dd  = float(sys.argv[2])

g = np.linspace(0.0001, np.pi/2*0.9999, 300)
v = np.linspace(0.0001, 1., 300)
G, V = np.meshgrid(g, v)


D = d_distvec(V, G, rs, dd)

#D0 = d_dist(V, G, 1.)

plt.pcolormesh(G, V, D, cmap='prism')#,shading='gouraud') #gist_ncar, nipy_spectral, gist_rainbow, plasma

plt.plot(g, v_min(g, rs), color='black', linewidth=4, linestyle='--')

plt.plot(g, v_speed(dd, g, rs), color='black', linewidth=3)
#plt.plot(g, v_speed(dd+0.01, g, rs), color='black', linewidth=3)
#plt.plot(g, v_speed(dd, g, 1.), color='black', linewidth=3)

plt.ylim(0,1)
plt.xlabel('Initial Zenith Angle (rad)')
plt.ylabel('Initial Speed (Vesc)')

plt.show()