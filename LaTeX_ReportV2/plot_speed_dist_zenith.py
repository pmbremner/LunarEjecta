#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

# https://pythonnumericalmethods.berkeley.edu/notebooks/chapter12.02-3D-Plotting.html
from mpl_toolkits import mplot3d

def vp1(g, d):
	return (1. - np.cos(2.*g) + np.sin(2.*g)/np.tan(d/2.))**(-1./2.)

def vp2(g, d, r):
	return (((1./r-np.cos(d))/(1.-np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.))**(-1./2.)


dx = np.linspace(1.E-5, 3.14, 1000)
gy = np.linspace(1.E-5, 3.14/2., 1000)

Dx, Gy = np.meshgrid(dx, gy)


V1 = vp1(Gy, Dx)

fig = plt.figure()
#ax = plt.axes(projection='3d')

plt.pcolormesh(Dx, Gy, V1)

#ax.set_zlim3d(0,1)
plt.show()