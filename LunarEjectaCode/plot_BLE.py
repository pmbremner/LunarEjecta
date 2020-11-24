#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def BLE(v, d):
	d0 = 0.25 * (v / 3.)**-0.63
	if d > d0:
		if v < 2.4:
			return (d/(0.252 * (2.4 / 3.)**-0.63))**-3.6 * (v/2.4)**-1.2

	return 0.

vBLE = vectorize(BLE)



diameter = np.linspace(0.01, 1., 100)
speed    = np.linspace(0.01, 4., 100)
diameter_grid, speed_grid = np.meshgrid(diameter, speed)

a = vBLE(speed_grid.flatten(), diameter_grid.flatten()).reshape(100, 100)
fig1, ax1 = plt.subplots(constrained_layout=True)
BLE_plot = ax1.pcolormesh(speed_grid, diameter_grid, a, cmap=plt.cm.jet)
cbar = fig1.colorbar(BLE_plot)

plt.xlabel("Impact Speed (km/s)")
plt.ylabel("Critical Particle Diameter (cm)")
plt.title("Relative Penetration Risk")

plt.show()