#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

a = 12.2
b = 18.
rho_max = 1.92

def dens(d_cm):
	return rho_max * (d_cm + a) / (d_cm + b)

def dens_avg(d_cm):
	return rho_max * (1. - (a - b) * np.log(b / (d_cm + b)) / d_cm)

def dens_avg_vol_weighted(d_cm):
	return rho_max / (4. * d_cm**3) * ( d_cm*(6.*a*b - 6.*b**2 - 3.*a*d_cm + 3.*b*d_cm + 4.*d_cm**2) + 6.*(a-b)*(b**2 - d_cm**2)*np.log(b / (d_cm + b)))


print(dens_avg_vol_weighted(60))


depth_cm = np.linspace(0., 300., 1000)

plt.plot(dens(depth_cm), depth_cm, label='Hyperbolic Fit to in situ Apollo Data')
plt.plot(dens_avg(depth_cm), depth_cm, label='Running Average by Depth')
plt.plot(dens_avg_vol_weighted(depth_cm), depth_cm, label='Running Average by Crater Volume')
plt.ylabel('Depth Below Lunar Surface (cm)')
plt.xlabel(r'Bulk Density of Lunar Soil (g/cm$^3$)')
plt.xlim(1., 2.5)
plt.gca().invert_yaxis() # https://stackoverflow.com/questions/2051744/reverse-y-axis-in-pyplot
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.legend()

plt.show()
