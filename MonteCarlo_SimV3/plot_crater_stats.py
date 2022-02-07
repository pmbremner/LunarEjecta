#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

p_r, p_m, R, vmax, M, p_type = np.loadtxt('crater_stats.txt', unpack=True)

# Creating bins
x_min = np.log10(np.min(vmax))
x_max = np.log10(np.max(vmax))

x_bins = np.logspace(x_min, x_max, 100)


plt.hist(vmax, bins=x_bins)
plt.yscale('log')
plt.xscale('log')
plt.figure()

p_r = p_r[vmax > 0]
p_m = p_m[vmax > 0]
R = R[vmax > 0]
M = M[vmax > 0]
p_type = p_type[vmax > 0]

# Creating bins
x_min = np.log10(np.min(R))
x_max = np.log10(np.max(R))
  
y_min = np.log10(np.min(M))
y_max = np.log10(np.max(M))
  
x_bins = np.logspace(x_min, x_max, 1000)
y_bins = np.logspace(y_min, y_max, 1000)

# https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
plt.hist2d(R, M, bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm())
plt.xlabel('Crater radius (units of impactor radius)')
plt.ylabel('Total ejecta mass (units of impactor mass)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
plt.scatter(R[p_type==0], M[p_type==0], s=1, color='red')
plt.scatter(R[p_type==1], M[p_type==1], s=1, color='green')
plt.scatter(R[p_type==2], M[p_type==2], s=1, color='blue')
plt.xlabel('Crater radius (units of impactor radius)')
plt.ylabel('Total ejecta mass (units of impactor mass)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python


plt.show()