#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

p_r, p_m, R, vmax, M, p_type, s_m = np.loadtxt('crater_stats.txt', unpack=True)

# Creating bins
x_min = np.log10(np.min(vmax))
x_max = np.log10(np.max(vmax))

x_bins = np.logspace(x_min, x_max, 100)


plt.hist(vmax, bins=x_bins)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Maximum ejecta speed [m/s]')
plt.vlines(2374., 0, 1E5, color='red')
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
  
x_bins = np.logspace(x_min, x_max, 100)
y_bins = np.logspace(y_min, y_max, 100)

# https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
plt.hist2d(R, M, bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm())
plt.xlabel('Crater radius (units of impactor radius)')
plt.ylabel('Total ejecta mass (units of impactor mass)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
plt.scatter(R[p_type==0], M[p_type==0], s=1, color='red', label='MEM HI')
plt.scatter(R[p_type==1], M[p_type==1], s=1, color='green', label='MEM LO')
plt.scatter(R[p_type==2], M[p_type==2], s=1, color='blue', label='NEO')
plt.xlabel('Crater radius (units of impactor radius)')
plt.ylabel('Total ejecta mass (units of impactor mass)')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
# Creating bins
x_min = np.log10(np.min(s_m[s_m > 0]))
x_max = np.log10(np.max(s_m))

x_bins = np.logspace(x_min, x_max, 100)

plt.hist(s_m, x_bins)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('fraction of ejected mass')

plt.figure()
# Creating bins
x_min = np.log10(np.min(p_r))
x_max = np.log10(np.max(p_r))
  
y_min = np.log10(np.min(s_m[s_m > 0]))
y_max = np.log10(np.max(s_m))
  
x_bins = np.logspace(x_min, x_max, 100)
y_bins = np.logspace(y_min, y_max, 100)

# https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
plt.hist2d(p_r[p_type==0], s_m[p_type==0], bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm())
plt.xlabel('Impactor radius')
plt.ylabel('Fraction of ejected mass (high density MEM)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
# https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
plt.hist2d(p_r[p_type==1], s_m[p_type==1], bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm())
plt.xlabel('Impactor radius')
plt.ylabel('Fraction of ejected mass (low density MEM)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
# https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
plt.hist2d(p_r[p_type==2], s_m[p_type==2], bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm())
plt.xlabel('Impactor radius')
plt.ylabel('Fraction of ejected mass (NEO)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
# https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
plt.hist2d(p_r, s_m, bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm())
plt.xlabel('Impactor radius')
plt.ylabel('Fraction of ejected mass (all)')
plt.yscale('log')
plt.xscale('log')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python



plt.show()