#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

N_env_v, N_env_zen, N_env_azm, N_env_size = np.loadtxt('ejecta_environment_flux.txt', unpack=True, max_rows=1, dtype=int)




flux = np.loadtxt('ejecta_environment_flux.txt', unpack=True, skiprows=1)

# these are just index arrays
v = np.arange(N_env_v)
g = np.arange(N_env_zen)
b = np.arange(N_env_azm)
a = np.arange(N_env_size)

gg, vv = np.meshgrid(g, v)

zz = gg * 0;

print(np.shape(flux))

for i in range(N_env_v):
	for j in range(N_env_zen):
		for k in range(N_env_azm):
			for l in range(N_env_size):
				idx = l + N_env_size * (k + N_env_azm * (j + N_env_zen * i))

				zz[i][j] += flux[idx]*1E10



plt.imshow(zz, cmap = plt.cm.nipy_spectral, norm=mpl.colors.LogNorm(), origin='lower')
plt.show()