#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

N_env_v, N_env_zen, N_env_azm, N_env_size = np.loadtxt('ejecta_environment_flux.txt', unpack=True, max_rows=1, dtype=int)




idx, flux = np.loadtxt('ejecta_environment_flux.txt', unpack=True, skiprows=1)

# these are just index arrays
v = np.linspace(100, 5000., N_env_v)
g = np.linspace(0., np.pi, N_env_zen)
b = np.linspace(0., 2.*np.pi, N_env_azm)
a = np.logspace(-6, -1, N_env_size)

print(np.sum(flux))

#gg, vv = np.meshgrid(g, v)

#zz = gg * 0;

flux1 = flux.reshape(N_env_v, N_env_zen, N_env_azm, N_env_size)

flux_sum = np.sum(flux1, axis=2)
flux_sum = np.sum(flux_sum, axis=2)


plt.imshow(flux_sum, cmap = plt.cm.nipy_spectral, origin='lower')#, norm=mpl.colors.LogNorm())#, norm=mpl.colors.LogNorm()) # norm=mpl.colors.LogNorm()

# azimuth
flux1 = flux.reshape(N_env_v, N_env_zen, N_env_azm, N_env_size)

flux_sum = np.sum(flux1, axis=3)
flux_sum = np.sum(flux_sum, axis=1)
flux_sum = np.sum(flux_sum, axis=0)

plt.figure()
plt.semilogy(b, flux_sum)



# size
flux1 = flux.reshape(N_env_v, N_env_zen, N_env_azm, N_env_size)

flux_sum = np.sum(flux1, axis=2)
flux_sum = np.sum(flux_sum, axis=1)
flux_sum = np.sum(flux_sum, axis=0)

plt.figure()
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.loglog(a, flux_sum)

plt.show()