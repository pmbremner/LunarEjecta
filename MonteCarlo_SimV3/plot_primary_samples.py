#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

flux_type, azm, zenith, speed, flux, dens, mass = np.loadtxt('primary_samples.txt', unpack=True)

#print(np.sum(flux_type[flux_type == 2]/2.))
print(np.sum(flux))

fig, axs = plt.subplots(3, 2)

axs[0, 0].hist((90 - azm)%360., 50, label='azimuth')
axs[0, 0].set_title('azimuth')

axs[0, 1].hist(zenith, 50, label='zenith')
axs[0, 1].set_title('zenith')

axs[1, 0].hist(speed, 80, label='speed')
axs[1, 0].set_title('speed')

axs[1, 1].hist(flux, 80, label='flux')
axs[1, 1].set_title('flux')

axs[2, 0].hist(dens, 158, label='dens')
axs[2, 0].set_title('dens')

plt.yscale('log')
axs[2, 1].hist(np.log10(mass), 80, label='mass')
axs[2, 1].set_title('mass')


plt.figure()

plt.hist2d(speed, 90-zenith, bins=35)

plt.show()