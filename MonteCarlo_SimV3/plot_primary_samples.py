#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

flux_type, azm, zenith, speed, flux = np.loadtxt('primary_samples.txt', unpack=True)

fig, axs = plt.subplots(2, 2)

axs[0, 0].hist((90 - azm)%360., 50, label='azimuth')
axs[0, 0].set_title('azimuth')
axs[0, 1].hist(zenith, 50, label='zenith')
axs[0, 1].set_title('zenith')
axs[1, 0].hist(speed, 80, label='speed')
axs[1, 0].set_title('speed')


plt.figure()

plt.hist2d(speed, 90-zenith, bins=35)

plt.show()