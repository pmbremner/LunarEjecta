#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

import glob

moon_x, moon_y, moon_z = np.loadtxt('moon.txt', unpack=True)

for filename in glob.glob('dists*.txt'):
	print(filename)
	if(filename != 'dists_11.txt'):
		#track_x, track_z          = np.loadtxt('tracks.txt', unpack=True, usecols=(0,1))
		track_x, track_z, track_dist, track_init_speed, track_ej_zang = np.loadtxt(filename, unpack=True)

		# plt.plot(moon_x, moon_z)
		# plt.plot(track_x, track_z, 'o', markersize=1)

		# plt.figure()
		plt.plot(track_ej_zang*180./np.pi, track_init_speed, 'o', markersize=0.5, label=filename)
	
plt.ylim([0, 1])
plt.xlim([-90, 90])
plt.legend()
plt.show()