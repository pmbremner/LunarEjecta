#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os


import glob

for file in glob.glob("*speed_vs_angle*.txt"):

	if os.stat(file).st_size > 0:
		sample_latp, sample_lonp, sample_azimuth_0, sample_zenith_0, sample_speed_0, sample_weight, sample_azimuth_f, sample_zenith_f, sample_speed_f = np.loadtxt(file, unpack=True)
		
		#idx = sample_speed_0 < 0.8

		#print(sample_zenith_f[idx], sample_speed_f[idx])

		plt.xlim(0., np.pi)
		plt.ylim(0., 5.)
		plt.scatter(sample_zenith_0, sample_speed_0, s=1, color='red')
		plt.scatter(sample_zenith_f, sample_speed_f, s=1, color='blue')
		plt.xlabel("Zenith Angle (rad)")
		plt.ylabel("Ejecta Speed (v_esc)")
		plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
		plt.savefig(file[:-4] + ".png")
		plt.show()
		plt.clf() # https://stackoverflow.com/questions/8213522/when-to-use-cla-clf-or-close-for-clearing-a-plot-in-matplotlib

		plt.figure()
		plt.scatter(sample_zenith_f, sample_weight, s=1, color='red')
		plt.ylim(np.min(sample_weight[sample_weight>0]), np.max(sample_weight))
		plt.yscale('log')
		plt.show()
		plt.clf()

		plt.figure()
		plt.scatter(sample_speed_f, sample_weight, s=1, color='blue')
		plt.ylim(np.min(sample_weight[sample_weight>0]), np.max(sample_weight))
		plt.yscale('log')
		plt.show()
		plt.clf()

