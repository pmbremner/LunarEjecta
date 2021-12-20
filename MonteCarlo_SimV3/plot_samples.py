#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os


import glob

for file in glob.glob("*speed_vs_angle*.txt"):

	if os.stat(file).st_size > 0:
		zenith, speed, weight = np.loadtxt(file, unpack=True)
		plt.xlim(0., np.pi/2.)
		plt.ylim(0., 5.)
		plt.scatter(zenith, speed, s=2)
		plt.xlabel("Zenith Angle (rad)")
		plt.ylabel("Ejecta Speed (v_esc)")
		plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
		plt.savefig(file[:-5] + ".png")
		plt.clf() # https://stackoverflow.com/questions/8213522/when-to-use-cla-clf-or-close-for-clearing-a-plot-in-matplotlib
