#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

import glob

for file in glob.glob("ejecta_vs_dist*.txt"):

	d, frac = np.loadtxt(file, unpack=True)

	N =np.shape(d)[0]

	# frac_tot = np.zeros(N)

	# for i in range(N):
	# 	frac_tot[i] = frac[i] + frac[N-1-i]


	plt.loglog(d/(2.*np.pi), frac, label=file)
	# plt.xlim([1E-6, np.pi])
	plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
	plt.legend()
	plt.xlabel("Lunar Circumference")
	plt.ylabel("Fraction of Ejecta")


plt.show()