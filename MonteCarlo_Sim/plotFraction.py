#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import glob

def double_power_law(x, A, b, c, d):
	return np.log10(A / ((10.**x)**b + (10.**x / c)**d))


glob_code = sys.argv[1] # h*r*

fit_vals = np.array([[]])
fdir_vec = np.array([])

for filename in glob.glob('./' + glob_code + '/FractionEjecta_vs_ImpactDistance.txt'):
	print(filename)
	fdir = filename[:-36]
	fdir_vec = np.append(fdir_vec, fdir)

	distance_m, fraction = np.loadtxt(filename, unpack=True)

	popt, pcov = curve_fit(double_power_law, np.log10(distance_m - 4.5), np.log10(fraction), bounds=(0, [1.e-1, 0.8, 1.e3, 2.5]), p0=(6.e-5, 0.54, 20., 1.67))

	fit_vals = np.append(fit_vals, [popt])


	plt.loglog(distance_m - 4.5, fraction, linestyle='-', marker='o', label=fdir)

fit_vals = fit_vals.reshape((np.shape(fdir_vec)[0], 4))
np.savetxt('fraction_fit_params.txt', fit_vals, delimiter=' ')

plt.xlabel('Impact Distance (m)')
plt.ylabel('Fraction of total ejecta')
plt.grid(b=True, which='both')
plt.legend()
plt.show()