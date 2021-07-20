#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import glob

def double_power_law(x, A, b, c, d, g, h):
	return np.log10(A / ((10.**x)**b + (10.**x / c)**d) + (10.**x / g)**h)

def comp_intFlux(x, y):

	x0 = x[0:-1]
	x1 = x[1:]
	y0 = y[0:-1]
	y1 = y[1:]

	bi = np.log(y1/y0) / np.log(x1/x0)
	
	flux = np.zeros(np.shape(bi)[0])

	flux[bi == -1.] = y0[bi == -1.] * x0[bi == -1.] * np.log(x1[bi == -1.] / x0[bi == -1.])
	flux[bi != -1.] = (x1[bi != -1.] * y1[bi != -1.] - x0[bi != -1.]* y0[bi != -1.]) / (bi[bi != -1.] + 1)

	flux = np.cumsum(flux)

	return flux[-1] - flux

glob_code = sys.argv[1] # h*r*

fit_vals = np.array([[]])
fdir_vec = np.array([])

fint_flux = np.array([])
fvmin     = np.array([])
fdir = ''


for filename in glob.glob('./' + glob_code + '/FractionEjecta_vs_ImpactDistance.txt'):
	print(filename)
	fdir = filename[:-36]
	fdir_vec = np.append(fdir_vec, fdir)

	distance_m, fraction = np.loadtxt(filename, unpack=True)

	distance_m = distance_m[fraction>0]
	fraction = fraction[fraction>0]

	popt, pcov = curve_fit(double_power_law, np.log10(distance_m - 4.5), np.log10(fraction), bounds=(0, [1.e-1, 0.8, 1.e3, 1.8, 1.e20, 10.5]), p0=(6.e-5, 0.54, 20., 1.1, 50., 1.67))

	fit_vals = np.append(fit_vals, [popt])


	plt.loglog(distance_m - 4.5, fraction, linestyle='-', marker='o', label=fdir)


	# compute integral flux

	int_flux = comp_intFlux(distance_m, fraction)
	fint_flux = np.append(fint_flux, int_flux[0])
	fvmin = np.append(fvmin, float(filename[7:11]))

	np.savetxt('IntegralFlux_' + fdir[2:-2] + '_.txt', np.c_[distance_m[:-1], int_flux])

	plt.loglog(distance_m[:-1] - 4.5, int_flux, linestyle='--', label=fdir+' integral flux')


fit_vals = fit_vals.reshape((np.shape(fdir_vec)[0], 6))
np.savetxt(glob_code[:-1] + '_fraction_fit_params.txt', fit_vals, delimiter=' ')

plt.xlabel('Impact Distance (m)')
plt.ylabel('Fraction of total ejecta')
plt.grid(b=True, which='both')
plt.legend()
plt.savefig(glob_code[:-1] + '_Fraction_vs_dist.png', dpi=1000, bbox_inches='tight', pad_inches=0.1)

plt.figure()
plt.loglog(fvmin, fint_flux, linestyle='-', marker='o')
plt.title('Lander height and radius = ' + fdir[15:27])
plt.xlabel('Minimum ejecta speed (m/s)')
plt.ylabel('Fraction of total ejecta')
plt.grid(b=True, which='both')
plt.savefig(glob_code[:-1] + '_IntFraction_vs_vmin.png', dpi=1000, bbox_inches='tight', pad_inches=0.1)

np.savetxt(glob_code[:-1] + '_IntFraction_vs_vmin.txt', np.c_[fvmin, fint_flux], delimiter=' ')

plt.show()