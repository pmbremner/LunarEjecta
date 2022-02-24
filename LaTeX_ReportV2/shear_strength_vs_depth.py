#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

def strength_vol_parabola(x):
	if x < 50:
		return 1.28805 + x/136.65
	elif x >= 50 and x < 250:
		return 1.14353 - 118.333/x**2 + 7.16294/x + x/120.63
	else:
		return 1.9423 + 16902.4/x**2 - 194.807/x + x/138.18


vstrength_vol_parabola = np.vectorize(strength_vol_parabola)

def int_ext_data(x_data, y_data, x):
	x0 = 0.
	y0 = 0.
	x1 = 0.
	y1 = 0.

	if x < x_data[0]:
		x0 = x_data[0]
		x1 = x_data[2]

		y0 = y_data[0]
		y1 = y_data[2]

	elif x > x_data[-1]:
		x0 = x_data[-3]
		x1 = x_data[-1]

		y0 = y_data[-3]
		y1 = y_data[-1]

	else:
		x0 = x_data[x >= x_data][-1]
		y0 = y_data[x >= x_data][-1]

		x1 = x_data[x <= x_data][0]
		y1 = y_data[x <= x_data][0]


	return (y1 - y0) * (x - x0) / (x1 - x0) + y0




d_cm_low, strength_low   = np.loadtxt('regolith_shear_strength_lower_bound_fig9.26.txt', unpack=True, delimiter=',')
d_cm_high, strength_high = np.loadtxt('regolith_shear_strength_upper_bound_fig9.26.txt', unpack=True, delimiter=',')


n = 600

d_x = np.linspace(0., 600, n)

f0 = np.zeros(n)
f1 = np.zeros(n)
f_avg = np.zeros(n)


for i in range(n):
	f0[i] = int_ext_data(d_cm_low, strength_low, d_x[i])
	f1[i] = int_ext_data(d_cm_high, strength_high, d_x[i])
	f_avg[i] = (f0[i] + f1[i]) / 2.
	
print(f_avg[1], f_avg[50], f_avg[250], f_avg[300])
print(1./((f_avg[50] - f_avg[1])/(d_x[50] - d_x[1])), f_avg[1])
print(1./((f_avg[250] - f_avg[50])/(d_x[250] - d_x[50])), f_avg[50])
print(1./((f_avg[300] - f_avg[250])/(d_x[300] - d_x[250])), f_avg[250])

# https://matplotlib.org/stable/gallery/color/named_colors.html
plt.plot(d_x[d_x < 335.], f0[d_x < 335.], label='Apollo and Lunokhod lower bound')
plt.plot(d_x[d_x >= 335.], f0[d_x >= 335.], color='tab:blue', linestyle='--')
plt.plot(d_x[d_x < 165.], f1[d_x < 165.], label='Apollo and Lunokhod upper bound')
plt.plot(d_x[d_x >= 165.], f1[d_x >= 165.], color='tab:orange', linestyle='--')
plt.plot(d_x[d_x < 165.], f_avg[d_x < 165.], label='Apollo and Lunokhod average')
plt.plot(d_x[(d_x >= 165.) & (d_x < 335)], f_avg[(d_x >= 165.) & (d_x < 335)], color='tab:green', linestyle='-.')
plt.plot(d_x[d_x >= 335.], f_avg[d_x >= 335.], color='tab:green', linestyle='--')

plt.plot(d_x, vstrength_vol_parabola(d_x), label='Volume-averaged (parabolic crater)')

plt.xlim(0., 600.)
plt.xlabel(r'Depth Below Lunar Surface, $z$ (cm)')
plt.ylabel(r'Shear Strength, $\tau$ (kPa)')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.legend()
plt.savefig('shear_strength_vs_depth.png', bbox_inches='tight', dpi=600)
plt.show()