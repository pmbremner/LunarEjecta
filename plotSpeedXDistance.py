#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def norm_distance_x(v, x):
	return np.arctan2(2. * v**2 * (1.-x)*np.sqrt(x*(2.-x)), 1. - 2. * v**2 *x*(2.-x)) / np.pi

def norm_speed_x(D, x):
	return (2.*(1.-x)*np.sqrt(x*(2.-x))/np.tan(D*np.pi) + 2.*x*(2.-x))**-0.5

# def cot_angle_p(v, D):
# 	b2 = v**2 / np.tan(D*np.pi)
# 	return b2 + np.sqrt(b2**2 + 2.0*v**2 - 1.0)

def norm_X_p(v, D):
	tan2 = (np.tan(D*np.pi))**2
	return 1. - np.sqrt( (v**2 + tan2*(2.*v**2 - 1.) + np.sqrt(v**4 + tan2*(2.*v**2-1.)))
		/ (2.*v**2*(1.+tan2)) )

def norm_X_m(v, D):
	tan2 = (np.tan(D*np.pi))**2
	return 1. - np.sqrt( (v**2 + tan2*(2.*v**2 - 1.) - np.sqrt(v**4 + tan2*(2.*v**2-1.)))
		/ (2.*v**2*(1.+tan2)) )

def v_min(D):
	return 1. - np.sqrt((1. - np.sin(D*np.pi)) / 2.0)

# units of escape speed
vmin = 0.
vmax = 1.
Nv = 200

xmin = 0.
xmax = 1.
Nx = 200

v = np.linspace(vmin, vmax, Nv)
x = np.linspace(xmin, xmax, Nx)

x_grid, v_grid = np.meshgrid(x, v)

D_grid = norm_distance_x(v_grid, x_grid)

fig, ax = plt.subplots()
plt.contourf(x, v, D_grid, 50, cmap='plasma') #inferno
plt.xlabel('Transformed Ejecta Zenith Angle [x = 1 - cos(Zenith Angle)]')
plt.ylabel('Ejecta Speed [Escape Speed]')
cbar = plt.colorbar()
cbar.set_label('Distance from Impact [Lunar Circumference]')
plt.ylim(0, 1)

# set of horizontal lines
y_lines = np.linspace(vmin, vmax, 10)
for i in range(0,10):
	plt.axhline(y=y_lines[i], linestyle='--', color='black')


# Example showing range of speeds and angles between two distances
D0 = 0.01
D1 = 0.015
D_array = np.array(((D0,D1)))

v_d0 = norm_speed_x(D0, x)
v_d1 = norm_speed_x(D1, x)

plt.plot(x, v_d0, 'k')
plt.plot(x, v_d1, 'k')

#####
plt.plot(norm_X_p(y_lines, D0), y_lines, 'ko')
plt.plot(norm_X_m(y_lines[y_lines < 0.5*np.sqrt(2)], D0), y_lines[y_lines < 0.5*np.sqrt(2)], 'ko')

plt.plot(norm_X_p(y_lines, D1), y_lines, 'ko')
plt.plot(norm_X_m(y_lines[y_lines < 0.5*np.sqrt(2)], D1), y_lines[y_lines < 0.5*np.sqrt(2)], 'ko')

v_min_array = v_min(D_array)
print(v_min_array)

plt.plot(v_min_array, norm_speed_x(D_array, v_min_array), 'go')

#print(norm_X_p(y_lines, D0))

plt.show()