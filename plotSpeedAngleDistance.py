#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def norm_distance(v, gamma):
	return np.arctan2(v**2 * np.sin(2.0*gamma) , (v**2 *(np.cos(2.0*gamma) - 1) + 1)) / np.pi

def norm_speed(D, gamma):
	return (1.0 + np.sin(2.0*gamma)/np.tan(D*np.pi) - np.cos(2.0*gamma))**-0.5

def opt_angle(D):
	return np.arctan2(np.cos(D*np.pi), 1.0 - np.sin(D*np.pi) )

# units of escape speed
vmin = 0
vmax = 1
Nv = 200

gammamin = 0
gammamax = np.pi/2.0
Ngamma = 200

v = np.linspace(vmin, vmax, Nv)
gamma = np.linspace(gammamin, gammamax, Ngamma)

gamma_grid, v_grid = np.meshgrid(gamma, v)

D_grid = norm_distance(v_grid, gamma_grid)

fig, ax = plt.subplots()
plt.contourf(gamma*180.0/np.pi, v, D_grid, 50, cmap='plasma') #inferno
plt.xlabel('Ejecta Zenith Angle [Degrees]')
plt.ylabel('Ejecta Speed [Escape Speed]')
cbar = plt.colorbar()
cbar.set_label('Distance from Impact [Lunar Circumference]')
plt.ylim(0, 1)

# optimal angle for a given velocity line

D_opt = np.linspace(0.0, 0.4999, Nv)
gamma_opt = opt_angle(D_opt)

plt.plot(gamma_opt*180.0/np.pi, norm_speed(D_opt, gamma_opt), '--c')
ax.annotate('Optimal Angle', xytext=(50, 0.05), xy=(50, 0.05), color='cyan')

# dashed line showing the antipodal point
v_antpodal = norm_speed(0.5, gamma)

plt.plot(gamma*180.0/np.pi, v_antpodal, '--r')
# plt.plot(gamma*180.0/np.pi, norm_speed(0.25, gamma), '--g')

ax.annotate('Antipodal Point', xytext=(50, 0.95), xy=(50, 0.95), color='red')

# Example showing range of speeds and angles between two distances
v_d0 = norm_speed(0.05, gamma)
v_d1 = norm_speed(0.06, gamma)

plt.plot(gamma*180.0/np.pi, v_d0, 'k')
plt.plot(gamma*180.0/np.pi, v_d1, 'k')

ax.annotate('D = 0.05', xytext=(33, 0.32), xy=(33, 0.32))
ax.annotate('D = 0.06', xytext=(33, 0.46), xy=(33, 0.46))

plt.savefig('Distance_vs_EjectaSpeed_and_ZenithAngle.png', dpi=400)
plt.show()
