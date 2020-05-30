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

def min_vel(D):
	return np.sqrt(np.tan(D*np.pi)**2 * (1./np.sin(D*np.pi) - 1.))

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
plt.xlim(0, 90)
# optimal angle for a given velocity line

#D_opt = np.linspace(0.0, 0.4999, Nv)
#gamma_opt = opt_angle(D_opt)

#plt.plot(gamma_opt*180.0/np.pi, norm_speed(D_opt, gamma_opt), '--c')
#ax.annotate('Optimal Angle', xytext=(50, 0.05), xy=(50, 0.05), color='cyan')
D0 = 0.01
D1 = D0*1.5


D_opt = np.linspace(D0, D1, 2)
gamma_opt = opt_angle(D_opt)

x0 = np.array((gamma_opt[0], gamma_opt[0]))*180.0/np.pi
x1 = np.array((gamma_opt[1], gamma_opt[1]))*180.0/np.pi
x2 = np.array((D_opt[0]*np.pi/2., D_opt[0]*np.pi/2.))*180.0/np.pi
y = np.array((0,1))

plt.plot(x0, y)
plt.plot(x1, y)
plt.plot(x2, y)

##################
x = np.array((0,90))
y0 = np.array((min_vel(D0), min_vel(D0)))
y1 = np.array((min_vel(D1), min_vel(D1)))
y2 = np.array((1./np.sqrt(2), 1./np.sqrt(2)))

plt.plot(x, y0)
plt.plot(x, y1)
plt.plot(x, y2)


# Example showing range of speeds and angles between two distances
v_d0 = norm_speed(D0, gamma)
v_d1 = norm_speed(D1, gamma)

plt.plot(gamma*180.0/np.pi, v_d0, 'k')
plt.plot(gamma*180.0/np.pi, v_d1, 'k')

#ax.annotate('D = 0.05', xytext=(33, 0.32), xy=(33, 0.32))
#ax.annotate('D = 0.06', xytext=(33, 0.46), xy=(33, 0.46))

#plt.savefig('Distance_vs_EjectaSpeed_and_ZenithAngle.png', dpi=400)
plt.savefig('Distance_vs_EjectaSpeed_and_ZenithAngle1.png', dpi=400)
plt.show()
