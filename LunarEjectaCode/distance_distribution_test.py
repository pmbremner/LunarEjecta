#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl # https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
import sys
import random as rd


r_ROI = float(sys.argv[1])
D0    = float(sys.argv[2])
D1    = float(sys.argv[3])
dphi  = float(sys.argv[4]) * np.pi / 180.

N_ROI   = int(sys.argv[5])
N_arc   = int(sys.argv[6])

p_ROI_x = np.zeros(N_ROI)
p_ROI_y = np.zeros(N_ROI)
p_arc_x = np.zeros(N_arc)
p_arc_y = np.zeros(N_arc)


for i in range(0, N_ROI):
	r   = rd.uniform(0., r_ROI)
	phi = rd.uniform(0., 2.*np.pi)

	p_ROI_x[i] = r * np.cos(phi)
	p_ROI_y[i] = r * np.sin(phi)


for i in range(0, N_arc):
	r   = rd.uniform(D0, D1)
	phi = rd.uniform(0., dphi)

	p_arc_x[i] = r * np.cos(phi)
	p_arc_y[i] = r * np.sin(phi)

# compute all mutual distances

d_array = np.zeros(N_ROI * N_arc)
ang_array = np.zeros(N_ROI * N_arc)

for i in range(0, N_ROI):
	
	for j in range(0, N_arc):
		x = p_ROI_x[i] - p_arc_x[j]
		y = p_ROI_y[i] - p_arc_y[j]

		idx = j + i*N_arc

		d = np.sqrt( x**2 + y**2 )
		d_array[idx] = d

		ang = np.abs((np.mod(np.arctan2(y, x) + 2*np.pi, 2*np.pi) - dphi/2.)*180./np.pi - 180.)
		ang_array[idx] = ang

		#print(d, np.mod(ang + 2.*np.pi, 2*np.pi))

#ang_array = ang_array*180./np.pi - 180.

x_bins = np.linspace(np.min(d_array), np.max(d_array), 50) 
y_bins = np.linspace(np.min(ang_array), np.max(ang_array), 50)

plt.hist2d(d_array, ang_array, bins =[x_bins, y_bins], norm=mpl.colors.LogNorm())
plt.figure()
plt.hist(d_array)
plt.figure()
plt.hist(ang_array)
plt.show()