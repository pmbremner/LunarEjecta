#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import random as rd


r_ROI = float(sys.argv[1])
D0    = float(sys.argv[2])
D1    = float(sys.argv[3])
dphi  = float(sys.argv[4]) * np.pi / 180.

# N_ROI_r   = int(sys.argv[5])
# N_ROI_phi = int(sys.argv[6])
# N_arc_r   = int(sys.argv[7])
# N_arc_phi = int(sys.argv[8])
N_ROI   = int(sys.argv[5])
N_arc   = int(sys.argv[6])

p_ROI_x = np.zeros(N_ROI)
p_ROI_y = np.zeros(N_ROI)
p_arc_x = np.zeros(N_arc)
p_arc_y = np.zeros(N_arc)

# #print("define ROI points x and y\n\n")
# # define ROI points x and y
# for i in range(0, N_ROI_r):
# 	r = (i + 1) * r_ROI / float(N_ROI_r)
# 	for j in range(0, N_ROI_phi):
# 		phi = j * 2.*np.pi / float(N_ROI_phi)

# 		idx = j + i * N_ROI_phi

# 		p_ROI_x[idx] = r * np.cos(phi)
# 		p_ROI_y[idx] = r * np.sin(phi)
		
# 		#print(p_ROI_x[idx], p_ROI_y[idx], np.sqrt(p_ROI_x[idx]**2 + p_ROI_y[idx]**2))

# #print("define arc points x and y\n\n")
# # define arc points x and y
# for i in range(0, N_arc_r):
# 	r = D0 + (D1 - D0) * i / float(N_arc_r-1.)
# 	for j in range(0, N_arc_phi):
# 		phi = dphi * j / float(N_arc_phi)

# 		idx = j + i * N_arc_phi

# 		p_arc_x[idx] = r * np.cos(phi)
# 		p_arc_y[idx] = r * np.sin(phi)

# 		#print(p_arc_x[idx], p_arc_y[idx], np.sqrt(p_arc_x[idx]**2 + p_arc_y[idx]**2))


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

		ang = np.mod(np.arctan2(y, x) + 2*np.pi, 2*np.pi)
		ang_array[idx] = ang

		#print(d, np.mod(ang + 2.*np.pi, 2*np.pi))

ang_array = ang_array*180./np.pi - 180.

x_bins = np.linspace(np.min(d_array), np.max(d_array), 20) 
y_bins = np.linspace(np.min(ang_array), np.max(ang_array), 20)

plt.hist2d(d_array, ang_array, bins =[x_bins, y_bins])
plt.figure()
plt.hist(d_array)
plt.figure()
plt.hist(ang_array)
plt.show()