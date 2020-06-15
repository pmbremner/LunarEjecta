#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://matplotlib.org/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
from matplotlib import ticker, cm

N_phi   = 36
N_theta = 0# defined later
N_v     = 80 #40
N_lat   = 37

d_phi   = 5
d_theta = 5
d_v     = 1 #2
d_lat   = 5

phi_min   = -90
theta_min = 0.
theta_max = 360.
v_min     = 0. #1
lat_min   = -90


#filename = sys.argv[1] # '/LatRunData/lat-90/HiDensity/flux_avg.txt'
preDirectory     = sys.argv[1] # LatRunData
densityDirectory = sys.argv[2] # HiDensity LoDensity


for ilat in range(0, N_lat):
	lat = int(lat_min + d_lat * ilat)
	figfilename = preDirectory + densityDirectory + str(ilat).zfill(3) + '_lat' + str(lat) + '.png'
	filename = preDirectory + '/lat' + str(lat) + '/' + densityDirectory + '/igloo_avg.txt'
	#data = np.loadtxt('RunData/SouthPole/HiDensity/flux_avg.txt', unpack=True) # Equator South45 SouthPole
	data = np.loadtxt(filename, unpack=True, skiprows = 834)

	print(np.shape(data))

	v = np.linspace(v_min, v_min + d_v * (N_v-1), N_v)
	phi = np.linspace(0, phi_min + d_phi * (N_phi-1), int(N_phi/2))

	data2d = np.zeros((int(N_phi/2), N_v))

	Ndata = np.shape(data)[1]

	i = 0
	j = 0
	print("\n")
	while i < Ndata:
		N_theta = int(np.round(360./data[6][i]))

		print(np.shape(data[9:, i:i+N_theta-1]), np.shape(data2d[:]))
		data2d[j] = data2d[j] + np.sum(data[9:, i:i+N_theta-1], axis=1)

		j = j + 1
		i = i + N_theta

	plt.contourf(v, phi, data2d/(data2d.max()), 20)
	plt.colorbar(label='Max flux [#/m^2/year] = ' + str("{0:.3E}".format(data2d.max())))
	plt.xlabel('Impact Speed [km/s]')
	plt.ylabel('Impact Angle from Horizon [degrees]')
	plt.title(filename)
	#plt.ylim(bottom=1e-5)
	print(figfilename)
	plt.savefig(figfilename)
	plt.clf()

