#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://matplotlib.org/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
from matplotlib import ticker, cm

N_phi   = 36
N_theta = 72 
N_v     = 80 #40
N_lat   = 37

d_phi   = 5
d_theta = 5
d_v     = 1 #2
d_lat   = 5

phi_min   = -90
theta_min = 0
v_min     = 0 #1
lat_min   = -90

#filename = sys.argv[1] # '/LatRunData/lat-90/HiDensity/flux_avg.txt'
preDirectory     = sys.argv[1] # LatRunData
densityDirectory = sys.argv[2] # HiDensity LoDensity


for ilat in range(0, N_lat):
	lat = int(lat_min + d_lat * ilat)
	figfilename = densityDirectory + str(ilat).zfill(3) + '_lat' + str(lat) + '.png'
	filename = preDirectory + '/lat' + str(lat) + '/' + densityDirectory + '/flux_avg.txt'
	#data = np.loadtxt('RunData/SouthPole/HiDensity/flux_avg.txt', unpack=True) # Equator South45 SouthPole
	data = np.loadtxt(filename, unpack=True) # Equator South45 SouthPole

	v = np.linspace(v_min, v_min + d_v * (N_v-1), N_v)
	phi = np.linspace(0, phi_min + d_phi * (N_phi-1), int(N_phi/2))
	#print(v)
	# data2 = np.reshape(data[:,int(N_phi*N_theta/2):], (N_v+2, N_theta, int(N_phi/2)))
	# print(np.shape(data2))
	# print(np.shape(np.sum(data2, axis=1)))
	# print(data2[:,:,0])
	# np.savetxt("test.txt", data2[:,:,0])

	# print(data2[0,0,:])
	# plt.plot(v, np.sum(data2[2:,:,17], axis=1))
	# plt.show()

	# plt.plot(v, data2[2:,1,0])
	# plt.show()

	data2d = np.zeros((int(N_phi/2), N_v))
	#print(np.shape(data2d))

	for i in range(int(N_phi/2), N_phi):
		##plt.plot(v, np.sum(data[2:,int(i*N_theta):int((i+1)*N_theta)], axis=1), label=str(data[0, int(i*N_theta)]))
		#print(np.shape(np.sum(data[2:,int(i*N_theta):int((i+1)*N_theta)], axis=1)))
		#print(np.shape(data2d[i-int(N_phi/2),:]))
		data2d[i-int(N_phi/2)] = np.sum(data[2:,int(i*N_theta):int((i+1)*N_theta)], axis=1)

	#plt.yscale('log')
	##plt.legend()
	##plt.show()
	#print(np.shape(v), np.shape(phi), np.shape(data2d))
	plt.contourf(v, phi, data2d/(data2d.max()), 20)
	plt.colorbar(label='Max flux [#/m^2/year] = ' + str("{0:.3E}".format(data2d.max())))
	plt.xlabel('Impact Speed [km/s]')
	plt.ylabel('Impact Angle from Horizon [degrees]')
	plt.title(filename)
	#plt.ylim(bottom=1e-5)
	print(figfilename)
	plt.savefig(figfilename)
	plt.clf()