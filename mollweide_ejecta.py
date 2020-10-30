#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def MEM(lat_phi, lon_theta): # lat and lon in radians, need to be in degrees for MEM file
	i = 0
	#not_found = 1

	phi   = lat_phi * 180. / np.pi
	theta = np.fmod(lon_theta * 180. /np.pi + 360., 360.) 


	while i < 1660:
		if phi >= data[3, i] and phi <= data[4, i] and theta >= data[5, i]  and theta <= data[6, i]:
			return np.log10(np.sum(data[9:, i])) ## sum over all the speeds
		i = i + 1
	return 0.

vMEM = vectorize(MEM)

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

lon = np.linspace(-np.pi, np.pi,360)
lat = np.linspace(0, np.pi/2.,180)
Lon,Lat = np.meshgrid(lon,lat)

#filename = sys.argv[1] # '/LatRunData/lat-90/HiDensity/flux_avg.txt'
preDirectory     = sys.argv[1] # LatRunData
densityDirectory = sys.argv[2] # HiDensity LoDensity


for ilat in range(0, N_lat):
	lat = int(lat_min + d_lat * ilat)
	figfilename = preDirectory + densityDirectory + str(ilat).zfill(3) + '_lat' + str(lat) + '.png'
	filename = preDirectory + '/lat' + str(lat) + '/' + densityDirectory + '/igloo_avg.txt'
	#data = np.loadtxt('RunData/SouthPole/HiDensity/flux_avg.txt', unpack=True) # Equator South45 SouthPole
	data = np.loadtxt(filename, unpack=True, skiprows = 834)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='mollweide')


	a = vMEM(Lat.flatten(), Lon.flatten()).reshape(180, 360)

	im = ax.pcolormesh(Lon, Lat, a, cmap=plt.cm.jet)
	plt.xlabel('Longitude from East CCW [degrees]')
	plt.ylabel('Impact Angle from Horizon [degrees]')
	plt.title(r'MEM number flux [#/m$^2$/yr]'+'\nmollweide_'+figfilename)

	print(figfilename)
	# https://stackoverflow.com/questions/35992492/python-savefig-cuts-off-title
	plt.savefig('mollweide_'+figfilename, dpi = 600, bbox_inches='tight')
	plt.clf()

	# https://stackoverflow.com/questions/9071084/polar-contour-plot-in-matplotlib-best-modern-way-to-do-it
	# https://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html
	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	ax.pcolormesh(np.fmod(Lon + np.pi*2., np.pi*2), (np.pi/2. - Lat)*180./np.pi, a)
	plt.xlabel('Longitude from East CCW [degrees]')
	plt.title(r'MEM number flux [#/m$^2$/yr]'+'\npolar_'+figfilename)
	#plt.ylabel('Impact Angle from Horizon [degrees]')

	plt.savefig('polar_'+figfilename, dpi = 600, bbox_inches='tight')
	plt.clf()
	plt.close('all')
