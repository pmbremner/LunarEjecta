#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def MEM(lat_phi, lon_theta, i_v): # lat and lon in radians, need to be in degrees for MEM file
	i = 0
	#not_found = 1

	phi   = lat_phi * 180. / np.pi
	theta = np.fmod(lon_theta * 180. /np.pi + 360., 360.) 


	while i < 1660:
		if phi >= data[3, i] and phi <= data[4, i] and theta >= data[5, i]  and theta <= data[6, i]:
			if i_v < 0:
				return np.log10(np.sum(data[9:, i])) ## sum over all the speeds
			else:
				return np.log10(data[9 + i_v, i]) #- np.log10(np.sum(data[9:, i])) #
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
eachSpeed        = int(sys.argv[3]) # 0 = sum over all, 1 = each speed


for ilat in range(0, N_lat):
	lat = int(lat_min + d_lat * ilat)
	
	if (eachSpeed == 1 and (ilat == 0 or ilat == 9 or ilat == 18)) or eachSpeed == 0:
		filename = preDirectory + '/lat' + str(lat) + '/' + densityDirectory + '/igloo_avg.txt'
		#data = np.loadtxt('RunData/SouthPole/HiDensity/flux_avg.txt', unpack=True) # Equator South45 SouthPole
		data = np.loadtxt(filename, unpack=True, skiprows = 834)

		a_norm = vMEM(Lat.flatten(), Lon.flatten(), -1).reshape(180, 360)

		for i_vel in range(2,50): # or i_vel = -1 for a sum over all speeds
			if eachSpeed == 0:
				i_vel = -1

			if i_vel < 0:
				figfilename = preDirectory + densityDirectory + str(ilat).zfill(3) + '_lat' + str(lat) + '.png'
			else:
				figfilename = preDirectory + densityDirectory + str(ilat).zfill(3) + '_lat' + str(lat) + '_' + str(v_min + i_vel*d_v) + '-km_s' + '.png'

			fig = plt.figure()
			ax = fig.add_subplot(111, projection='mollweide')


			a = np.divide(vMEM(Lat.flatten(), Lon.flatten(), i_vel).reshape(180, 360), a_norm)

			im = ax.pcolormesh(Lon, Lat, a, cmap=plt.cm.jet)
			plt.xlabel('Longitude from East CCW [degrees]')
			plt.ylabel('Impact Angle from Horizon [degrees]')
			if i_vel < 0:
				plt.title(r'MEM number flux [#/m$^2$/yr]'+'\nmollweide_'+figfilename)
			else:
				plt.title(r'MEM number flux [#/m$^2$/yr] for '+ str(v_min + i_vel*d_v) + ' km/s' +'\nmollweide_'+figfilename)

			print(figfilename)
			# https://stackoverflow.com/questions/35992492/python-savefig-cuts-off-title
			plt.savefig('mollweide_'+figfilename, dpi = 600, bbox_inches='tight')
			plt.clf()

			# https://stackoverflow.com/questions/9071084/polar-contour-plot-in-matplotlib-best-modern-way-to-do-it
			# https://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html
			fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
			ax.pcolormesh(np.fmod(Lon + np.pi*2., np.pi*2), (np.pi/2. - Lat)*180./np.pi, a)
			plt.xlabel('Longitude from East CCW [degrees]')
			if i_vel < 0:
				plt.title(r'MEM number flux [#/m$^2$/yr]'+'\npolar_'+figfilename)
			else:
				plt.title(r'MEM number flux [#/m$^2$/yr] for '+ str(v_min + i_vel*d_v) + ' km/s'+'\npolar_'+figfilename)
			#plt.ylabel('Impact Angle from Horizon [degrees]')

			plt.savefig('polar_'+figfilename, dpi = 600, bbox_inches='tight')
			plt.clf()
			plt.close('all')

			if eachSpeed == 0:
				break