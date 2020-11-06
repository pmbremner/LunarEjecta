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


def downStream(impactZenith): # in degrees out and in
	return 20. + impactZenith * (1.5206 + impactZenith * (-0.036 + impactZenith * 0.0003))

def upStream(impactZenith): # in degrees out and in
	return 20. + impactZenith * (0.129 + impactZenith * (0.0236 + impactZenith * -0.00042))

def genStream(impactZenith, beta_azm):
	return ( downStream(impactZenith * 180./np.pi) + (upStream(impactZenith * 180./np.pi) - downStream(impactZenith * 180./np.pi)) * np.fabs(beta_azm) / np.pi ) * np.pi / 180.

def excZone(impactZenith): # impactZenith in rad, output in rad
	A = downStream(impactZenith * 180./np.pi)
	B = upStream(impactZenith * 180./np.pi)

	if B >= 0.:
		return 0.
	else:
		return -np.pi * B / (A-B)

vexcZone = vectorize(excZone)

def a_power(ejectaZenith, beta_azm):
	zenithMax = genStream(ejectaZenith, beta_azm)
	return np.sqrt(np.cos(zenithMax) / 2.) / np.fabs(np.sin(zenithMax / 2.))

# beta_azm = [-pi, pi]
def ejectaDistribution(impactZenith, ejectaZenith, beta_azm):
	excAngle = vexcZone(impactZenith)

	# compute azimuth term
	if beta_azm > np.pi - excAngle or beta_azm < -np.pi + excAngle:
		azmDist = 0.
	else:
		if impactZenith < np.pi/3.:
			azmDist = (1. + np.cos(beta_azm) * 3.*impactZenith / (2.*np.pi - 3.*impactZenith)) / (2.*np.pi)
		else:
			azmDist = (1. + np.cos(beta_azm) ) / (2.*np.pi)

	# compute zenith term
	a = a_power(ejectaZenith, beta_azm)
	#print(a)
	zenDist =  (1. - np.cos(ejectaZenith)) ** (1./a) * (np.cos(ejectaZenith)) ** a

	return zenDist * azmDist

vejectaDistribution = vectorize(ejectaDistribution)

#N_phi   = 18
#N_theta = 0# defined later
#N_v     = 12 #40
#N_lat   = 37

# d_phi   = 5
# d_theta = 10
# d_v     = 2.37/(N_v-1.) #2
# d_lat   = 5

# phi_min   = -90
# theta_min = 0.
# theta_max = 360.
# v_min     = 0.1 #1
# lat_min   = -90

lon = np.linspace(-np.pi, np.pi, 360)
lat = np.linspace(0, np.pi/2., 180)
Lon,Lat = np.meshgrid(lon,lat)

#filename = sys.argv[1] # '/LatRunData/lat-90/HiDensity/flux_avg.txt'
# preDirectory     = sys.argv[1] # LunarEjectaCode
# ejectaFileName   = sys.argv[2] # run_equator_DSNE_test0.txt

#impZen = float(sys.argv[1]) * np.pi / 180. # impact zenith in degrees to rad

for impZen in np.linspace(0.01, np.pi/2., 15):
	figfilename = 'ejectaDistAt_' + str(impZen * 180/np.pi) + '_degrees_ImpactZenith.png'

	# figfilename = preDirectory + '_'+ ejectaFileName[:-4] + '.png'
	# filename = preDirectory + '/' + ejectaFileName
	# #data = np.loadtxt('RunData/SouthPole/HiDensity/flux_avg.txt', unpack=True) # Equator South45 SouthPole
	# data = np.loadtxt(filename, unpack=True, skiprows = 206)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='mollweide')


	#a = vMEM(Lat.flatten(), Lon.flatten()).reshape(180, 360)

	a = vejectaDistribution(impZen, np.pi/2. - Lat.flatten(), Lon.flatten()).reshape(180, 360)

	im = ax.pcolormesh(Lon, Lat, a, cmap=plt.cm.jet)
	plt.xlabel('Longitude from East CCW [degrees]')
	plt.ylabel('Ejecta Angle from Horizon [degrees]')
	plt.title(r'Ejecta number flux [#/m$^2$/yr]'+'\nmollweide_'+figfilename)

	print(figfilename)
	#plt.show()
	# https://stackoverflow.com/questions/35992492/python-savefig-cuts-off-title
	plt.savefig('mollweide_'+figfilename, dpi = 600, bbox_inches='tight')
	plt.clf()

	# https://stackoverflow.com/questions/9071084/polar-contour-plot-in-matplotlib-best-modern-way-to-do-it
	# https://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html
	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	ax.pcolormesh(np.fmod(Lon + np.pi*2., np.pi*2), (np.pi/2. - Lat)*180./np.pi, a)
	plt.xlabel('Longitude from East CCW [degrees]')
	plt.title(r'Ejecta number flux [#/m$^2$/yr]'+'\npolar_'+figfilename)
	#plt.ylabel('Ejecta Angle from Horizon [degrees]')
	#plt.show()
	plt.savefig('polar_'+figfilename, dpi = 600, bbox_inches='tight')
	plt.clf()
	plt.close('all')
