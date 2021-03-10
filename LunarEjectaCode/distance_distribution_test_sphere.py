#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl # https://stackoverflow.com/questions/23309272/matplotlib-log-transform-counts-in-hist2d
import sys
import random as rd


def get_lat_from_equator(D, bearing):
	return np.arcsin(np.sin(D) * np.cos(bearing))

def get_lon_from_equator(D, bearing):
	return np.arctan2(np.sin(D) * np.sin(bearing), np.cos(D))


def get_lat_from_pole(D, bearing):
	return np.pi/2. - D

def get_lon_from_pole(D, bearing):
	return bearing

# D in units of radius, angles in radians
def get_lat_from_ROI(D, bearing, lat_R, lon_R):
	return np.arcsin(np.sin(lat_R)*np.cos(D) + np.cos(lat_R)*np.sin(D)*np.cos(bearing))

def get_lon_from_ROI(D, bearing, lat_R, lon_R):
	return lon_R + np.arctan2(np.sin(bearing)*np.sin(D), np.cos(D)*np.cos(lat_R) - np.sin(D)*np.sin(lat_R)*np.cos(bearing))


r_ROI   = float(sys.argv[1])
lat_ROI = float(sys.argv[2]) * np.pi / 180.
lon_ROI = float(sys.argv[3]) * np.pi / 180.

D0     = float(sys.argv[4])
D1     = float(sys.argv[5])
dphi   = float(sys.argv[6]) * np.pi / 180.
phi_offset = float(sys.argv[7]) * np.pi / 180.

N_ROI   = int(sys.argv[8])
N_arc   = int(sys.argv[9])
radius  = float(sys.argv[10])

p_ROI_lat = np.zeros(N_ROI)
p_ROI_lon = np.zeros(N_ROI)
p_arc_lat = np.zeros(N_arc)
p_arc_lon = np.zeros(N_arc)

# compute lat and lon points in ROI (assuming north pole for simplicity)
for i in range(0, N_ROI):
	z   = rd.uniform(0., 2.*np.sin(r_ROI/(2.*radius))**2) # units of radius
	D   = 2.*np.arcsin(np.sqrt(z/2.)) # units of radius
	phi = rd.uniform(0., 2.*np.pi) # bearing

	p_ROI_lat[i] = get_lat_from_ROI(D, phi, lat_ROI, lon_ROI)
	p_ROI_lon[i] = get_lon_from_ROI(D, phi, lat_ROI, lon_ROI)

	#print(p_ROI_lat[i], p_ROI_lon[i])

# compute lat and lon points in arc
for i in range(0, N_arc):
	z   = rd.uniform(2.*np.sin(D0/(2.*radius))**2, 2.*np.sin(D1/(2.*radius))**2) # units of radius
	D   = 2.*np.arcsin(np.sqrt(z/2.)) # units of radius
	phi = rd.uniform(-dphi/2. + phi_offset, dphi/2. + phi_offset) # bearing

	p_arc_lat[i] = get_lat_from_pole(D, phi)
	p_arc_lon[i] = get_lon_from_pole(D, phi)

	#print(p_arc_lat[i], p_arc_lon[i])

# compute all mutual distances

d_array = np.zeros(N_ROI * N_arc)
ang_array = np.zeros(N_ROI * N_arc)

for i in range(0, N_ROI):
	for j in range(0, N_arc):
		
		idx = j + i*N_arc

		dlat = p_ROI_lat[i] - p_arc_lat[j]
		dlon = p_ROI_lon[i] - p_arc_lon[j]

		a = np.sin(dlat/2.)**2 + np.cos(p_ROI_lat[i])*np.cos(p_arc_lat[j])*np.sin(dlon/2.)**2

		d_array[idx] = 2.*np.arcsin(np.sqrt(a))

		ang_array[idx] = np.arctan2(np.sin(dlon)*np.cos(p_arc_lat[j]), np.cos(p_ROI_lat[i])*np.sin(p_arc_lat[j]) - np.sin(p_ROI_lat[i])*np.cos(p_arc_lat[j])*np.cos(dlon))

ang_array = ang_array * 180./np.pi


x_N = 10 #int(radius*(np.max(d_array) - np.min(d_array)) / (D1-D0) * 5)
y_N = 10#int((np.max(ang_array) - np.min(ang_array)) / dphi)

print(x_N, y_N)

# bearing vs distance
x_bins = np.linspace(np.min(d_array), np.max(d_array), x_N) 
y_bins = np.linspace(np.min(ang_array), np.max(ang_array), y_N)

plt.hist2d(d_array, ang_array, bins =[x_bins, y_bins], norm=mpl.colors.LogNorm())

# lat lon coords
x_bins = np.linspace(-np.pi/2, np.pi/2, 180) 
y_bins = np.linspace(-np.pi, np.pi, 360)

plt.figure()
plt.hist2d(np.concatenate([p_ROI_lon, p_arc_lon]), np.concatenate([p_ROI_lat, p_arc_lat]), bins =[y_bins, x_bins], norm=mpl.colors.LogNorm())
# plt.figure()
# plt.hist(d_array)
# plt.figure()
# plt.hist(ang_array)
plt.show()
