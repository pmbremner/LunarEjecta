#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os

def dist_pdf(a, h, Rm, mu, d, c):
	if d <= c*h:
		return np.sin((a + c*h) / Rm)
	else:
		return np.sin((a + d) / Rm) *( (d / (c*h)) ** -(3.*mu/2. + 1.) + 1.E-3)


def lat_from_dist_bearing(lat1, D_rm, azm):
	return np.arcsin(np.sin(lat1) * np.cos(D_rm) + np.cos(lat1) * np.sin(D_rm) * np.cos(azm))

def lon_from_dist_bearing(lat1, lat2, lon1, D_rm, azm):
	return lon1 + np.arctan2(np.sin(azm) * np.sin(D_rm) * np.cos(lat1), np.cos(D_rm) - np.sin(lat1) * np.sin(lat2))


vdist_pdf = np.vectorize(dist_pdf)

a = float(sys.argv[1]) # m
h = float(sys.argv[2]) # m
N = int(sys.argv[3])   # number of samples to take

asset_lat = float(sys.argv[4]) / 180.*np.pi
asset_lon = float(sys.argv[5]) / 180.*np.pi


Rm = 1737.E3 # m
mu = 0.4
c = 1. / np.tan(1./180.*np.pi)


d = np.logspace(-1, np.log10(Rm*np.pi - a), 10000)

pdf = vdist_pdf(a, h, Rm, mu, d, c)

# https://alpynepyano.github.io/healthyNumerics/posts/sampling_arbitrary_distributions_with_python.html

cdf = np.zeros(np.shape(pdf)[0]+1)
cdf[1:] = np.cumsum(pdf) # prepended with a zero, last element will be 1

cdf /= np.max(cdf) # normalize correctly, so last element is 1



# plt.loglog(d, pdf)
# plt.loglog(d, cdf[1:])
# plt.show()

d_sample   = np.zeros(N)
azm_sample = np.zeros(N)

lat_sample = np.zeros(N)
lon_sample = np.zeros(N)
 

#print(np.argmax(cdf >= 0.9) - 1, d[np.argmax(cdf >= 0.9) - 1])

for i in range(N):
	idx = np.argmax(cdf >= np.random.uniform(0, 1)) - 1
	d_sample[i] = np.random.uniform(d[idx], d[idx+1])

	azm_sample[i] = np.random.uniform(0, 2.*np.pi)

# plt.scatter(d_sample / (Rm*np.pi), azm_sample)
# plt.figure()
# plt.hist(d_sample / (Rm*np.pi), 250)
# plt.yscale('log', nonposy='clip')
# plt.show()
	

lat_sample = lat_from_dist_bearing(asset_lat, d_sample / Rm, azm_sample)

lon_sample = lon_from_dist_bearing(asset_lat, lat_sample, asset_lon, d_sample / Rm, azm_sample)


lon_sample = np.fmod( lon_sample + np.pi, 2.*np.pi ) - np.pi

x_bins = np.linspace(-np.pi, np.pi, 5)
y_bins = np.linspace(-np.pi/2., np.pi/2., 38)


h, xedges, yedges, image = plt.hist2d(lon_sample, lat_sample, bins =[x_bins, y_bins],  norm=mpl.colors.LogNorm(), cmap = plt.cm.gist_ncar)#, cmap = plt.cm.nipy_spectral)

print(np.shape(h))

plt.colorbar()
plt.tight_layout() 
plt.figure()

lat_sample = lat_sample[lon_sample>0]
lon_sample = lon_sample[lon_sample>0]

lat_sample = lat_sample[lon_sample<np.pi/2]
lon_sample = lon_sample[lon_sample<np.pi/2]


lon_sample = lon_sample[lat_sample>0]
lat_sample = lat_sample[lat_sample>0]

lon_sample = lon_sample[lat_sample<5./180.*np.pi]
lat_sample = lat_sample[lat_sample<5./180.*np.pi]

print(lon_sample, lat_sample)

plt.scatter(lon_sample, lat_sample)
plt.ylim(-np.pi/2., np.pi/2.)
plt.xlim(-np.pi, np.pi)

plt.show()