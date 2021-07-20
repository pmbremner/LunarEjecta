#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import glob

#moon_x, moon_y, moon_z = np.loadtxt('moon.txt', unpack=True)

##fig, (ax1, ax2) = plt.subplots(1, 2)

Rm = 1737.1E3

dist_array = np.loadtxt('batch_run.bat', unpack=True, skiprows=1, usecols=(2))

N_run = np.shape(dist_array)[0]
hits_array = np.zeros(N_run)

track_array = np.array([])
speed_array = np.array([])

print('N_run = ', N_run)

for filename in glob.glob('dists*.txt'):
	print(filename)

	Ni = int(filename[6:-4])
	print(Ni)
	#track_x, track_z          = np.loadtxt('tracks.txt', unpack=True, usecols=(0,1))
	track_x, track_z, track_dist, track_init_speed, track_ej_zang, sum_miss, sum_hit, N_miss, N_hit = np.loadtxt(filename, unpack=True)
	N = np.shape(sum_miss)[0]
	i = np.linspace(0, N, N)

	N_tot = N_miss + N_hit

	hits_array[Ni] = sum_hit[-1]/N_tot[-1]

	### track_array = np.append(track_array, track_dist)
	### speed_array = np.append(speed_array, track_init_speed)

	# plt.plot(moon_x, moon_z)
	# plt.plot(track_x, track_z, 'o', markersize=1)

	# plt.figure()
	##ax1.plot(track_ej_zang*180./np.pi, track_init_speed, 'o', markersize=0.5, label=filename)

	##ax2.loglog(i, sum_hit/N_tot,  label=filename + ' sum hit')
	##ax2.loglog(i, sum_miss/N_tot, label=filename + ' sum miss')
	#ax2.loglog(i, (sum_hit+sum_miss)/N_tot, label=filename + ' sum total')
	
##ax1.set_ylim([0, 1])
##ax1.set_xlim([-90, 90])
##ax1.legend()
##ax2.legend(loc='lower right')
##ax2.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()

np.savetxt('FractionEjecta_vs_ImpactDistance.txt', np.c_[dist_array, hits_array])
###np.savetxt('distance_vs_speed.txt', np.c_[track_array, speed_array])


plt.loglog(dist_array-4.5, hits_array, linestyle='-', marker='o')
plt.xlabel('Impact Distance (m)')
plt.ylabel('Fraction of total ejecta')
plt.grid(b=True, which='both')
plt.savefig('Fraction_vs_dist.png', dpi=800, bbox_inches='tight', pad_inches=0.1)

# plt.figure()
# plt.hist2d(track_array*Rm-4.5, speed_array, norm=mpl.colors.LogNorm(), bins=35)
# plt.xlabel('Impact Distance (m)')
# plt.ylabel(r'Ejecta Speed (v$_{esc}$)')
# plt.grid(b=True, which='both')
# plt.savefig('distance_vs_speed.png', dpi=800, bbox_inches='tight', pad_inches=0.1)

#plt.show()