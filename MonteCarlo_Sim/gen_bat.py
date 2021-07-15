#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

# python ../gen_bat.py 13 4.75 3162. 100000 9 4.5 6.25 0.
# python ../gen_bat.py 13 4.75 3162. 100000 9 4.5 12.5 0.
# python ../gen_bat.py 13 4.75 3162. 100000 9 4.5 25. 0.
# python ../gen_bat.py 13 4.75 3162. 100000 9 4.5 50. 0.
# python ../gen_bat.py 13 4.75 3162. 100000 9 4.5 100. 0.
# python ../gen_bat.py 13 4.75 31620. 100000 12 4.5 200. 0.
# python ../gen_bat.py 13 4.75 316200. 100000 18 4.5 400. 0.

# python ../gen_bat.py 39 4.75 31620. 20000 3 4.5 50. 1000.
# python ../gen_bat.py 13 4.75 3162. 10000 9 4.5 50. 100.
# python ../gen_bat.py 39 4.75 3162. 10000 3 4.5 50. 10.

N_dist      = int(sys.argv[1]) # per proc
dmin        = float(sys.argv[2]) # m
dmax        = float(sys.argv[3]) # m
N_integrate = int(sys.argv[4])
N_procs     = int(sys.argv[5])
lander_radius = sys.argv[6]
lander_height = sys.argv[7]
vmin          = sys.argv[8] # m/s

d = np.logspace(np.log10(dmin), np.log10(dmax), N_dist*N_procs)

with open('batch_run.bat', 'w') as bat_file:

	#print('@echo off', file=bat_file)
	print('g++ -O2 ../MC_sim.cpp -o MC_sim.exe', file=bat_file)

	for i in range(N_dist*N_procs):
		#print('start /b MC_sim.exe ' + str(N_integrate) + ' ' + str(d[i]) + ' ' + str(i), file=bat_file)
		print('MC_sim.exe ' + str(np.maximum(int(N_integrate/(d[i]/d[0])**1), 100)) + ' ' + str(d[i]) + ' ' + str(i) + ' ' + lander_radius + ' ' + lander_height + ' ' + vmin, file=bat_file)

for j in range(N_procs):
	with open('batch_run' + str(j) + '.bat', 'w') as bat_file:
		if j == 0:
			print('g++ -O2 ../MC_sim.cpp -o MC_sim.exe', file=bat_file)

		for i in range(j*N_dist, (j+1)*N_dist):
			#print('start /b MC_sim.exe ' + str(N_integrate) + ' ' + str(d[i]) + ' ' + str(i), file=bat_file)
			print('MC_sim.exe ' + str(np.maximum(int(N_integrate/(d[i]/d[0])**1), 100)) + ' ' + str(d[i]) + ' ' + str(i) + ' ' + lander_radius + ' ' + lander_height + ' ' + vmin, file=bat_file)