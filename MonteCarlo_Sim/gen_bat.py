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

# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 1000.
# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 316.
# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 100.
# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 32.
# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 10.
# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 1.
# python ../gen_bat.py 100 4.75 5.45E6 20000 1 4.5 50. 0.

def delta_d(r0, d0, ri):
	return r0 * (np.sqrt((1. + d0/r0)**2 - (1. - ri**2)) - ri - d0/r0)

N_run_min = 200
run_exp = 0.65

N_dist      = int(sys.argv[1]) # per proc
dmin        = float(sys.argv[2]) # m
dmax        = float(sys.argv[3]) # m
N_integrate = int(sys.argv[4])
N_procs     = int(sys.argv[5])
lander_radius = float(sys.argv[6])
lander_height = sys.argv[7]
vmin          = sys.argv[8] # m/s
zang_dist_type  = sys.argv[9] # 0 = iso, 1 = 45 beta, 2 = user-def beta
zang_max        = sys.argv[10]
Nr             = int(sys.argv[11])

lh = float(sys.argv[7])

Nd0 = int(N_dist*N_procs/2)
Nd1 = N_dist - Nd0

d0 = np.logspace(np.log10(dmin), np.log10(3*lh), Nd0, endpoint=False)
d1 = np.logspace(np.log10(3*lh), np.log10(dmax), Nd1)

d = np.append(d0, d1)

r = np.linspace(1., 0., Nr, endpoint=False)



with open('batch_run.bat', 'w') as bat_file:

	#print('@echo off', file=bat_file)
	print('g++ -O2 ../MC_sim.cpp -o MC_sim.exe', file=bat_file)

	for i in range(N_dist*N_procs):
		#print('start /b MC_sim.exe ' + str(N_integrate) + ' ' + str(d[i]) + ' ' + str(i), file=bat_file)
		for k in range(Nr):
			print('MC_sim.exe ' + str(np.maximum(int(N_integrate/(d[i]/d[0])**run_exp), N_run_min)) + ' ' + str(d[i] + delta_d(lander_radius, d[i], r[k])) + ' ' + str(i) + ' ' + str(lander_radius*r) + ' ' + lander_height + ' ' + vmin + ' ' + zang_dist_type + ' ' + zang_max, file=bat_file)

for j in range(N_procs):
	with open('batch_run' + str(j) + '.bat', 'w') as bat_file:
		if j == 0:
			print('g++ -O2 ../MC_sim.cpp -o MC_sim.exe', file=bat_file)

		for i in range(j*N_dist, (j+1)*N_dist):
			#print('start /b MC_sim.exe ' + str(N_integrate) + ' ' + str(d[i]) + ' ' + str(i), file=bat_file)
			for k in range(Nr):
				print('MC_sim.exe ' + str(np.maximum(int(N_integrate/(d[i]/d[0])**run_exp), N_run_min)) + ' ' + str(d[i] + delta_d(lander_radius, d[i], r[k])) + ' ' + str(i) + ' ' + str(lander_radius*r) + ' ' + lander_height + ' ' + vmin + ' ' + zang_dist_type + ' ' + zang_max, file=bat_file)