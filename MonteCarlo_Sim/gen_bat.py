#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

# python gen_bat.py 12 5. 5457.26E3 1000

N_dist      = int(sys.argv[1])
dmin        = float(sys.argv[2]) # m
dmax        = float(sys.argv[3]) # m
N_integrate = int(sys.argv[4])

d = np.logspace(np.log10(dmin), np.log10(dmax), N_dist)

with open('batch_run.bat', 'w') as bat_file:

	print('@echo off', file=bat_file)
	print('g++ -O2 MC_sim.cpp', file=bat_file)

	for i in range(N_dist):
		print('start /b MC_sim.exe ' + str(N_integrate) + ' ' + str(d[i]) + ' ' + str(i), file=bat_file)
