#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt


# isotropic zenith angle distribution
# python ./vmin_study_gen.py 200 4.55 5.45E6 50000 1 4.5 50. 1. 1000. 7 0 0. vmin

# peak zenith angle centered around 45 degrees
# python ./vmin_study_gen.py 200 4.55 5.45E6 50000 1 4.5 50. 1. 1000. 7 1 1.5708 a-45

N_dist      = sys.argv[1] # per proc
dmin        = sys.argv[2] # m
dmax        = sys.argv[3] # m
N_integrate = sys.argv[4]
N_procs     = sys.argv[5]
lander_radius = sys.argv[6]
lander_height = sys.argv[7]
vmin          = float(sys.argv[8])
vmax          = float(sys.argv[9])
Nv            = int(sys.argv[10])
zang_dist_type  = sys.argv[11] # 0 = iso, 1 = 45 beta, 2 = user-def beta
zang_max        = sys.argv[12]
run_name      = sys.argv[13]
Nr            = sys.argv[14]

v = np.logspace(np.log10(vmin), np.log10(vmax), Nv)
v_cur = 0.

with open(run_name + '_study_run.bat', 'w') as bat_file:
	for i in range(Nv+1):
		if i == 0:
			v_cur = 0.
		else:
			v_cur = v[i-1]

		dir_name = run_name + '-' + str(int(v_cur)).zfill(4) + 'mps-h' + lander_height.zfill(4) + 'm-r' + lander_radius + 'm'


		print('mkdir ' + dir_name, file=bat_file)
		print('cd ./' + dir_name, file=bat_file)
		print('python ../gen_bat.py ' + N_dist + ' ' + dmin + ' ' + dmax + ' ' + N_integrate + ' ' + N_procs + ' ' + lander_radius + ' ' + lander_height + ' ' + str(v_cur) + ' ' + zang_dist_type + ' ' + zang_max + ' ' + Nr, file=bat_file)

		for j in range(int(N_procs)):
			print('start batch_run' + str(j) + '.bat', file=bat_file)

		print('cd ../', file=bat_file)
	