#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import glob

glob_code = sys.argv[1] # vmin*, a-45*

with open('batch_plotFraction.bat', 'w') as bat_file:

	for cur_dir in glob.glob('./' + glob_code + '/'):
		print('cd ' + cur_dir, file=bat_file)
		print('python ../plotTraj.py', file=bat_file)
		print('cd ../', file=bat_file)

	print('python plotFraction.py ' + glob_code, file=bat_file)