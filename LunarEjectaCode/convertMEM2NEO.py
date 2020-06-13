#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# read in speed distribution
NEOfractionFile = "../NEA_Brown/vmass.txt";

# we will assume the speeds range from 1 to 69
NEOfacData = np.loadtxt(NEOfractionFile, unpack=True)

# read file to renormalize (assumes we have lat files from -90 to 90 in 5 degree increments, including 0)
for fi in range(0, 37):
	fromfilename = "../LatRunData/lat" + str(int(-90 + 180*fi/36)) + "/HiDensity/igloo_avg.txt"

	tofilename = "../LatNEOData/lat" + str(int(-90 + 180*fi/36)) + "/HiDensity/igloo_avg.txt"

	print(fromfilename)

	MEMheader = ''
	i = 0                  
	with open (fromfilename, 'rt') as myfile:
	    for myline in myfile:
	    	i = i + 1
	    	if i <= 8:
	        	MEMheader = MEMheader + myline
	#print(MEMheader)                          

	data = np.loadtxt(fromfilename, unpack=True)

	# for each speed column, assuming we have 35 rows in the vmass file
	# (speed centers ranging from 1 to 69, or speed domain of 0 to 70 km/s)
	for iv in range(0,35):
		vfrac = NEOfacData[1][iv]
		MEMsum1 = np.sum(data[9 + 2*iv    ])
		MEMsum2 = np.sum(data[9 + 2*iv + 1])

		if MEMsum1 == 0.:
			vfrac = 0.;
			MEMsum1 = 1.;
		if MEMsum2 == 0.:
			vfrac = 0.;
			MEMsum2 = 1.;

		# renormalize the speed columns
		renorm1 = vfrac / MEMsum1
		renorm2 = vfrac / MEMsum2
		
		# Need to divide by 2 because the two bins should add to the total fraction of the NEO bin
		data[9 + 2*iv    ] = renorm1 * data[9 + 2*iv    ] / 2.
		data[9 + 2*iv + 1] = renorm2 * data[9 + 2*iv + 1] / 2.
	
	# The last part of the speed distribution is set to zero (no NEO's here)
	for iv in range(70,80):
		data[9 + iv] = 0. * data[9 + iv]

	# save the new file
	MEM_fmt = '%d', '%d', '%d', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e', '%e'

	np.savetxt(tofilename, data.T, header=MEMheader[:-2], comments='', fmt=MEM_fmt)
	print(tofilename)