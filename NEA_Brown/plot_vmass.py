#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

speed, frac = np.loadtxt('vmass.txt', unpack=True)

plt.bar(speed, frac, width=2)
plt.xlabel('v (km/s)')
plt.ylabel('fraction per km/s')
plt.show()
