#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/2507808/how-to-check-whether-a-file-is-empty-or-not
import os


idx, cdf = np.loadtxt('cdf_example.txt', unpack=True)

N = cdf.size

print(N)

plt.plot(cdf[0:N//3-1], label='MEM_hi_fluxes')
plt.plot(cdf[N//3:2*N//3-1], label='MEM_lo_fluxes')
plt.plot(cdf[2*N//3:N-1], label='NEO_fluxes')

plt.legend()
plt.show()