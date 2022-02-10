#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


idx, idx_a, idx_b, idx_g, idx_v, f, azm, zenith, speed = np.loadtxt('idx.txt', unpack=True)

plt.hist2d(zenith, speed, 100)
plt.xlim([0, np.pi])
plt.ylim([0, 5000.])
plt.show()