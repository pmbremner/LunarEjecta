#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

mx, my, mgen = np.loadtxt("phase_miss_file.txt", unpack=True)
hx, hy, hgen = np.loadtxt("phase_hit_file.txt", unpack=True)

plt.style.use('dark_background')

plt.scatter(mx, my, c=max(mgen) - mgen, cmap=plt.cm.Blues, s=5)
plt.scatter(hx, hy, c=max(hgen) - hgen, cmap=plt.cm.Reds, s=5)
plt.show()