#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

zenith, speed, weight = np.loadtxt("samples.txt", unpack=True)

plt.scatter(zenith, speed, s=2)
plt.show()