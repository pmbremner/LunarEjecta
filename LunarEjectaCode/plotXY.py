#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('test.txt', unpack=True)

plt.scatter(x, y, s=1)
plt.show()