#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

r = float(sys.argv[1])
N = int(sys.argv[2])

s = np.random.poisson(r, N)

plt.hist(s, ec="k")
plt.show()