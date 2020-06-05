#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Test.txt", unpack=True)

plt.contourf(data)
plt.show()