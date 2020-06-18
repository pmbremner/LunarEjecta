#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("meshgrid_3_D0D1_0.007485_0.008500.txt", unpack=True)

plt.contourf(data,1000)
plt.show()