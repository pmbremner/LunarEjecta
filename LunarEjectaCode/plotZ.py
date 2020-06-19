#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("meshgrid_2_D0D1_0.024797_0.025360.txt", unpack=True)

plt.contourf(data,1000)
plt.show()