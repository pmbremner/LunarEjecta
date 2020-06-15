#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("run0.txt", unpack=True)

print(np.shape(data)[0]/20)

data2 = data.reshape(int(np.shape(data)[0]/20), 20)

np.savetxt("run0_edit.txt", data2)