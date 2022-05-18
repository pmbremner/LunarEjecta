#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

def gp(d, vp, pm):
	return np.arctan2(1., vp**2 / np.tan(d/2.) + pm*np.sqrt(vp**4 / np.tan(d/2.)**2 + 2.*vp**2 - 1.))


d  = float(sys.argv[1])

v = np.linspace(0., 1., 1000)

g_p = gp(d, v, 1.)
g_m = gp(d, v, -1.)


plt.plot(g_p, v, label='plus')
plt.plot(g_m, v, label='minus')

plt.legend()
plt.xlim((0., np.pi/2.))
plt.show()
