#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://matplotlib.org/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
#from matplotlib import ticker, cm

def beta(a, m, n):
	return 9. * m * n / ( a*(3.*m - 2.*n) + 6.*n)

def gamma(a, m, n):
	return -1. + a/3.*(1. - 3.*m / (2.*n))

def beta2(a, g, n):
	return -n/g * (3./a * (1.+g) - 1.)

def gamma2(a, b, n):
	return -(3. - a) / (3. - a*b/n)




filename = "Zook1967_Fig2.txt" #sys.argv[1]

Ve, G = np.loadtxt(filename, unpack=True)

# truncate first couple points
Ve = Ve[5:]
G = G[5:]

# compute approx of derivative
dVe = (Ve[1:] - Ve[:-1])
dG  = (G[1:] - G[:-1])
dG_dVe = dG / dVe

# compute approx of derivative, using center approx, capping at ends with one-sided
dG_dVe_c = 1*dG_dVe # need the 1*, otherwise a soft copy will be made
dG_dVe_c[1:] += dG_dVe_c[:-1]
dG_dVe_c /= 2.
dG_dVe_c = np.append(dG_dVe_c, [[dG_dVe[-1]]])

# If we force the Zook 1967 model to match HH11, then
N = 10
a = np.linspace(0.3, 0.7, N).reshape(1, N)


beta_Ve  = beta(a, 0.4, np.repeat(dG_dVe_c.reshape(np.size(dG_dVe_c), 1), N, axis=1)/2.)
gamma_Ve = gamma2(a, np.average(beta_Ve, axis=0), np.repeat(dG_dVe_c.reshape(np.size(dG_dVe_c), 1), N, axis=1)/2.)


beta_Ve2  = beta2(a, -1.25, np.repeat(dG_dVe_c.reshape(np.size(dG_dVe_c), 1), N, axis=1)/2.)

#print(np.shape(Ve.reshape(1, np.size(Ve) )))

plt.plot(Ve, G)

plt.figure()
plt.plot(Ve[:-1], dG_dVe)
plt.plot(Ve, dG_dVe_c)

print(np.average(dG_dVe_c[:-1], weights=dVe))

plt.figure()
#plt.plot(Ve, beta_Ve)
plt.plot(Ve, beta_Ve2)
#plt.plot(Ve, -gamma_Ve)

plt.figure()
#plt.plot(Ve, gamma_Ve )
plt.plot(a.T, np.average(beta_Ve, axis=0) )
plt.plot(a.T, np.average(beta_Ve2, axis=0) )
#plt.plot(a.T, np.average(-gamma_Ve, axis=0) )
plt.show()