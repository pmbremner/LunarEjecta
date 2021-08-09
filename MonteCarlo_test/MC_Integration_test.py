#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

def p_norm(a, xmin, xmax):
	return (xmax**(a+1.) - xmin**(a+1.)) / (a+1.)

def p_law(x, a, xmin, xmax):
	return x**a / p_norm(a, xmin, xmax)

def inv_CDF(r, a, xmin, xmax): # r is multiplied by a normalization factor
	return (xmax**(a+1.) - r)**(1./(a+1.))
	#return ((1.-r)*xmax**(a+1.) + r*xmin**(a+1.))**(1./(a+1.))


N = int(sys.argv[1])
M = int(sys.argv[2])
a = float(sys.argv[3])

xmin = 0.001
xmax = 1.

sum_v0 = np.zeros(M)
sum_v1 = np.zeros(M)
sum_v2 = np.zeros(M)


for j in range(M):

	for i in range(N):

		r = np.random.rand()

		x0 = xmin + (xmax - xmin) * r  # same as inv_CDF(x, 0, xmin, xmax, r)
		r = np.random.rand()
		x1 = inv_CDF(r *p_norm(0, xmin, xmax), 1., xmin, xmax) # p_norm(0, xmin, xmax)
		r = np.random.rand()
		x2 = inv_CDF(r *p_norm(0, xmin, xmax), 0.5, xmin, xmax) # p_norm(0, xmin, xmax)

		sum_v0[j] += p_law(x0, a, xmin, xmax)
		sum_v1[j] += p_law(x1, a, xmin, xmax) / p_law(x1, 1., xmin, xmax) 
		sum_v2[j] += p_law(x2, a, xmin, xmax) / p_law(x2, 0.5, xmin, xmax) 


sum_v0 /= float(N) 
sum_v1 /= float(N) 
sum_v2 /= float(N) 

plt.plot(range(M), sum_v0, label='std')
plt.plot(range(M), sum_v1, label='importance | a = 1.')
plt.plot(range(M), sum_v2, label='importance | a = 0.5')
plt.legend()

plt.figure()
plt.hist(sum_v0, ec="k", alpha=0.5, label='std')
plt.hist(sum_v1, ec="k", alpha=0.5, label='importance | a = 1.')
plt.hist(sum_v2, ec="k", alpha=0.5, label='importance | a = 0.5')
plt.legend()
plt.show()