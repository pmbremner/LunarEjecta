#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import random as rng


N = int(sys.argv[1])
k = float(sys.argv[2])
t = int(sys.argv[3])

sum_i = np.zeros(N)
sum_r_i = np.zeros(N)
sum_r2_i = np.zeros(N)

sum_avg = np.zeros(N)

hit = 0
miss = 0

for i in range(N):
	x = rng.uniform(0, 1.)
	y = rng.uniform(0, 1.)

	if x**2 + y**2 <= 1.:
		hit += 1
	else:
		miss += 1

	t = rng.randint(2, i//2+2)

	sum_i[i] = float(hit) / float(miss+hit)
	sum_r_i[i] = (t**k * sum_i[i] - sum_i[i//t]) / (t**k - 1.)
	sum_r2_i[i] = (t**(2*k) * sum_r_i[i] - sum_r_i[i//t]) / (t**(2*k) - 1.)

	sum_avg[i] = (2.*sum_i[i] + sum_r_i[i] + sum_r2_i[i]) / 4.

print('MC sum = ', 4.*sum_i[-1])
print('RMC sum = ', 4.*sum_r_i[-1])
print('RMC2 sum = ', 4.*sum_r2_i[-1])
print('Avg sum = ', 4.*sum_avg[-1])
plt.loglog(np.abs(4.*sum_r2_i/np.pi-1.), label='Richardson 2 MC')
plt.loglog(np.abs(4.*sum_r_i/np.pi-1.), label='Richardson MC')
plt.loglog(np.abs(4.*sum_avg/np.pi-1.), label='Avg MC')
plt.loglog(np.abs(4.*sum_i/np.pi-1.), label='Normal MC')

plt.legend()

plt.figure()


plt.hist(4.*sum_r2_i, bins=100, range=(3.1,3.2))
plt.hist(4.*sum_r_i, bins=100, range=(3.1,3.2))
plt.hist(4.*sum_avg, bins=100, range=(3.1,3.2))
plt.hist(4.*sum_i, bins=100, range=(3.1,3.2))


plt.show()