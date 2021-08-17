#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import random as rng


# N     = int(sys.argv[1])
# alpha = float(sys.argv[2])
# dr    = float(sys.argv[3])
# M     = int(sys.argv[4])

Nai = int(sys.argv[1])
Ndi = int(sys.argv[2])
N   = int(sys.argv[3])
M   = int(sys.argv[4])


amin = 0.05
amax = 1.0

dmax = 0.2
dmin = 0.01

var_vec = np.zeros((Nai, Ndi))
count_vec = np.zeros((Nai, Ndi))

alpha = np.linspace(amin, amax, Nai)
dr = np.logspace(np.log10(dmin), np.log10(dmax), Ndi)

for ai in range(Nai):
	for di in range(Ndi):
		
		ans = np.zeros(M)
		c_ans = np.zeros(M)
		c_miss = np.zeros(M)

		for j in range(M):


			xi = np.zeros(N)
			yi = np.zeros(N)


			w = 0.
			i = 0

			hit = 0
			miss = 0

			h_count = 0
			m_count = 0

			f_sum = 0.
			f_count = 0
			i_count = np.zeros(N)
			i_sum = np.zeros(N)

			while i < N:
				
				if(rng.uniform(0, 1.) < alpha[ai] or i == 0):
					f_count += 1

					x = rng.uniform(0, 1.)
					y = rng.uniform(0, 1.)
					w = 1./alpha[ai]
				else:
					jj = rng.randint(0, i-1)

					i_count[jj] += 1
					
					x = rng.uniform(max(0., xi[jj]-dr[di]), min(1., xi[jj]+dr[di]))
					y = rng.uniform(max(0., yi[jj]-dr[di]), min(1., yi[jj]+dr[di]))
					w = 1./(1.-alpha[ai]) * (min(1., xi[jj]+dr[di]) - max(0., xi[jj]-dr[di])) * (min(1., yi[jj]+dr[di]) - max(0., yi[jj]-dr[di]))

				#if x**2 + y**2 <= 1.:
				if np.abs(x - 0.3) < 0.01 and y < 0.5:
					xi[i] = x
					yi[i] = y

					

					hit += w 
					h_count += 1
					i += 1
				else:
					miss += w
					m_count += 1

			ans[j] = hit / (h_count + m_count)
			c_ans[j] = h_count
			c_miss[j] = m_count
		
		#var_vec[ai, di] = np.var(ans)
		var_vec[ai, di] = np.abs(1. - np.average(ans)/0.01)
		count_vec[ai, di] = np.average(c_ans/(c_ans + c_miss)) * 100.
		print('alpha = ', alpha[ai], ' dr = ', dr[di], 'percent done = ',f'{float(ai*Ndi + di)*100./float(Ndi*Nai):.4f}' , 'hit perc = ', count_vec[ai, di], ' err = ', var_vec[ai, di], ' ans = ', np.average(ans), ' total iters = ', m_count+h_count)


plt.pcolor(var_vec)
plt.colorbar()


d = np.linspace(dmin, dmax, Ndi)
a = np.linspace(amin, amax, Nai)

plt.figure()
plt.plot(d, np.average(count_vec, axis=0))
plt.title('dr')

plt.figure()
plt.plot(a, np.average(count_vec, axis=1))
plt.title('alpha')
#plt.legend()

plt.show()