#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import random as rng


def hit_func(x, y):
	hit = 0
	if np.abs(x - 0.3) < 0.01 and y < 0.5: #
	#if x**2 + y**2 < 1.:
	#if y >= x**2 and y <= np.sqrt(x):
		hit = 1
	return hit

def get_dens(xi, yi, x, y, dr, prob_dens, pop_full_flag, cur_pop_idx, N_max_pop):
	N = 0
	dens = 0.

	if pop_full_flag == 1:
		N = N_max_pop
	else:
		N = cur_pop_idx

	if N >= 0:
		for i in range(0, N):
			if np.abs(x - xi[i]) <= dr and np.abs(y - yi[i]) <= dr:
				dens += prob_dens[i]
		dens /= float(N+1)

	return dens

N_tot_hits = int(sys.argv[1])
N_max_pop  = int(sys.argv[2])
alpha      = float(sys.argv[3])
dr         = float(sys.argv[4])

hit_flag = 0
search_flag = 0
pop_full_flag = 0

cur_pop_idx = 0
prob_dens = np.zeros(N_max_pop)
xi = np.zeros(N_max_pop)
yi = np.zeros(N_max_pop)

weight = 0.
x = 0.
y = 0.
jj = 0

hit_integral = 0.
total_integral = 0.
hit_count = 0
attempt_count = 0
search_count = 0
destroy_count = 0

while hit_count < N_tot_hits:
	
	# search
	if(rng.uniform(0, 1.) < alpha or hit_count == 0):
		search_flag = 1
		search_count += 1

		x = rng.uniform(0, 1.)
		y = rng.uniform(0, 1.)

		if hit_count > 0:
			weight = alpha
			if alpha < 1:
				weight += get_dens(xi, yi, x, y, dr, prob_dens, pop_full_flag, cur_pop_idx, N_max_pop) * (1. - alpha)
			weight = 1./weight
		else:
			weight = 1.

	# destroy
	else:
		search_flag = 0
		destroy_count += 1
		jj = rng.randint(0, min(hit_count, N_max_pop)-1)

		x = rng.uniform(max(0., xi[jj]-dr), min(1., xi[jj]+dr))
		y = rng.uniform(max(0., yi[jj]-dr), min(1., yi[jj]+dr))

		weight = 1./( alpha + get_dens(xi, yi, x, y, dr, prob_dens, pop_full_flag, cur_pop_idx, N_max_pop) * (1. - alpha))


	# check hit region
	hit_flag = hit_func(x, y)
	attempt_count += 1
	total_integral += weight

	
	#print('weight = ', weight, ' hit = ', hit_flag, x, y)

	if hit_flag == 1:

		if search_flag == 1:
			plt.scatter(x,y, c='purple')
		else:
			plt.scatter(x,y, c='red')

		hit_integral += weight

		# store new hit location
		xi[cur_pop_idx] = x
		yi[cur_pop_idx] = y

		prob_dens[cur_pop_idx] = 1. / ((min(1., xi[jj]+dr) - max(0., xi[jj]-dr)) * (min(1., yi[jj]+dr) - max(0., yi[jj]-dr)))

		
		cur_pop_idx += 1
		hit_count += 1
		#print('hit_count = ', hit_count)

		# roll over
		if cur_pop_idx >= N_max_pop:
			pop_full_flag = 1
			cur_pop_idx = 0
	else:
		if search_flag == 1:
			plt.scatter(x,y, c='blue')
		else:
			plt.scatter(x,y, c='orange')

total_integral /= float(attempt_count)
hit_integral /= float(attempt_count)

print('integral = ', hit_integral, attempt_count)

plt.show()