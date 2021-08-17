#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import random as rng


def hit_func(x, y):
	hit = 0
	#if np.abs(x - 0.3) < 0.01 and y < 0.5:
	if x**2 + y**2 < 1.:
		hit = 1
	return hit



N_tot_hits = int(sys.argv[1])
N_max_pop  = int(sys.argv[2])
alpha      = float(sys.argv[3])
dr         = float(sys.argv[4])

hit_flag = 0
search_flag = 0
pop_full_flag = 0
current_pop_integral = np.zeros(N_max_pop)
current_pop_count    = np.zeros(N_max_pop)
current_pop_volume   = np.zeros(N_max_pop)

xi = np.zeros(N_max_pop)
yi = np.zeros(N_max_pop)

search_integral = 0.
search_count = 0.

cur_pop_idx = 0

hit_count = 0

weight = 0.
x = 0.
y = 0.
jj = 0

total_integral = 0.

while hit_count < N_tot_hits:
	
	# search
	if(rng.uniform(0, 1.) < alpha or hit_count == 0):
		search_flag = 1
		search_count += 1.

		x = rng.uniform(0, 1.)
		y = rng.uniform(0, 1.)
		weight = 1./alpha

	# destroy
	else:
		search_flag = 0
		jj = rng.randint(0, min(hit_count, N_max_pop)-1)

		current_pop_count[jj] += 1.

		x = rng.uniform(max(0., xi[jj]-dr), min(1., xi[jj]+dr))
		y = rng.uniform(max(0., yi[jj]-dr), min(1., yi[jj]+dr))

		if current_pop_count[jj] == 1:
			current_pop_volume[jj] = (min(1., xi[jj]+dr) - max(0., xi[jj]-dr)) * (min(1., yi[jj]+dr) - max(0., yi[jj]-dr))

		weight = 1./(1.-alpha) * float(min(hit_count, N_max_pop))

	# check hit region
	hit_flag = hit_func(x, y)

	if hit_flag == 1:

		# integrate/count the hit
		if search_flag == 1:
			search_integral += weight
		else:
			current_pop_integral[jj] += weight

		# replace old pop, saving integral
		if pop_full_flag == 1 and current_pop_count[cur_pop_idx] > 0:
			total_integral += current_pop_volume[cur_pop_idx] * current_pop_integral[cur_pop_idx] / float(current_pop_count[cur_pop_idx])
			current_pop_volume[cur_pop_idx] = 0.
			current_pop_integral[cur_pop_idx] = 0.
			current_pop_count[cur_pop_idx] = 0

		# store new hit location
		xi[cur_pop_idx] = x
		yi[cur_pop_idx] = y

		
		cur_pop_idx += 1
		hit_count += 1

		# roll over
		if cur_pop_idx >= N_max_pop:
			pop_full_flag = 1
			cur_pop_idx = 0

# finished, count current pop remaining
for i in range(N_max_pop):
	if current_pop_count[i] > 0:
		total_integral += current_pop_volume[i] * current_pop_integral[i] / float(current_pop_count[i])

# count up search space integral
total_integral += search_integral / float(search_count)

print(total_integral)