import numpy as np
import matplotlib.pyplot as plt

N_phi   = 36
N_theta = 72 
N_v     = 40

d_phi   = 5
d_theta = 5
d_v     = 2

phi_min   = -90
theta_min = 0
v_min     = 1

data = np.loadtxt('RunData/SouthPole/HiDensity/flux_avg.txt', unpack=True) # Equator South45 SouthPole

v = np.linspace(v_min, v_min + d_v * (N_v-1), N_v)
#print(v)
# data2 = np.reshape(data[:,int(N_phi*N_theta/2):], (N_v+2, N_theta, int(N_phi/2)))
# print(np.shape(data2))
# print(np.shape(np.sum(data2, axis=1)))
# print(data2[:,:,0])
# np.savetxt("test.txt", data2[:,:,0])

# print(data2[0,0,:])
# plt.plot(v, np.sum(data2[2:,:,17], axis=1))
# plt.show()

# plt.plot(v, data2[2:,1,0])
# plt.show()

for i in range(int(N_phi/2), N_phi):
	plt.plot(v, np.sum(data[2:,int(i*N_theta):int((i+1)*N_theta)], axis=1), label=str(data[0, int(i*N_theta)]))

#plt.yscale('log')
plt.legend()
plt.show()