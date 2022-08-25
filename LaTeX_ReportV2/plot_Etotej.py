#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
import scipy.special as sc

# # define constants (for SFA)
# k     = 0.3
# C1    = 0.55
# rho   = 1.3    # target density
# delta_i = 4.     # impactor density
# n1    = 1.2
# n2    = 1.
# p     = 0.3
# mu    = 0.4
# nu    = 0.4

# define constants (for WCB, Weakly Cemented Basalt)
k     = 0.3
C1    = 0.18
rho   = 1.3    # target density
delta_i = 4.     # impactor density
n1    = 1.2
n2    = 1.
p     = 0.3
mu    = 0.46
nu    = 0.4

def Etotej(x, delta):
	#return 9.*k*C1**2/(8.*np.pi) * (rho/delta)**(-2.*nu/mu+1) * (x/n1)**(2/mu-3.) * sc.betainc(2.*p + 1, -2./mu + 3, 1.-x)
	return 9.*k*C1**2/(8.*np.pi) * (rho/delta)**(-2.*nu/mu+1) * (x/n1)**(2/mu-3.) * (1.-x)**(2.*p+1.) * sc.hyp2f1(2.*p + 1, 2./mu - 2, 2.*p + 2, 1.-x) / (2.*p+1.) / 0.5


x = np.logspace(-5, 0, 1000)


delta_i = np.linspace(0.5, 5., 10)

for d_i in delta_i:
	plt.plot(x, Etotej(x, d_i), label=r"$\delta =$ " + f"{d_i:.1f}")
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.title("Lunar Regolith Target (weakly cemented basalt)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +r" g/cc, Impactor Density $\delta$")
plt.xlabel(r"Effective Impactor to Crater Radius $\frac{n_1 a}{n_2 R}$", fontsize=14)
plt.ylabel(r"Total Ejecta Kinetic Energy $E_{tot,ej} / E_{imp}$", fontsize=14)
plt.legend()
plt.savefig('TotalEjectaKE_vs_EffectiveCraterSize_3.png', bbox_inches='tight', dpi=600)
plt.show()