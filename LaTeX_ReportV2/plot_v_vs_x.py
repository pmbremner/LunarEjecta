#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
import scipy.special as sc

# https://pundit.pratt.duke.edu/wiki/Python:Finding_roots
import scipy.optimize as opt

g = 1.625 # m/s^-2, lunar gravity
Y = 4.E3 # Pa (Si units), from SFA, use for all
K = 5 # crater depth-to-diameter ratio

# define constants (for SFA)
k     = 0.3
C1    = 0.55
rho   = 1.3    #g/cm^3, target density
delta_i = 4.   #g/cm^3, impactor density
n1    = 1.2
n2    = 1.
p     = 0.3
mu    = 0.4
nu    = 0.4
H1    = 0.59
H2    = 0.4


# # define constants (for WCB, Weakly Cemented Basalt)
# k     = 0.3
# C1    = 0.18
# rho   = 1.3    #g/cm^3, target density
# delta_i = 4.   #g/cm^3, impactor density
# n1    = 1.2
# n2    = 1.
# p     = 0.3
# mu    = 0.46
# nu    = 0.4
# H1    = 0.5
# H2    = 0.38

# define power-law constants
alpha = 0.56 # basalt
beta = 2.
delta_ind = 9. * mu / alpha



def v_func(x, U, delta_dens, R, a):
	return C1 * U * ((x/a) * (rho/delta_dens)**(nu))**(-1./mu) * (1. - x/(n2*R))**(p)

def v_func_0(x, v, U, delta_dens, R, a):
	return v_func(x, U, delta_dens, R, a) - v/U

def v_max(U, delta_dens, R, a):
	return v_func(n1, U, delta_dens, R, a)

def get_x(v, U, delta_dens, R, a):
	vmax_i = v_max(U, delta_dens, R, a)
	if v > vmax_i:
		return np.nan
	else:
		xmin = n1
		xmax = n2*R
		return opt.brentq(lambda x: v_func_0(x, v, U, delta_dens, R, a), xmin, xmax)

get_xvec = vectorize(get_x)


m_imp        = float(sys.argv[1])         # g
U            = float(sys.argv[2]) * 1000. # [U] = km/s, to m/s
delta_dens_i = float(sys.argv[3])         # g/cm^3


N = 1000
M = 10
ar = np.linspace(0.01, 1./n1, M) # n1*a / (n2*R)


for i in range(M):
	xa = np.linspace(n1, n1/ar[i])

	plt.loglog(xa/n1, v_func(xa, 1., delta_dens_i, n1/(ar[i]*n2), 1.), label=r'$\frac{n_1 a}{n_2 R} = $' + f'{ar[i]:.2f}')

plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.title("Lunar Regolith Target (sand fly ash)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens_i:.1f} g/cc")
plt.xlabel(r'Location in Crater $x/a$')
plt.ylabel(r'Ejecta Speed $v/U$')

plt.legend()
#plt.savefig('v_vs_x_hidens.png', bbox_inches='tight', dpi=600)
#plt.savefig('v_vs_x_lodens.png', bbox_inches='tight', dpi=600)

plt.figure()

ar = np.logspace(-3., np.log10(1.), N)

# SFA target, high density impactor
delta_dens = 4.
plt.plot(ar, v_max(1., delta_dens, n1/(ar*n2), 1.), label=r'Target = SFA, $\delta = $' + f'{delta_dens:.1f} g/cc')

# SFA target, low density impactor
delta_dens = 1.
plt.plot(ar, v_max(1., delta_dens, n1/(ar*n2), 1.), label=r'Target = SFA, $\delta = $' + f'{delta_dens:.1f} g/cc')

# define constants (for WCB, Weakly Cemented Basalt)
k     = 0.3
C1    = 0.18
rho   = 1.3    #g/cm^3, target density
#delta_i = 4.   #g/cm^3, impactor density
n1    = 1.2
n2    = 1.
p     = 0.3
mu    = 0.46
nu    = 0.4
H1    = 0.5
H2    = 0.38

# WCB target, high density impactor
delta_dens = 4.
plt.plot(ar, v_max(1., delta_dens, n1/(ar*n2), 1.), label=r'Target = WCB, $\delta = $' + f'{delta_dens:.1f} g/cc')

# WCB target, low density impactor
delta_dens = 1.
plt.plot(ar, v_max(1., delta_dens, n1/(ar*n2), 1.), label=r'Target = WCB, $\delta = $' + f'{delta_dens:.1f} g/cc')

plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.legend(loc='upper right')
plt.xlabel(r"Effective Impactor to Crater Radius $\frac{n_1 a}{n_2 R}$", fontsize=14)
plt.ylabel(r'Maximum Ejected Speed $v_{max}/U$', fontsize=14)

#plt.savefig('vmax_vs_x.png', bbox_inches='tight', dpi=600)

plt.figure()
delta_dens = 4.
ar = np.linspace(0.01, 1./n1, M) # n1*a / (n2*R)

vi = np.logspace(-6, 1., N)
for i in range(M):
	xi = get_xvec(vi, 1., delta_dens, n1/(ar[i]*n2), 1.)

	plt.loglog(vi, xi, label=r'$\frac{n_1 a}{n_2 R} = $' + f'{ar[i]:.2f}')

plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.title("Lunar Regolith Target (sand fly ash)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens:.1f} g/cc")
plt.ylabel(r'Location in Crater $x/a$')
plt.xlabel(r'Ejecta Speed $v/U$')
plt.legend()

plt.show()