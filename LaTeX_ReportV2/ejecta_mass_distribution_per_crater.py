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

# cm
def radius(m_imp, delta_dens): 
	return (3.*m_imp / (4.*np.pi*delta_dens))**(1./3.)

# cm
def R_crater(m_imp, U, delta_dens, a):
	return np.minimum(H2 * (m_imp / rho)*(1./3.) * (rho / delta_dens)**((1.-3.*nu)/3.) * (Y / (rho * U**2))**(-mu/2.), H1 * (m_imp / rho)*(1./3.) * (rho / delta_dens)**((2.+mu-6.*nu)/(3.*(2.+mu))) * (g*a / U**2)**(-mu/(2.+mu)))

def M_tot(m_imp, U, delta_dens, R, a):
	return 3.*k/(4.*np.pi) * (rho/delta_dens) *((n2*R/a)**3 - n1**3) * m_imp

def x_func(x, v, m_imp, U, delta_dens, R, a):
	return C1 * U * ((x/a) * (rho/delta_dens)**(nu))**(-1./mu) * (1. - x/(n2*R))**(p) - v

# v in m/s
def x_val(v, m_imp, U, delta_dens, R, a):
	xmin = n1 * a
	xmax = n2 * R

	return opt.secant(lambda x: x_func(x, v, m_imp, U, delta_dens, R, a), xmin, xmax)


def f_C(m_ej, m_imp, U, delta_dens):
	a = radius(m_imp, delta_dens) # cm
	R = R_crater(m_imp, U, delta_dens, a) # cm

	Mtot = M_tot(m_imp, U, delta_dens, R, a) # g
	mb = 0.01 * Mtot

	vmin = 2. * np.sqrt(g * R / (K * 100.)) # m/s

	vmax0 = U*C1 * (n1 * (rho/delta_dens)**nu)**(-1./mu) * (1. - n1*a/(n2*R))**p # m/s

	vmax1 = vmin * (m_ej/mb)**(-1./delta_ind) # m/s

	vmax_m = np.minimum(vmax0, vmax1) # m/s

	xmax = x_val(vmax_m, m_imp, U, delta_dens, R, a)
	xmin = x_val(vmin, m_imp, U, delta_dens, R, a)


	return 3.*k*beta / (4.*np.pi) * (rho/delta_dens)**(-beta*delta_ind*nu/mu + 1.) * (m_imp/mb) * (m_ej/mb)**(beta/3. - 1.) * (C1*U/vmin)**(beta*delta_ind/3.) * (a / (n2*R))**(beta*delta_ind/(3.*mu) - 3.) * (sc.betainc(-beta*delta_ind/(3.*mu)+3., beta*delta_ind*p/3. + 1, xmax) - sc.betainc(-beta*delta_ind/(3.*mu)+3., beta*delta_ind*p/3. + 1, xmin))




m_ej_i = np.logspace(-10,2, 1000)


m_imp        = float(sys.argv[1])         # g
U            = float(sys.argv[2]) * 1000. # [U] = km/s, to m/s
delta_dens_i = float(sys.argv[3])         # g/cm^3


#f_C_i = f_C(m_ej_i, m_imp, U, delta_dens_i)


# plt.loglog(m_ej_i, f_C_i)
# plt.show()

N = 1000
M = 10
ar = np.linspace(0.01, 1./n1, M)


for i in range(M):
	xa = np.linspace(n1, n1/ar[i])

	plt.loglog(xa/n1, x_func(xa, 0., m_imp, 1., delta_dens_i, n1/(ar[i]*n2), 1.), label=r'$\frac{n_1 a}{n_2 R} = $' + f'{ar[i]:.2f}')

plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.title("Lunar Regolith Target (sand fly ash)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens_i:.1f} g/cc")
plt.xlabel(r'Location in Crater $x/a$')
plt.ylabel(r'Ejecta Speed $v/U$')

plt.legend()
plt.show()