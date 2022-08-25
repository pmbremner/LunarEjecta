#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
import scipy.special as sc

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
def R_crater(m_imp, U, delta_dens):
	a = radius(m_imp, delta_dens)
	return np.min(H2 * (m_imp / rho)*(1./3.) * (rho / delta_dens)**((1.-3.*nu)/3.) * (Y / (rho * U**2))**(-mu/2.), H1 * (m_imp / rho)*(1./3.) * (rho / delta_dens)**((2.+mu-6.*nu)/(3.*(2.+mu))) * (g*a / U**2)**(-mu/(2.+mu)))

def M_tot(m_imp, U, delta_dens):
	return 3.*k/(4.*np.pi) * (rho/delta_dens) *((n2*R/a)**3 - n1**3) * m_imp


def f_C(m_ej, m_imp, U, delta_dens):
	R = R_crater(m_imp, U, delta_dens) # cm
	a = radius(m_imp, delta_dens) # cm

	Mtot = M_tot(m_imp, U, delta_dens) # g
	mb = 0.01 * Mtot

	vmin = 2. * np.sqrt(g * R / (K * 100.)) # m/s

	vmax0 = U*C1 * (n1 * (rho/delta_dens)**nu)**(-1./mu) * (1. - n1*a/(n2*R))**p # m/s

	vmax1 = vmin * (m_ej/mb)**(-1./delta_ind) # m/s

	vmax_m = np.min(vmax0, vmax1) # m/s

	xmax = x_val(vmax_m)
	xmin = x_val(vmin)


	return 3.*k*beta / (4.*np.pi) * (rho/delta_dens)**(-beta*delta_ind*nu/mu + 1.) * (m_imp/mb) * (m_ej/mb)**(beta/3. - 1.) * (C1*U/vmin)**(beta*delta_ind/3.) * (a / (n2*R))**(beta*delta_ind/(3.*mu) - 3.) * (sc.betainc(-beta*delta_ind/(3.*mu)+3., beta*delta_ind*p/3. + 1, xmax) - sc.betainc(-beta*delta_ind/(3.*mu)+3., beta*delta_ind*p/3. + 1, xmin))