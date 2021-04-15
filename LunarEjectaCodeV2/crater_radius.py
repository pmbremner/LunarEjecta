import numpy as np
import matplotlib.pyplot as plt
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement

rho   = 3100. * 0.45 # kg/m^3, regolith density times porosity
delta = 2000. # kg/m^3, target density

mu = 0.4
nu = 0.4

H1 = 0.8
H2 = 0.4

g = 1.625 # m/s^2
Y = 2.E+3 # pa

def target_a(m, rho1):
	return (3.*m/(4.*np.pi*rho1))**(1./3.)

def R_grav(m, U):
	return H1 * (m/rho)**(1./3.) * (rho/delta)**((2.+mu-6.*nu)/(3.*(2.+mu))) * (g*target_a(m,rho)/U**2)**(-mu/(2.+mu))

def R_strg(m, U):
	return H2 * (m/rho)**(1./3.) * (rho/delta)**((1.-3.*nu)/3.) * (Y/(rho*U**2))**(-mu/2.)

vR_grav = vectorize(R_grav)
vR_strg = vectorize(R_strg)


R0 = 0.5
R1 = 1.5

def R_both(m, U):
	R = rho*g*target_a(m,rho)/Y

	if R < R0:
		return R_strg(m, U)
	elif R < R1:
		return ((R - R0) * R_grav(m,U) + (R1 - R) * R_strg(m, U)) / (R1 - R0)
	else:
		return R_grav(m,U)


vR_both = vectorize(R_both)

U1 = 1000.

m = np.logspace(-6, 13, 1000)

#print(target_a(m,rho))

U = np.linspace(100., 40000., 1000)

plt.loglog(m, vR_grav(m, U1), label='grav')
plt.loglog(m, vR_strg(m, U1), label='strg')
plt.loglog(m, vR_both(m, U1), label='both')
plt.legend()
plt.show()