#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

EPS = 1.E-4

# units of vesc
v_domain_max = 1.
v_domain_min = 100. / 2376.

mu = 0.4


def vp_disc(d, g, rs):
	return ((1./rs - np.cos(d))/(1. - np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.)

def vp(d, g, rs):
	vp_disc_i = vp_disc(d, g, rs)
	if vp_disc_i <= 0:
		return v_domain_max
	else:
		return vp_disc_i**-(1/2)


def vmax(d, g, rs, h, a):
	return np.minimum(
		       np.max(np.array([vp(d,        g, rs),
			                    vp(d + 2.*a, g, rs),
			                    vp(d,        g, rs + h),
			                    vp(d + 2.*a, g, rs + h),
			                    v_domain_min])),
		       v_domain_max)



def vmin(d, g, rs, h, a):
	return np.min(np.array([vp(d,        g, rs),
		                    vp(d + 2.*a, g, rs),
		                    vp(d,        g, rs + h),
		                    vp(d + 2.*a, g, rs + h),
		                    v_domain_max]))

vmaxvec = vectorize(vmax)
vminvec = vectorize(vmin)


def dM(d, g, rs, h, a):
	return (vminvec(d, g, rs, h, a)**(-3.*mu) - vmaxvec(d, g, rs, h, a)**(-3.*mu)) * np.sin(g) * 2.*a#* np.sin(d) #/(np.pi*(a+h)) #


def int_dM(d, g, rs, h, a):

	int_dM_i = np.zeros(len(d))

	for i in range(len(d)):
		dM_i = dM(d[i], g, rs, h, a)
		int_dM_i[i] = np.sum( (g[1:] - g[:-1]) * (dM_i[:-1] + dM_i[1:]) / 2.)

	return int_dM_i 


def iint_dM_h(d, g, rs, h, a):

	iint_dM_i = np.zeros(len(h))

	for i in range(len(h)):
		int_dM_i = int_dM(d, g, rs, h[i], a)
		iint_dM_i[i] = np.sum((d[1:] - d[:-1]) * (int_dM_i[:-1] + int_dM_i[1:]) / 2.)

	return iint_dM_i


def iint_dM_r(d, g, rs, h, a):

	iint_dM_i = np.zeros(len(rs))

	for i in range(len(rs)):
		int_dM_i = int_dM(d, g, rs[i], h, a)
		iint_dM_i[i] = np.sum((d[1:] - d[:-1]) * (int_dM_i[:-1] + int_dM_i[1:]) / 2.)

	return iint_dM_i

def iint_dM_a(d, g, rs, h, a):

	iint_dM_i = np.zeros(len(a))

	for i in range(len(a)):
		int_dM_i = int_dM(d, g, rs, h, a[i])
		iint_dM_i[i] = np.sum((d[1:] - d[:-1]) * (int_dM_i[:-1] + int_dM_i[1:]) / 2.)

	return iint_dM_i


rr = float(sys.argv[1])
af = float(sys.argv[2])
hh = float(sys.argv[3])
#dd = float(sys.argv[4])

Ng = 250
Nd = 100
Nh = 20
Na = 10


gp = np.linspace(EPS, np.pi/2.-EPS, Ng)

dd = np.logspace(-4., np.log10(2.*np.pi - 2.*af), Nd)

rv = np.linspace(0., 1., 10) + rr
#aa = np.linspace(0., af, Na)

print(rv-1)

# plt.plot(gp, vmaxvec(dd, gp, rr, hh, aa))
# plt.plot(gp, vminvec(dd, gp, rr, hh, aa))

# plt.loglog(dd, int_dM(dd, gp, rr, 0., aa))
# plt.loglog(dd, int_dM(dd, gp, rr, 2.*hh, aa))
# plt.loglog(dd, int_dM(dd, gp, rr, 4.*hh, aa))
# plt.loglog(dd, int_dM(dd, gp, rr, 8.*hh, aa))

for ri in rv:
	plt.loglog(dd, int_dM(dd, gp, ri, hh, af))

plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

plt.figure()
plt.plot(rv-1, iint_dM_r(dd, gp, rv, hh, af))

# plt.loglog(rr-1., iint_dM_r(dd, gp, rr, hh, aa))
# plt.loglog(rr-1., iint_dM_r(dd, gp, rr, hh, 2.*aa))
# plt.loglog(rr-1., iint_dM_r(dd, gp, rr, hh, 4.*aa))

#plt.plot(aa, iint_dM_a(dd, gp, rr, hh, aa))


plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

#plt.legend()

plt.show()


