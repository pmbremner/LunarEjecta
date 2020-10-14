#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def secEjt_v0(m_gram):
	return 10.**(-10.79-1.2*np.log10(m_gram))

def secEjt_v1(m_gram):
	return 10.**(-11.88-1.2*np.log10(m_gram))

def secEjt_v2(m_gram):
	return 10.**(-13.41-1.2*np.log10(m_gram))


m = np.logspace(-8, 2, 1000)

plt.loglog(m, secEjt_v0(m), label=r'$v_{ej}$ 0 to 0.1 km/s')
plt.loglog(m, secEjt_v1(m), label=r'$v_{ej}$ 0.1 to 0.25 km/s')
plt.loglog(m, secEjt_v2(m), label=r'$v_{ej}$ 0.25 to 1.0 km/s')
plt.title('Secondary Flux-Mass from NASA SP-8013')
plt.ylabel(r'Particles of mass $>m$ [#/m$^2$/s]')
plt.xlabel(r'Mass $m$ [g]')
plt.legend()
plt.grid()
plt.savefig('secFlux_SP8013.png',dpi=600)
plt.show()