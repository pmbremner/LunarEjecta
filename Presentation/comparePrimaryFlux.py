#!/usr/bin/python

import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt

def primary_sp8013(m_gram):
	if m_gram < 1E-6:
		return 10.**(-14.339 - 1.584*np.log10(m_gram) - 0.063*(np.log10(m_gram))**2)
	else:
		return 10.**(-14.37 - 1.213*np.log10(m_gram))

vprimary_sp8013 = vectorize(primary_sp8013)

def primary_grun85(m_gram):
	return (2.2E3*m_gram**0.306 + 15.)**-4.38 + 1.3E-9*(m_gram + 1.E11*m_gram**2 + 1.E27*m_gram**4)**-0.36 + 1.3E-16*(m_gram + 1.E6*m_gram**2)**-0.85

def primary_zook1967(m_gram):
	if m_gram < 1E-8:
		return 10.**-8 * m_gram**-0.3
	elif m_gram < 10.**-2.15:
		return 10.**-13.6 * m_gram**-1
	else:
		return 10.**-14.33 * m_gram**-1.34

vprimary_zook1967 = vectorize(primary_zook1967)

m = np.logspace(-12, 0, 1000)

plt.loglog(m, vprimary_sp8013(m), label='NASA SP-8013')
plt.loglog(m, primary_grun85(m), label=r'Gr$\"{u}$n et al. 1985')
plt.loglog(m, vprimary_zook1967(m), label='DS-21 Rev A, see Zook 1967')
plt.title('Comparison of Primary Flux-Mass')
plt.ylabel(r'Particles of mass $>m$ [#/m$^2$/s]')
plt.xlabel(r'Mass $m$ [g]')
plt.legend()
plt.grid()
plt.savefig('comparePrimaryFlux.png',dpi=600)
plt.show()