#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

iso_filename = 'vmin_IntFraction_vs_vmin.txt'
a45_filename = 'a-45_IntFraction_vs_vmin.txt'

iso_vmin, iso_flux = np.loadtxt(iso_filename, unpack=True)
a45_vmin, a45_flux = np.loadtxt(a45_filename, unpack=True)

plt.loglog(iso_vmin, iso_flux, label='Isotropic angular distribution', linestyle='-', marker='o')
plt.loglog(a45_vmin, a45_flux, label='45 degree-centered angular distribution', linestyle='-', marker='o')
plt.legend()
plt.xlabel('Minimum ejecta speed (m/s)')
plt.ylabel('Fraction of total ejecta')
plt.grid(b=True, which='both')
plt.savefig('IntFraction_vs_vmin_comparison.png', dpi=1000, bbox_inches='tight', pad_inches=0.1)
plt.show()