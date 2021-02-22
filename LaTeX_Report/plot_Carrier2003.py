import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt
import scipy.special as sp

# x: particle diameter in mm
def Carrier2003_MassWeightedCDF(x, mu, sigma):
	return 0.5 * sp.erfc((np.log(x) - mu)/(np.sqrt(2.)*sigma))

def Carrier2003_MassWeightedPDF(x, mu, sigma):
	return 1. / (sigma*np.sqrt(2.*np.pi)*x) * np.exp(-(np.log(x) - mu)**2 / (2.*sigma**2))

# rho: particle density in g/cm^3
def Carrier2003_NumberWeightedCDF(x, mu, sigma ,rho):
	return 6. / (np.pi*rho * 1.E-3 * sigma*np.sqrt(2.*np.pi)) * sp.erfc((np.log(x)-mu+3.*sigma**2)/(np.sqrt(2.)*sigma)) * np.exp(-3.*mu + 9.*sigma**2/2.)

def Carrier2003_NumberWeightedPDF(x, mu, sigma ,rho):
	return 6. / (np.pi*rho * 1.E-3 * sigma*np.sqrt(2.*np.pi) * x**4) * np.exp(-(np.log(x) - mu)**2 / (2.*sigma**2))


# read in digitized data from Carrier 2003 Fig. 1
C03_size_mm, C03_percentByMass = np.loadtxt('Carrier2003_digitized.txt', unpack=True)

C03_mu    = -2.649
C03_sigma = 1.786
DSNE_rho  = 3.1 # g/cc

x_mm = np.logspace(-3, 1, 200)
y_massWeightedCDF = Carrier2003_MassWeightedCDF(x_mm, C03_mu, C03_sigma)

fig = plt.subplot(1)
plt.semilogx(C03_size_mm, 100. - C03_percentByMass, '.', label='Carrier 2003 Fig. 1 digitized')
plt.semilogx(x_mm, y_massWeightedCDF*100., label=r'Log-normal fit using $\mu=$'+str(C03_mu)+r' and $\sigma=$'+str(C03_sigma))
plt.xlabel('Particle Diameter (mm)')
plt.ylabel('Percent (mass-weighted) > x')
plt.title('Mass-Weighted Particle CDF')
plt.grid(True, which='both')
plt.legend()
fig.savefig('Carrier2003_MassWeightedCDF.png', dpi=600)

###########################################################################################
y_massWeightedPDF = Carrier2003_MassWeightedPDF(x_mm, C03_mu, C03_sigma)

fig = plt.figure(2)
plt.semilogx(x_mm[:-1], y_massWeightedPDF[:-1] * (x_mm[1:]-x_mm[:-1]))
plt.xlabel('Particle Diameter (mm)')
plt.ylabel('Fractional Percent (PDF * dx)')
plt.title('Mass-Weighted Particle PDF * dx')
plt.grid(True, which='both')
fig.savefig('Carrier2003_MassWeightedPDF.png', dpi=600)

###########################################################################################
y_numberWeightedCDF = Carrier2003_NumberWeightedCDF(x_mm, C03_mu, C03_sigma, DSNE_rho)

fig = plt.figure(3)
plt.loglog(x_mm, y_numberWeightedCDF)
plt.xlabel('Particle Diameter (mm)')
plt.ylabel('Number of Particles in 1 kg > x')
plt.title('Number-Weighted Particle CDF')
plt.grid(True, which='both')
fig.savefig('Carrier2003_NumberWeightedCDF.png', dpi=600)

###########################################################################################
y_numberWeightedPDF = Carrier2003_NumberWeightedPDF(x_mm, C03_mu, C03_sigma, DSNE_rho)

#x_mm = np.logspace(-10,1,200)

fig = plt.figure(4)
plt.loglog(x_mm, y_numberWeightedCDF)
plt.xlabel('Particle Diameter (mm)')
plt.ylabel('Number of Particles in 1 kg (1/mm)')
plt.title('Number-Weighted Particle PDF')
plt.grid(True, which='both')
#plt.minorticks_on()
fig.savefig('Carrier2003_NumberWeightedPDF.png', dpi=600)

plt.show()