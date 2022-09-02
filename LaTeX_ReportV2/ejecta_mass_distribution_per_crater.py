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

# # define constants (for SFA)
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


# define constants (for WCB, Weakly Cemented Basalt)
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

# [m_imp]     = g
# [delta_imp] = g/cm^3
# [radius_a]  = cm
def radius_a(m_imp, delta_imp):
	return (3.*m_imp / (3.*np.pi*delta_imp))*(1./3.)

# [delta_imp]  = g/cm^3
# [U]          = km/s
# [R_bar_str]  = n2*R / (n1*a)
def R_bar_str(delta_imp, U):
	return n2*H2/n1 *(rho/delta_imp)**(-nu) * (1.E-9 * Y/(rho*U**2))**(-mu/2.)

# [x_bar]     = a
# [delta_imp] = g/cm^3
# [M_bar]     = m_imp
def M_bar(x_bar, delta_imp):
	return 3.*k/(4.*np.pi) * (rho/delta_imp) * ((x_bar)**3 - n1**3)

# [delta_imp] = g/cm^3
# R_bar       = n2*R / (n1*a)
# [M_bar_tot] = m_imp
def M_bar_tot(delta_imp, R_bar):
	return M_bar(n1 * R_bar, delta_imp)

# [delta_imp] = g/cm^3
# R_bar       = n2*R / (n1*a)
# [m_imp]     = g
# [mb]        = g
def mb(delta_imp, R_bar, m_imp):
	return 0.2 * (M_bar_tot(delta_imp, R_bar) * m_imp)**0.8

# [x_bar]     = a
# [delta_imp] = g/cm^3
# R_bar       = n2*R / (n1*a)
# [v_bar_func] = U
def v_bar_func(x_bar, delta_imp, R_bar):
	return C1 * (x_bar * (rho/delta_imp)**(nu))**(-1./mu) * (1. - x_bar/(n1*R_bar))**(p)

# [delta_imp] = g/cm^3
# R_bar       = n2*R / (n1*a)
# [v_bar_max] = U
def v_bar_max(delta_imp, R_bar):
	return v_bar_func(n1, delta_imp, R_bar)

# R_bar       = n2*R / (n1*a)
# [U]         = km/s
# [a]         = cm
# [v_bar_min] = U
def v_bar_min(R_bar, U, a):
	return 2./(1000.*U) *np.sqrt((g*n1*a*R_bar) / (100.*n2*K))

# [x_bar]     = a
# [delta_imp] = g/cm^3
# R_bar       = n2*R / (n1*a)
# [v_bar]     = U
def v_bar_0(x_bar, v_bar, delta_imp, R_bar):
	return v_bar_func(x_bar, delta_imp, R_bar) - v_bar

# [v_bar]     = U
# [delta_imp] = g/cm^3
# R_bar       = n2*R / (n1*a)
# [x_bar]     = a
def get_x_bar(v_bar, delta_imp, R_bar):
	v_bar_max_i = v_bar_max(delta_imp, R_bar)
	if np.iscomplex(v_bar_max_i) or v_bar > v_bar_max_i:
		return np.nan
	else:
		x_min = n1
		x_max = n1*R_bar
		return opt.bisect(lambda x: v_bar_0(x, v_bar, delta_imp, R_bar), x_min, x_max, xtol=1E-26)#, xtol=1E-20) #brentq, bisect, toms748

get_x_barvec = vectorize(get_x_bar)

# [m_imp]     = g
# [m_bar_ej]  = mb
# [delta_imp] = g/cm^3
# [U]         = km/s
def ejecta_mass_distribution(m_imp, m_bar_ej, delta_imp, U, inc_beta_flag):

	R_bar = R_bar_str(delta_imp, U)    # n2*R / (n1*a)
	a     = radius_a(m_imp, delta_imp) # [cm]

	v_bar_min_i = v_bar_min(R_bar, U, a) # [U]
	v_bar_max_i = np.minimum(v_bar_min_i * m_bar_ej**(-1./delta_ind), v_bar_max(delta_imp, R_bar)) # [U]

	x_bar_min_i = get_x_barvec(v_bar_max_i, delta_imp, R_bar) # [a]
	x_bar_max_i = get_x_barvec(v_bar_min_i, delta_imp, R_bar) # [a]


	mb_i = mb(delta_imp, R_bar, m_imp) # g

	if inc_beta_flag == 0:
		return x_bar_min_i, x_bar_max_i, 3.*k*beta / (4.*np.pi)  * (rho/delta_imp)**(-beta*delta_ind*nu/(3.*mu) + 1.) * (m_imp/mb_i) * (m_bar_ej)**(beta/3. - 1.) * (v_bar_min_i/C1)**(-beta*delta_ind/3.) * (n1*R_bar)**(-beta*delta_ind/(3.*mu))
	else:
		af = -beta*delta_ind/(3.*mu)+3
		bf = beta*delta_ind*p/3. + 1.

		x0 = x_bar_max_i/(n1*R_bar)
		x1 = x_bar_min_i/(n1*R_bar)

		# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.hyp2f1.html
		# FA = (x_bar_max_i/(n1*R_bar))**(af) * sc.hyp2f1(af, 1.-bf, af+1., x_bar_max_i/(n1*R_bar)) / (af)
		# FB = (x_bar_min_i/(n1*R_bar))**(af) * sc.hyp2f1(af, 1.-bf, af+1., x_bar_min_i/(n1*R_bar)) / (af)

		FB = (1.-x0)**bf * x0**af * sc.hyp2f1(1, af+bf, bf+1, 1.-x0) / bf
		FA = (1.-x1)**bf * x1**af * sc.hyp2f1(1, af+bf, bf+1, 1.-x1) / bf

		return x1, x0, FA, FB, mb_i, 3.*k*beta / (4.*np.pi)  * (rho/delta_imp)**(-beta*delta_ind*nu/(3.*mu) + 1.) * (m_imp/mb_i) * (m_bar_ej)**(beta/3. - 1.) * (v_bar_min_i/C1)**(-beta*delta_ind/3.) * (n1*R_bar)**(-beta*delta_ind/(3.*mu)) * (FA - FB) 





N = 1000

m_min_exp = -10.
m_max_exp = 0.

m_bar_ej = np.logspace(m_min_exp, m_max_exp, N)

m_imp        = float(sys.argv[1]) # g
U            = float(sys.argv[2]) # [U] = km/s
delta_dens_i = float(sys.argv[3]) # g/cm^3


fig = plt.figure()
ax = plt.subplot(111)

for m_imp_i in np.logspace(-6., 6., 13):
	#plt.figure()
	xmin, xmax, FA, FB, mb_i, ej_dist = ejecta_mass_distribution(m_imp_i, m_bar_ej, delta_dens_i, U, 1)

	ej_cum = ej_dist[:-1] * (m_bar_ej[1:] - m_bar_ej[:-1]) * mb_i / m_imp_i
	ej_cum = np.cumsum(ej_cum)

	plt.loglog(m_bar_ej[:-1], ej_cum[-1] - ej_cum, label=r'$m_{imp} = $'+f'{m_imp_i:.1e}' + r', $m_b = $' + f'{mb_i:.2e}'+' g')


plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

# https://www.delftstack.com/howto/matplotlib/set-matplotlib-grid-interval/
major_ticks_topx = np.logspace(m_min_exp, m_max_exp, 11)
minor_ticks_topx = np.logspace(m_min_exp, m_max_exp, 21)

major_ticks_topy = np.logspace(-23, -1, 23)
minor_ticks_topy = np.logspace(-23, -1, 23)

ax.set_xticks(major_ticks_topx)
ax.set_yticks(major_ticks_topy)
ax.set_xticks(minor_ticks_topx,minor=True)
ax.set_yticks(minor_ticks_topy,minor=True)

ax.grid(which="major",alpha=0.6)
ax.grid(which="minor",alpha=0.3)

plt.title("Lunar Regolith Target (sand fly ash)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens_i:.1f} g/cc, Impactor Speed {U:.0f} km/s")
#plt.title("Lunar Regolith Target (weakly cemented basalt)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens_i:.1f} g/cc, Impactor Speed {U:.0f} km/s")
plt.xlabel(r'Ejecta Particle Mass $\bar{m}_{ej} = m_{ej}/m_b$', fontsize=14)
plt.ylabel(r'Ejecta Yield for Particle Mass $>\bar{m}_{ej}$', fontsize=14)

# https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('IntegralEjectedMass_vs_ParticleMass.png', bbox_inches='tight', dpi=600)
plt.show()

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

Nm = 100
NU = 100

N = 1000

m_min_exp = -50.
m_max_exp = 0.

m_bar_ej = np.logspace(m_min_exp, m_max_exp, N)

m_ej_user = 1.E-6 # g


m_imp_v = np.logspace(-6., 6., Nm)
U_v     = np.logspace(-2, 2, NU)

Xi, Yj = np.meshgrid(m_imp_v, U_v)
Zij = np.zeros(np.shape(Xi))

for i, m_imp_i in enumerate(m_imp_v): # g
	for j, U_j in enumerate(U_v): # km/s
	
		#print(m_imp_i, m_bar_ej, delta_dens_i, U_j)

		#plt.figure()
		xmin, xmax, FA, FB, mb_i, ej_dist = ejecta_mass_distribution(m_imp_i, m_bar_ej, delta_dens_i, U_j, 1)

		ej_cum = ej_dist[:-1] * (m_bar_ej[1:] - m_bar_ej[:-1]) * mb_i / m_imp_i
		ej_cum = np.cumsum(ej_cum)


		Y_ij = ej_cum[-1] - ej_cum
		x_mask = m_bar_ej[:-1][m_ej_user/mb_i >= m_bar_ej[:-1]]

		if x_mask.size > 0:
			Zij[j][i] = Y_ij[np.argmax(x_mask)]
		else:
			Zij[j][i] = np.nan


plt.figure()

plt.pcolor(np.log10(Xi), np.log10(Yj), Zij, cmap=plt.cm.nipy_spectral)

plt.title("Lunar Regolith Target (sand fly ash)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens_i:.1f} g/cc")
plt.title("Lunar Regolith Target (weakly cemented basalt)\n" r"Regolith Bulk Density $\rho = $" + f"{rho:.1f}" +f" g/cc, Impactor Density {delta_dens_i:.1f} g/cc")
#plt.xlabel(r'Impactor Mass [g] $\log(m_{imp})$', fontsize=14)
plt.ylabel(r'Impactor Speed [km/s] $\log(U)$', fontsize=14)
cbar = plt.colorbar(label='Ejecta Yield for Ejecta Particle Mass > ' + f'{m_ej_user:.1e} g')
cbar.ax.tick_params(labelsize=12)
plt.grid(b=True, which='both', color='grey', linewidth=1.5) # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.savefig('EjectaYield_vs_ImpactorSpeed_vs_Mass.png', bbox_inches='tight', dpi=600)
plt.show()