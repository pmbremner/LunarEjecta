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

Nazm = 5

def g_p_ap(v, rs):
	return np.arccos(np.sqrt((rs - 1.) * (rs / v**2 - (rs - 1.))))

def F(v, g, rs):
	sign = 1.
	if g > g_p_ap(v, rs):
		sign = -1.

	return sign * np.sqrt(1. + (rs - 1.) * (rs + 1. - rs / v**2) / np.cos(g)**2)

Fvec = vectorize(F, otypes=[float])


def r_final(d, v, g):
	return (2. * (v * np.sin(g))**2) / (1. + (v**2 - 1.) * np.cos(d) - v**2 * np.cos(d - 2.*g))


def F_and_r(v, g, rs, h, d0): 
	sign = 1.
	r_cur = rs
	r_at_d = r_final(d0, v, g)

	if r_at_d > rs + h: # hit top
		sign = -1. 
		r_cur = rs + h
	elif r_at_d < rs: # hit bottom
		sign = 1.
		r_cur = rs
	else: # hit side
		if g > g_p_ap(v, r_at_d): # above local horizon
		 	sign = -1.
		 	r_cur = r_at_d
		else: # below local horizon
		 	sign = 1. 
		 	r_cur = r_at_d

	return sign * np.sqrt(1. + (r_cur - 1.) * (r_cur + 1. - r_cur / v**2) / np.cos(g)**2), r_cur

F_and_rvec = vectorize(F_and_r)


def dist_new(v, g, rs, h, d0):
	Fi, r_cur = F_and_rvec(v, g, rs, h, d0)

	return np.mod(2.*np.pi + 2.*np.arctan2(v**2 * np.sin(g) * np.cos(g) * (1. - Fi/r_cur),  (1. - (v * np.sin(g))**2 * (1. + 1./r_cur) )), 2.*np.pi)



def dist(v, g, rs, s=0):
	Fi = Fvec(v, g, rs)
	if s == -1:
		Fi = -np.abs(Fi)
	elif s == 1:
		Fi = np.abs(Fi)

	if rs == 1.:
		return 2. * np.arctan2(v**2 * np.sin(2.*g), 1. - 2.* (v * np.sin(g))**2 )
	else:
		#return 2. * np.arctan2(v**2 * np.sin(2.*g) * (1. - Fi) + (rs - 1.) * (2.*v**2 - 1.) * np.tan(g), (rs - Fi) - 2.* (v * np.sin(g))**2 * (1. - Fi))
		return np.mod(2.*np.pi + 2.*np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) , 2.*np.pi)
		#return np.mod(2.*np.pi + 2.*np.arctan2(v**2 * np.sin(g)*np.cos(g) * (1. - Fi/rs), 1. - (v*np.sin(g))**2*(1.+1./rs)), 2.*np.pi)

	# Fi = F(v, g, rs)
	# if g > g_ap(d0, rs) and d0 < np.pi:
	# 	Fi *= -1
	# 	if v**2 > 1./(1. + 1./rs):
	# 		return 2.*np.mod(np.pi + np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) , 2.*np.pi)
	# 	else:
	# 		return 2.*np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) 
	# return 2.*np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) 






def vp_disc(d, g, rs):
	return ((1./rs - np.cos(d))/(1. - np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.)

def vp(d, g, rs):
	vp_disc_i = vp_disc(d, g, rs)
	if vp_disc_i <= 0:
		return v_domain_max
	else:
		return vp_disc_i**-(1/2)

vp_vec = vectorize(vp)


def vmax(d, g, rs, h, a):
	return np.minimum(
		       np.max(np.array([vp(d,        g, rs),
			                    vp(d + 2.*a, g, rs),
			                    vp(d,        g, rs + h),
			                    vp(d + 2.*a, g, rs + h),
			                    v_domain_min])),
		       v_domain_max)


def vmaxN(d, g, rs, h, a, N):

	ar = np.concatenate([vp_vec(np.linspace(d, d + 2.*a, N+2, endpoint=True), g, rs),    # bottom
			                    vp_vec(d, g, np.linspace(rs, rs + h, N+2, endpoint=True)),      # side
			                    vp_vec(np.linspace(d, d + 2.*a, N+2, endpoint=True), g, rs + h)# top
			                    ]).ravel()
	ar = np.append(ar, v_domain_min)

	#print(np.minimum( np.max(ar), v_domain_max))

	return np.minimum( np.max(ar), v_domain_max)


def vmin(d, g, rs, h, a):
	return np.maximum(
		   np.min(np.array([vp(d,        g, rs),
		                    vp(d + 2.*a, g, rs),
		                    vp(d,        g, rs + h),
		                    vp(d + 2.*a, g, rs + h),
		                    v_domain_max])),
		   v_domain_min)

def vminN(d, g, rs, h, a, N):

	ar = np.concatenate([vp_vec(np.linspace(d, d + 2.*a, N+2, endpoint=True), g, rs),    # bottom
			                    vp_vec(d, g, np.linspace(rs, rs + h, N+2, endpoint=True)),      # side
			                    vp_vec(np.linspace(d, d + 2.*a, N+2, endpoint=True), g, rs + h)# top
			                    ]).ravel()
	ar = np.append(ar, v_domain_max)

	#print(np.minimum( np.max(ar), v_domain_max))

	return np.maximum( np.min(ar), v_domain_min)


vmaxvec = vectorize(vmax)
vminvec = vectorize(vmin)

vmaxNvec = vectorize(vmaxN)
vminNvec = vectorize(vminN)


def bmax(d, rs):
	b_disc = np.sin(d) * np.sin(d + 2.*rs)

	# factor = 1.
	# if d+rs > np.pi:
	# 	factor = -1.

	if b_disc >= 0.:
		return np.arctan2(np.sin(rs), np.sqrt(b_disc))
	else:
		return np.pi/2.

bmaxvec = vectorize(bmax)

def cosa(d, rs, b):
	sign = 1.
	# if d+rs > np.pi:
	# 	sign = -1.
	return -sign * np.sqrt(1. - ((np.sin(d+rs) * np.sin(b))/(np.sin(rs)))**2)

def ai(d, rs, b):
	return rs*np.abs(cosa(d, rs, b))

def dazm(d, rs, b):
	if d+2.*rs < np.pi:
		return np.arctan2(2.*np.sin(d)*np.sin(d+2.*rs), np.sin(2.*(d+rs))*np.cos(b) - np.sin(2.*rs)*cosa(d, rs, b))
	else:
		return np.pi + np.arctan2(2.*np.sin(d)*np.sin(d+2.*rs), np.sin(2.*(d+rs))*np.cos(b) - np.sin(2.*rs)*cosa(d, rs, b))

# def dazm(d, rs, b):
# 	return np.arctan2(2.*np.sin(d)*np.abs(np.sin(d+2.*rs)), np.abs(np.sin(2.*(d+rs)))*np.cos(b) - np.sin(2.*rs)*cosa(d, rs, b))

dazmvec = vectorize(dazm)


def dM(d, g, rs, h, a):
	return (vminvec(d, g, rs, h, a)**(-3.*mu) - vmaxvec(d, g, rs, h, a)**(-3.*mu)) * np.sin(g)#* np.sin(d) #/(np.pi*(a+h)) #


def int_dM(d, g, rs, h, a):

	int_dM_i = np.zeros(len(d))

	for i in range(len(d)):
		dM_i = dM(d[i], g, rs, h, a) * bmax(d[i], rs) / np.pi * a
		int_dM_i[i] = np.nansum( (g[1:] - g[:-1]) * (dM_i[:-1] + dM_i[1:]) / 2.)

		dM_i = dM(d[i]+np.pi, g, rs, h, a) * bmax(d[i]+np.pi, rs) / np.pi * a
		int_dM_i[i] += np.nansum( (g[1:] - g[:-1]) * (dM_i[:-1] + dM_i[1:]) / 2.)


		# plt.figure()
		# plt.plot(g, dM_i)
		# plt.show()

	return int_dM_i 


def int_dM_azm(d, g, rs, h, a):

	int_dM_i = np.zeros(len(d))

	for i in range(len(d)):
		azm = np.linspace(0., bmax(d[i], rs), Nazm)

		for azmi in azm[:-1]:
			#print(azm/bmax(d[i], rs), ai(d[i], rs, azmi)/rs, dazm(d[i], rs, azmi)/d[i])
			dM_i = dM(dazm(d[i], rs, azmi), g, rs, h, ai(d[i], rs, azmi)) * azm[1] / np.pi * a
			#dM_i = dM(d[i], g, rs, h, a) * bmax(d[i], rs) / np.pi * a
			int_dM_i[i] += np.nansum( (g[1:] - g[:-1]) * (dM_i[:-1] + dM_i[1:]) / 2.)

			dM_i = dM(dazm(d[i]+np.pi, rs, azmi), g, rs, h, ai(d[i]+np.pi, rs, azmi)) * azm[1] / np.pi * a
			#dM_i = dM(d[i]+np.pi, g, rs, h, a) * bmax(d[i]+np.pi, rs) / np.pi * a
			int_dM_i[i] += np.nansum( (g[1:] - g[:-1]) * (dM_i[:-1] + dM_i[1:]) / 2.)


	return int_dM_i / float(Nazm-1)


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
df = float(sys.argv[4])
#gg = float(sys.argv[5])


#vv = vp(df, gg, rr)

# print('v = ', vv)
# print('r = ', r_final(df, vv, 0.51951786))
# print('d = ', dist(vv, gg, rr, -1))


Ng = 2000
Nd = 100
Nh = 20
Na = 10
Nv = 10

Np = 20


gp = np.linspace(EPS, np.pi/2.-EPS, Ng)

#dd = np.logspace(-4., np.log10(np.pi - 10.*af), Nd)
#dd = np.linspace(0., np.pi - af*10, Nd)

#rv = np.linspace(0., 1., 10) + rr
#aa = np.linspace(0., af, Na)

#print(rv-1)
#cdf = np.zeros(Ng+1)
cdf = np.cumsum(vmaxvec(df, gp, rr, hh, af) - vminNvec(df, gp, rr, hh, af, Np))
cdf /= cdf[-1]

#print(cdf)

N = 100000

g_sample = np.zeros(N)
v_sample = np.zeros(N)

x = np.random.uniform(size=N)

for i in range(N):
	idx = np.argmax(cdf >= np.random.uniform(0, 1)) - 1
	g_sample[i] = np.random.uniform(gp[idx], gp[idx+1])
	v_sample[i] = np.random.uniform(vminNvec(df, g_sample[i], rr, hh, af, Np), vmaxvec(df, g_sample[i], rr, hh, af))

# height at the side of the asset
r_si = r_final(df, v_sample, g_sample)

# plt.plot(gp, vmaxNvec(df, gp, rr, hh, af, Np))
# plt.plot(gp, vminNvec(df, gp, rr, hh, af, Np))
# plt.plot(gp, vmaxvec(df, gp, rr, hh, af))
# plt.plot(gp, vminvec(df, gp, rr, hh, af))
# plt.xlim([0., np.pi/2.])
# plt.ylim([0., 1.])



plt.figure()
plt.plot(gp*180./np.pi, vmaxvec(df, gp, rr, hh, af), color='k')
plt.plot(gp*180./np.pi, vminNvec(df, gp, rr, hh, af, Np), color='k')
plt.xlim([0., 90.])
plt.ylim([0., 1.])
plt.xlabel(r'Ejecta Zenith Angle $\gamma_p$ (deg)')
plt.ylabel(r'Ejecta Speed $v_p$ ($v_{esc}$)')
plt.title(r'Asset Size Parameters: a = ' + str(af) + ', h = ' + str(hh) + ', D = ' + str(df))
#plt.plot(gp, cdf)

# side hits
g_side = g_sample[(r_si >= rr) & (r_si <= rr + hh)]
v_side = v_sample[(r_si >= rr) & (r_si <= rr + hh)]

# top hits
g_top = g_sample[r_si > rr + hh]
v_top = v_sample[r_si > rr + hh]

# bottom hits
g_bottom = g_sample[r_si < rr]
v_bottom = v_sample[r_si < rr]


plt.scatter(g_side*180./np.pi, v_side, s=3, color='b', label='Side')
plt.scatter(g_bottom*180./np.pi, v_bottom, s=3, color='r', label='Bottom')
plt.scatter(g_top*180./np.pi, v_top, s=3, color='g', label='Top')
plt.scatter(g_top[dist_new(v_top, g_top, rr, hh, df) < df]*180./np.pi, v_top[dist_new(v_top, g_top, rr, hh, df) < df], s=1, color='r')
plt.scatter(g_top[dist_new(v_top, g_top, rr, hh, df) > df + 2.*af]*180./np.pi, v_top[dist_new(v_top, g_top, rr, hh, df) > df + 2.*af], s=1, color='y')
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.legend()

# https://stackoverflow.com/questions/339007/how-to-pad-zeroes-to-a-string
plt.savefig(f'asset_speed_zenith_plot_{rr:.3e}_{af:.3e}_{hh:.3e}_{df:.3e}.png', bbox_inches='tight', dpi=600)


plt.figure()
for di in np.linspace(df, df + 2.*af, 100):
	plt.plot(gp, vp_vec(di, gp, rr+hh), color='g')
	#print(di, dist_new(vp_vec(di, gp, rr+hh), gp, rr, hh, df))
for di in np.linspace(df, df + 2.*af, 100):
	plt.plot(gp, vp_vec(di, gp, rr), color='r')
for ri in np.linspace(rr, rr+hh, 100):
	plt.plot(gp, vp_vec(df, gp, ri), color='b')
plt.xlim([0., np.pi/2.])
plt.ylim([0., 1.])
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python


#print(g_bottom[0], v_bottom[0], dist(g_bottom[0], v_bottom[0], rr), F(v_bottom[0], g_bottom[0], rr))


############
# plt.figure()
# plt.scatter((np.ones(N) * df)[(r_si >= rr) & (r_si <= rr + hh)], r_si[(r_si >= rr) & (r_si <= rr + hh)], s=1.5, color='b')
# if rr > 1.:
# 	plt.scatter(dist_new(v_bottom, g_bottom, rr, hh, df), (np.ones(np.size(g_bottom)) * rr), s=1.5, color='r')
# plt.scatter(dist_new(v_top, g_top, rr, hh, df), (np.ones(np.size(g_top)) * (rr + hh)), s=1.5, color='g')
###########


# plt.axhline(y = rr, color='r', linestyle='-') # https://stackoverflow.com/questions/33382619/plot-a-horizontal-line-using-matplotlib
# plt.axhline(y = rr + hh, color='g', linestyle='-') 

##########################################################3
# F0 = int_dM(dd, gp, rr, hh, af)

# plt.loglog(dd, F0, label=r'$v\in $' + f'({v_domain_min:.2f}, {v_domain_max:.2f})')

# v = np.linspace(v_domain_min, v_domain_max, Nv)
# for i in range(Nv-1):
# 	v_domain_max = v[i+1]
# 	v_domain_min = v[i]
# 	F0 = int_dM(dd, gp, rr, hh, af)
# 	plt.loglog(dd, F0, label=r'$v\in $' + f'({v_domain_min:.2f}, {v_domain_max:.2f})')
##########################################################3

# plt.loglog(dd, int_dM(dd, gp, rr, hh, af))
# plt.loglog(dd, int_dM_azm(dd, gp, rr, hh, af))

#print(np.nansum((dd[1:]-dd[:-1])*(F0[:-1]+F0[1:])/2.) / np.nansum((dd[1:]-dd[:-1])*(F1[:-1]+F1[1:])/2.))


############
# for ri in rv:
# 	plt.loglog(dd/np.pi, int_dM(dd, gp, ri, hh, af))
plt.legend()
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python

# plt.figure()
# plt.plot(rv-1, iint_dM_r(dd, gp, rv, hh, af))


# plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
# plt.figure()


# plt.plot(dd, dazmvec(dd, rr, bmaxvec(dd, rr)))

plt.show()


