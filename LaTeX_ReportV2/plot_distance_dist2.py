#!/usr/bin/python

import sys
import numpy as np
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement
import matplotlib.pyplot as plt
from matplotlib import cm

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
import colorsys

rgb_to_hsv_vec = vectorize(colorsys.rgb_to_hsv)

###################################################################

def colorFader3(h, s, Nv, v0, v1):
	cgrad = [colorsys.hsv_to_rgb(h, s, v) for v in np.linspace(v0, v1, num=Nv, endpoint=True, dtype=np.float)]

    # https://numpy.org/doc/stable/reference/generated/numpy.column_stack.html
	cgrad = np.column_stack((cgrad, np.ones(Nv))) # for alpha column
	return cgrad

def cm_varval(cm_name, Nc, Nv, v0=0.35, v1=1.):
	# get an array of Nc number of rbga colors from the cm_name built-in colormap
	cur_cmap = cm.get_cmap(cm_name, Nc) # rgba, effectively down-samples to Nc from cmap.N

	cur_cm = cur_cmap(np.arange(0, cur_cmap.N))

	# convert the rgba space to hsva space
	# https://docs.python.org/3/library/colorsys.html#module-colorsys
	cur_cm[:,0:3] = np.transpose(rgb_to_hsv_vec(cur_cm[:,0], cur_cm[:,1], cur_cm[:,2]))

	newcolormap = []
	for i, cur_cmi in enumerate(cur_cm):
		if i == 0:
			newcolormap = colorFader3(cur_cmi[0], cur_cmi[1], Nv, v0*cur_cmi[2], v1*cur_cmi[2])
		else:
			# https://towardsdatascience.com/creating-colormaps-in-matplotlib-4d4de78a04b8
			newcolormap = np.vstack((newcolormap, colorFader3(cur_cmi[0], cur_cmi[1], Nv, v0*cur_cmi[2], v1*cur_cmi[2]) ))

	return ListedColormap(newcolormap, name = cm_name + '_varval')

## Example ##
#new_cm = cm_varval('tab20', 20, 50) # 'viridis' plasma

###################################################################



def v_speed(d, g, rs):
	return (((1./rs - np.cos(d))/(1. - np.cos(d)))*(1. - np.cos(2.*g)) + np.sin(2.*g)/np.tan(d/2.))**-(1/2)
	#return ((((1. - rs)/rs + 2.*np.sin(d/2.))/(np.sin(d/2.)))*(np.sin(g)) + np.sin(2.*g)/np.tan(d/2.))**-(1/2)

def gamma_pp(d, v, rs):
	return np.arctan(1./(v**2/np.tan(d/2) + np.sqrt(v**4 / np.tan(d/2)**2 + 2.*v**2 *((1./rs - np.cos(d)) / (1. - np.cos(d))) - 1.) ))

def gamma_pm(d, v, rs):
	return np.arctan(1./(v**2/np.tan(d/2) - np.sqrt(v**4 / np.tan(d/2)**2 + 2.*v**2 *((1./rs - np.cos(d)) / (1. - np.cos(d))) - 1.) ))

def F(v, g, rs):
	return -np.sqrt(1. + (rs - 1.)/(np.cos(g)**2) * (rs + 1. - rs/v**2))

def g_ap(d, rs):
	return np.arctan2(rs*np.tan(d/2.), rs - 1.)

def d_dist(v, g, rs, d0):
	Fi = F(v, g, rs)
	if g > g_ap(d0, rs) and d0 < np.pi:
		Fi *= -1
		if v**2 > 1./(1. + 1./rs):
			return 2.*np.mod(np.pi + np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) , 2.*np.pi)
		else:
			return 2.*np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) 
	return 2.*np.arctan2(2.*v**2 * np.sin(g)*np.cos(g) + ((rs-1.)/(1.-Fi))*(2.*v**2-1.)*np.tan(g), ((rs-Fi)/(1.-Fi))-2.*(v*np.sin(g))**2 ) 


# def d_dist(v, g, rs, d0):
# 	Fi = F(v, g, rs)
# 	if g > g_ap(d0, rs):
# 		Fi *= -1
# 	# if d0 > np.pi:
# 	# 	return 2.*np.mod(2.*np.pi + np.arctan2(v**2 * np.sin(g)*np.cos(g) * (1. - Fi/rs), 1. - (v*np.sin(g))**2*(1.+1./rs)), 2.*np.pi)
# 	# else:
# 	return np.mod(2.*np.pi + 2.*np.arctan2(v**2 * np.sin(g)*np.cos(g) * (1. - Fi/rs), 1. - (v*np.sin(g))**2*(1.+1./rs)), 2.*np.pi)



d_distvec = vectorize(d_dist)

def v_min(g, rs):
	return np.sqrt((2.*rs*(rs-1.)) / (np.cos(2.*g)-1. + 2.*rs**2))

def Gg(rs, d):
	return 1. - 1/rs + 2.*np.sin(d/2) * (np.sqrt(1. - (1.-1./rs)*(1.-(1.-1./rs)/(4.*np.sin(d/2)**2))) - 1.)


def g_opt(d, rs):
	return np.arctan2(1., np.tan(np.pi/4. - d/4.) + Gg(rs, d)/np.sin(d))


def g_final(v, g, rs, d0):
	Fi = F(v, g, rs)
	if g > g_ap(d0, rs) and d0 < np.pi:
		Fi *= -1
	return np.arctan2(np.tan(g), -Fi)

g_finalvec = vectorize(g_final)


def v_final(v, rs):
	return np.sqrt(1./rs + v**2 - 1.)


rs  = float(sys.argv[1])
dd  = float(sys.argv[2])
ii  = int(sys.argv[3])


print(d_dist(0.8, 80./180.*np.pi, rs, dd))

N = 700

g = np.linspace(0.0001, np.pi/2*0.9999, N)
v = np.linspace(0.0001, 1., N)
G, V = np.meshgrid(g, v)


D = d_distvec(V, G, rs, dd)

#D0 = d_dist(V, G, 1.)

fig = plt.figure(figsize=(12, 8))
axs = plt.axes()


cm_val1 = cm_varval('tab20', 20, 128)

plt.pcolormesh(G/np.pi*180., V, D/np.pi, cmap=cm_val1)#,shading='gouraud') #gist_ncar, nipy_spectral, gist_rainbow, plasma, prism


vi = v_speed(dd, g, rs)

#print(vi)

gfinal_i = g_finalvec(vi[~np.isnan(vi)], g[~np.isnan(vi)], rs, dd)

vfinal_i = v_final(vi[~np.isnan(vi)], rs)

plt.plot(g/np.pi*180., vi, color='black', linewidth=10, label=f'Impact at height {rs-1.:.3e}' + r' $r_m$' + f' | distance {dd/np.pi:.3f}' + r' $\pi r_m$', zorder=4)


# https://matplotlib.org/3.5.0/gallery/lines_bars_and_markers/multicolored_line.html
# https://stackoverflow.com/questions/19390895/matplotlib-plot-with-variable-line-width
#########################################################

x = g[~np.isnan(vi)]/np.pi*180.
y = vi[~np.isnan(vi)]

yc = gfinal_i/np.pi*180.

lw = vfinal_i / y * 10
#yc = np.cos(0.5 * (x[:-1] + x[1:]))  # first derivative

# Create a set of line segments so that we can color them individually
# This creates the points as a N x 1 x 2 array so that we can stack points
# together easily to get the segments. The segments array for line collection
# needs to be (numlines) x (points per line) x 2 (for x and y)
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

#fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)

# Create a continuous norm to map from data points to colors
cm_val2 = cm_varval('Set1', 9, 4)

norm = plt.Normalize(0., 180.)
lc = LineCollection(segments, linewidths=lw, cmap=cm_val2, norm=norm, zorder=5) #seismic, bwr, brg
# Set the values used for colormapping
lc.set_array(yc)
#lc.set_linewidth(3)
line = axs.add_collection(lc)
fig.colorbar(line, ax=axs, orientation="vertical").set_label(label=r'Final Incoming Zenith Angle $\gamma_s$ (deg)' + '\n' + r'Line Thickness = Final/Initial Speeds ($v_s/v_p$)', size=15)



############################################################


#plt.plot(g, v_speed(dd+1.E-5, g, rs), color='black', linewidth=1)
#plt.plot(g, v_speed(dd, g, rs+0.0001), color='black', linewidth=1)

di = np.linspace(0.000001, np.pi*0.99, 1000)

g_opti = g_opt(di, rs)

plt.plot(g_opti/np.pi*180., v_speed(di, g_opti, rs), color='black', linewidth=4, linestyle=":", label='Optimal Zenith Angle $\gamma_{p,opt}$')
plt.plot(g/np.pi*180., v_min(g, rs), color='black', linewidth=4, linestyle='--', label=r'Minimum Speed $v_{min}$')
plt.ylim(0, 1)
plt.xlim(0, 90.) #np.pi/2.
#plt.title()
plt.xlabel(r'Initial Zenith Angle $\gamma_p$ (deg)', size=15)
plt.ylabel(r'Initial Speed $v_p$ ($V_{esc}$)', size=15)
plt.grid(b=True, which='both') # https://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
plt.legend(handlelength=4, loc='upper left', bbox_to_anchor=(0., 1.13), borderaxespad=0.0) # https://stackoverflow.com/questions/39824599/python-matplotlib-legend-linestyle
plt.colorbar(orientation="vertical").set_label(label=r'Distance $D$ ($\pi r_m$)', size=15)
plt.clim(0,2) #https://stackoverflow.com/questions/3373256/set-colorbar-range-in-matplotlib


# https://stackoverflow.com/questions/339007/how-to-pad-zeroes-to-a-string
plt.savefig(f'dist_speed_zenith_plot_{ii:03}_{rs-1.:.3e}_{dd:.3f}.png', bbox_inches='tight', dpi=600)

# plt.figure()
# plt.plot(g[~np.isnan(vi)]/np.pi*180., gfinal_i/np.pi*180.)


plt.show()