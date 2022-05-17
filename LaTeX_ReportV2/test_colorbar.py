import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import colorsys
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement


rgb_to_hsv_vec = vectorize(colorsys.rgb_to_hsv)



def colorFader2(h, v, Ns, s0, s1):
    cgrad = [colorsys.hsv_to_rgb(h, s, v) for s in np.linspace(s0, s1, num=Ns, endpoint=True, dtype=np.float)]

    # https://numpy.org/doc/stable/reference/generated/numpy.column_stack.html
    cgrad = np.column_stack((cgrad, np.ones(Ns))) # for alpha column
    return cgrad


def colorFader3(h, s, Nv, v0, v1):
	cgrad = [colorsys.hsv_to_rgb(h, s, v) for v in np.linspace(v0, v1, num=Nv, endpoint=True, dtype=np.float)]

    # https://numpy.org/doc/stable/reference/generated/numpy.column_stack.html
	cgrad = np.column_stack((cgrad, np.ones(Nv))) # for alpha column
	return cgrad


def cm_varsat(cm_name, Nc, Ns, s0=0.5, s1=1.):
	# get an array of Nc number of rbga colors from the cm_name built-in colormap
	cur_cmap = cm.get_cmap(cm_name, Nc) # rgba, effectively down-samples to Nc from cmap.N

	cur_cm = cur_cmap(np.arange(0, cur_cmap.N))

	# convert the rgba space to hsva space
	# https://docs.python.org/3/library/colorsys.html#module-colorsys
	cur_cm[:,0:3] = np.transpose(rgb_to_hsv_vec(cur_cm[:,0], cur_cm[:,1], cur_cm[:,2]))

	newcolormap = []
	for i, cur_cmi in enumerate(cur_cm):
		if i == 0:
			newcolormap = colorFader2(cur_cmi[0], cur_cmi[2], Ns, s0*cur_cmi[1], s1*cur_cmi[1])
		else:
			# https://towardsdatascience.com/creating-colormaps-in-matplotlib-4d4de78a04b8
			newcolormap = np.vstack((newcolormap, colorFader2(cur_cmi[0], cur_cmi[2], Ns, s0*cur_cmi[1], s1*cur_cmi[1]) ))

	return ListedColormap(newcolormap, name = cm_name + '_varsat')


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


new_cm = cm_varval('tab20', 20, 50) # 'viridis' plasma


X, Y = np.meshgrid(np.linspace(0,1000.,1000), np.linspace(0,1000.,1000))

plt.pcolormesh(X*Y, cmap = new_cm)
plt.colorbar()
plt.show()