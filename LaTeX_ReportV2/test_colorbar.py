import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import colorsys
from numpy import vectorize # https://stackoverflow.com/questions/8036878/function-of-numpy-array-with-if-statement


rgb_to_hsv_vec = vectorize(colorsys.rgb_to_hsv)



def colorFader2(h, v, Ns, s0, s1):
    cgrad = [colorsys.hsv_to_rgb(h,s,v) for s in np.linspace(s0, s1, num=Ns, endpoint=True, dtype=np.float)]

    # https://numpy.org/doc/stable/reference/generated/numpy.column_stack.html
    cgrad = np.column_stack((cgrad, np.ones(Ns))) # for alpha column
    return cgrad


def cm_varsat(cm_name, Nc, Ns, s0=0.5, s1=1.):
	# get an array of Nc number of rbga colors from the cm_name built-in colormap
	cur_cmap = cm.get_cmap(cm_name, Nc) # rgba, effectively down-samples to Nc from cmap.N

	cur_cm = cur_cmap(np.arange(0, cur_cmap.N))
	print(np.shape(rgb_to_hsv_vec(cur_cm[:,0], cur_cm[:,1], cur_cm[:,2])))
	# convert the rgba space to hsva space
	# https://docs.python.org/3/library/colorsys.html#module-colorsys
	cur_cm[:,0:3] = np.transpose(rgb_to_hsv_vec(cur_cm[:,0], cur_cm[:,1], cur_cm[:,2]))
	print(cur_cm)
	# 
	#
	newcolormap = []
	for i, cur_cmi in enumerate(cur_cm):
		if i == 0:
			newcolormap = colorFader2(cur_cmi[0], cur_cmi[2], Ns, s0, s1)
		else:
			# https://towardsdatascience.com/creating-colormaps-in-matplotlib-4d4de78a04b8
			newcolormap = np.vstack((newcolormap, colorFader2(cur_cmi[0], cur_cmi[2], Ns, s0, s1)))


	#print(np.shape(newcolormap))

	return ListedColormap(newcolormap, name = cm_name + '_varsat')


new_cm = cm_varsat('viridis', 18, 4)


X, Y = np.meshgrid(np.linspace(0,100.,100), np.linspace(0,100.,100))

plt.pcolormesh(X*Y, cmap = new_cm)
plt.colorbar()
plt.show()