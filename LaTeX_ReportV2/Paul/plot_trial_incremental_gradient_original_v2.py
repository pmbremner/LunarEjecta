'''
Test methods to stich together scalebars with
a gradient within each discrete color increment.
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from matplotlib import gridspec
import colorsys



#=========================================================================
#
#       CREATE TEST DATA (td)
#
#=========================================================================

numvals = 10

td_x1 = np.linspace(0.0,100.0,num=numvals, endpoint=True, dtype=np.float)
td_x2 = np.linspace(101.0,200.0,num=numvals, endpoint=True, dtype=np.float)
td_x3 = np.linspace(201.0,300.0,num=numvals, endpoint=True, dtype=np.float)

td_y1 = 1.0 * np.ones(numvals, dtype=np.float)
td_y2 = 10.0 * np.ones(numvals, dtype=np.float)
td_y3 = 20.0 * np.ones(numvals, dtype=np.float)

#====================================================================
#
#       ATTEMPT TO MAKE CUSTOM COLOR SCALE
#
#====================================================================

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(colors.to_rgb(c1))
    c2=np.array(colors.to_rgb(c2))
    return colors.to_hex((1-mix)*c1 + mix*c2)

c1='#1f77b4' #blue
c2='green' #green #hsv=(120,1.0,0.5)
n=300

for x in range(n+1):
    clrs_1=colorFader(c1,c2,x/n)


def colorFader2(h,v,ninc):
    cgrad = [colorsys.hsv_to_rgb(h,s,v) for s in np.linspace(0.25,1.0,num=ninc, endpoint=True, dtype=np.float)]
    return cgrad

greenfade = colors.LinearSegmentedColormap.from_list('greenfade', colorFader2(2./6,0.5,30), N=10)
bluefade = colors.LinearSegmentedColormap.from_list('bluefade', colorFader2(4./6,0.5,30), N=10)
purplefade = colors.LinearSegmentedColormap.from_list('purplefade', colorFader2(5./6,0.5,30), N=10)

#====================================================================
#
#       SETUP THE FIGURE LAYOUT
#
#====================================================================

# Set the figure fonts
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica-Bold']})
rc('text', usetex=True)

# Create a list of color maps
#cmap_names = ['winter', 'cool','Wistia']
cmap_names = [greenfade, bluefade, purplefade]


totalfig = plt.figure(num=1,figsize=(10,10))

gs_plots = gridspec.GridSpec(1,1)
gs_plots.update(left=0.10, right=0.85, top=0.95, bottom=0.12, hspace=0.04)

gs_bar = gridspec.GridSpec(3,1, hspace=0.0, wspace=0.0)
#gs_bar = gridspec.GridSpec(3,2, hspace=0.0, wspace=0.0)
gs_bar.update(left=0.86, right=0.96,top=0.95,bottom=0.12)

# Plot Panels
mm_ax = totalfig.add_subplot(gs_plots[0,0])
#clrbar_label = totalfig.add_subplot(gs_bar[:3,0])
clrbar_ax1 = totalfig.add_subplot(gs_bar[2,0])
clrbar_ax2 = totalfig.add_subplot(gs_bar[1,0])
clrbar_ax3 = totalfig.add_subplot(gs_bar[0,0])


#====================================================================
#
#       PROCESS THE DATA
#
#====================================================================

PLOTNAME = f'incremental_gradient_test'

# Plot current dataset
plt1 = mm_ax.scatter(td_x1, td_y1, marker='o', c=td_x1, cmap=cmap_names[0], alpha=1.0, linewidth=3.0, s=100)
plt2 = mm_ax.scatter(td_x2, td_y2, marker='o', c=td_x2, cmap=cmap_names[1], alpha=1.0, linewidth=3.0, s=100)
plt3 = mm_ax.scatter(td_x3, td_y3, marker='o', c=td_x3, cmap=cmap_names[2], alpha=1.0, linewidth=3.0, s=100)


mm_ax.set_ylabel('y-values', fontsize=24, fontweight='bold')
mm_ax.set_xlabel('x-values', fontsize=24, fontweight='bold')
mm_ax_yticks = mm_ax.get_yticks() #(minor=False)
mm_ax.set_yticklabels(mm_ax_yticks, fontsize=18, fontweight='regular')
mm_ax_xticks = mm_ax.get_xticks() #(minor=False)
mm_ax.set_xticklabels(mm_ax_xticks, fontsize=18, fontweight='regular',rotation=90)


#=========================================================================
#
#       MAKE THE COLOR BAR
#
#=========================================================================

clrbar_ax2.set_ylim([0,10])
clrbar_ax2.set_xlim([0,10])

# Create and label the colorbar
totalfig.colorbar(plt1,ax=clrbar_ax1)
totalfig.colorbar(plt2,ax=clrbar_ax2)
totalfig.colorbar(plt3,ax=clrbar_ax3)

clrbar_ax2.text(10.0,5.0,'The Label',fontsize=30,fontweight='heavy',color='k',rotation=90,verticalalignment='center',horizontalalignment='right')

# Turn off the axes
clrbar_ax1.axis('off')
clrbar_ax2.axis('off')
clrbar_ax3.axis('off')
#______________________________________________________________

plt.show()
# #for fmt in ['png','eps','jpg','svg']:
# #for fmt in ['png','eps','svg']:
# for fmt in ['png']: # temporarily removed jpg from list due to matplotlib bug
#     totalfig.savefig(str(cwd.joinpath('plots',f'{PLOTNAME}.{fmt}')),
#                 facecolor='w', edgecolor='w',
#                 orientation='portrait', papertype=None,
#                 format=fmt, transparent=False,
#                 bbox_inches='tight', pad_inches=0.1
#                 )
# plt.close(totalfig)

print (f'\nDone!\n')