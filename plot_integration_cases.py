from matplotlib import pyplot as plt
import numpy as np

fig, axs = plt.subplots(3, 4)

x = np.array((1,1, 0.5,1, 0.5,1, 0.25,0.75, 0,1, 0,1, 0,0.5, 0,1, 0,1, 0,0.5, 0,0))

y = np.array((1,1, 1,0.5, 1,0, 1,0, 1,0.5, 1,0, 1,0, 0.75,0.25, 0.5,0, 0.5,0, 0,0))



for i in range(0,12):
	if(i<11):
		axs[i//4, i%4].plot((0,1),(1,1),'b')
		axs[i//4, i%4].plot((0,1),(0,0),'b')
		axs[i//4, i%4].plot((0,0),(0,1),'b')
		axs[i//4, i%4].plot((1,1),(0,1),'b')

		axs[i//4, i%4].plot(x[2*i:2*i+2], y[2*i:2*(i+1)], '-ro', markersize=8)

		axs[i//4, i%4].set_ylim(-0.25, 1.25)
		axs[i//4, i%4].set_xlim(-0.25, 1.25)
		axs[i//4, i%4].set_title('Case ' + str(i))
	axs[i//4, i%4].tick_params(
	    axis='both',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom=False,      # ticks along the bottom edge are off
	    top=False,         # ticks along the top edge are off
	    left=False,
	    right=False,
	    labelleft=False,
	    labelbottom=False) # labels along the bottom edge are off

#plt.plot(1,1, 'r.', markersize=12)



plt.show()