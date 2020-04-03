import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def api(ai):
	return 0.0003*ai**3 - 0.036*ai**2 + 1.5206*ai + 20

def a0(ai):
	return -0.00042*ai**3 + 0.0236*ai**2 + 0.129*ai + 20

def amax(b, ai):
	return api(ai)*np.sin(b/2)**2 + a0(ai)*np.cos(b/2)**2

N = 200

b = np.linspace(0, 180, N)

ai = np.linspace(0, 90, N)

azm = np.zeros(N)

for i in range(0, N):
	if amax(0, ai[i]) * amax(np.pi, ai[i]) < 0:
		azm[i] = 180./np.pi*optimize.bisect(amax, 0, np.pi, args=(ai[i]))
	else:
		azm[i] = 0

print(azm)

fig, ax = plt.subplots()
ax.plot(ai, azm)
ax.set_xlabel('Impact Zenith Angle [Degrees]')
ax.set_ylabel('Exlusion Azimuth Range at Peak Ejecta Angle [Degrees]')

ax.set_ylim(0,90)
ax.set_xlim(0,90)
ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(which='major', linewidth='1')
ax.grid(which='minor', linewidth='0.5', linestyle=':')

plt.savefig('ExlusionRange_vs_ImpactAngle.png', dpi=400)
plt.show()

# plt.plot(b, amax(b/180.0*np.pi, 70))
# plt.show()