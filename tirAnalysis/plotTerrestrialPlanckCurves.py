from tirAnalysis.EMRadiation import planckFunction, weinDisplacementLaw

import numpy as np
import matplotlib.pyplot as plt

plt.close()
mu = np.multiply(np.linspace(2, 30, 2801) + 1e-3, 1e-6)
T = np.array([265, 273.15, 280, 290, 300])

fig, ax = plt.subplots()

for T_ in T:
     ax.plot(mu, planckFunction(mu, T_)*1e-6, '-', label=r'T = %6.2f K' % T_)


def xtickFn(x):
     x *= 1e6                                                          
     return ["%d" % np.round(x_) for x_ in x]


ax.set_xticklabels(xtickFn(ax.get_xticks())) 
ax.set_xlabel(r'Wavelength/[$\mu$m]')
ax.set_ylabel(r'Radiance/[W/(m$^2$ sr $\mu$m)]')
ax.grid(True)

_T = np.linspace(265, 300, 101)
planckPeakLocs = [planckFunction(mu, T_).max()*1e-6 for T_ in _T]
ax.plot(weinDisplacementLaw(_T), planckPeakLocs, 'k--', 
        label='Wein Displacement Law')

ax.legend(loc=0)

# Save figure
fig.tight_layout()
fig.savefig('/home/david/Documents/programming/python/EMRadiation_planckCurvesTerrestrialTemps.pdf', bbox_inches='tight')
