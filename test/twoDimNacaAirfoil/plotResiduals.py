################################################################################
# plotResiduals.py
#   Creator: Ho Mun Onn, Joel
#
# Plots residual output from BES
# Written for python3 with qt module installed
################################################################################

import matplotlib
matplotlib.use('qt5agg') # Comment this out if qt is not available and tkinter is installed

import matplotlib.pyplot as plt
import numpy as np

residuals = np.genfromtxt('output/residual.csv', delimiter=',', skip_footer=1)

m, n = residuals.shape
if n < 9:
  residualHeaders = ['rho', 'rhoU', 'rhoV', 'rhoE']
else:
  residualHeaders = ['rho', 'rhoU', 'rhoV', 'rhoW','rhoE']
  
f = plt.figure()
ax = f.add_subplot(111)
for i, header in enumerate(residualHeaders):
  ax.plot(residuals[:, 0], residuals[:, i+1], label=header)

ax.set_title('Residuals')
ax.set_xlabel('Iteration no.')
ax.set_ylabel('Absolute residual')
ax.set_yscale('log')
ax.grid(True)
ax.legend()
plt.show()
