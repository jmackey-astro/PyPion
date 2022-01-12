import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator

import numpy as np
from numpy import *
import math

E = np.arange(0.2, 2.1, 0.1)
N_H = np.array([2.7e21, 3.7e21, 4.7e21, 5.7e21, 6.7e21, 7.7e21])


sigma = 2.27e-22 * E**(-2.485)

tau_2p7 = N_H[0] * sigma
tau_3p7 = N_H[1] * sigma
tau_4p7 = N_H[2] * sigma
tau_5p7 = N_H[3] * sigma
tau_6p7 = N_H[4] * sigma
tau_7p7 = N_H[5] * sigma

x_2p7 = np.zeros(len(tau_2p7))
x_3p7 = np.zeros(len(tau_3p7))
x_4p7 = np.zeros(len(tau_4p7))
x_5p7 = np.zeros(len(tau_5p7))
x_6p7 = np.zeros(len(tau_6p7))
x_7p7 = np.zeros(len(tau_7p7))

for i in range(len(tau_2p7)):
    x_2p7[i] = math.exp(tau_2p7[i])
    x_3p7[i] = math.exp(tau_3p7[i])
    x_4p7[i] = math.exp(tau_4p7[i])
    x_5p7[i] = math.exp(tau_5p7[i])
    x_6p7[i] = math.exp(tau_6p7[i])
    x_7p7[i] = math.exp(tau_7p7[i])


fig, ax1 = plt.subplots()

ax1.plot(E, x_2p7, color='blue', linestyle='solid', label='N_H = 2.7e21')
ax1.plot(E, x_3p7, color='yellow', linestyle='solid', label='N_H = 3.7e21')
ax1.plot(E, x_4p7, color='red', linestyle='solid', label='N_H = 4.7e21')
ax1.plot(E, x_5p7, color='purple', linestyle='solid', label='N_H = 5.7e21')
ax1.plot(E, x_6p7, color='cyan', linestyle='solid', label='N_H = 6.7e21')
ax1.plot(E, x_7p7, color='black', linestyle='solid', label='N_H = 7.7e21')

ax1.set_ylabel(r'$e^\tau$')
ax1.set_xlabel('E (kev)')

# ax1.set_yscale('log')
ax1.set_ylim((1, 5))
ax1.legend(loc='upper right', ncol=1, framealpha=1)
ax1.set_xlim((0.2, 2))
ax1.grid(True)

# Save the figures as .png files:
plt.savefig("/home/green/bubble_sims/N_H.png",
            bbox_inches='tight', dpi=300)

#plt.show()