# Author: Sam Green
# Created: 20-08-2019

# Script that plots 3D data from Nested_Grid_Pion.

from Plotting_Classes import Plotting3d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
import time

import warnings
import matplotlib.cbook

import sys
sys.path.append('.')

warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

imag = sys.argv[1]

arr = []
for i in range(2, len(sys.argv)):
    arr.append(sys.argv[i])


var1 = ["Density", -21, -27, "viridis", 'log', 'y', 127]

# for files in InputValues().time_dicts:
t0 = time.time()

fig = plt.figure()
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.05], width_ratios=[1, 1])
gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)
imagefile = "%sfig%s.png" % (imag, arr[0][len(arr[0]) - 13:len(arr[0]) - 6])

lim = Plotting3d(arr[0])
lim_min = lim.xmin().to(u.pc)
lim_max = lim.xmax().to(u.pc)
lim.close()
del lim

for f in arr:
    level = Plotting3d(f)
    level.plotsilo_3dslice(lim_min, lim_max, fig, gs, var1)

    level.close()
    del level

plt.savefig(imagefile, bbox_inches='tight', dpi=300)
plt.close()

del fig, gs, imagefile

t1 = time.time()
print("Running Time: {0}".format(t1 - t0))
