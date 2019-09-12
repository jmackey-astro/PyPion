# Author: Sam Green and Maggie Goulden
# Created: 15-02-2019

# Script that plots 2D data from the Nested_Grid_Pion.

# Edits:
# - 15-02-2019 MM: Fixed axis bug.
# - 18-02-2019 SG: Added in new variables, var1 and var2, to be transferred to the Plotting_Classes.py script to give
# information on what data to be plotted. This should allow the user to pick what data needs to be plotted from this
# script alone.
# - 18-02-2019 SG: Added comments to try and explain things.
# - 12-08-2019 SG: Code now using multiple cores.

from Plotting_Classes import Plotting2d
import matplotlib.pyplot as plt
import astropy.units as u
import time

import sys
import warnings
import matplotlib.cbook

sys.path.append('.')
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

imag = sys.argv[1]

arr = []
for i in range(2, len(sys.argv)):
    arr.append(sys.argv[i])

# ################# This is where you change all the parameters. #######################
# Expected layout: ["variable", max colour value, min colour value, colormap]
var1 = ["Density", -21, -27, "viridis", 'log']
var2 = ["Temperature", 8.2, 4, "gist_heat", 'log']
# Variable choices: "Density", "DivB", "MagneticFieldX", "MagneticFieldY", "Pressure", "Temperature", "VelocityX",
# "VelocityY", "VelocityZ", "glmPSI".
# #######################################################################################
t0 = time.time()

# First loop goes through each set of files for each timestep.
# Function to take all the layer files and combine into one figure
fig = plt.figure()
imagefile = "%szeta_%s.png" % (imag, arr[0][len(arr[0]) - 22:len(arr[0]) - 6])

# Define the max extents of the figure grid.
lim = Plotting2d(arr[0])
lim_min = lim.xmin().to(u.pc)
lim_max = lim.xmax().to(u.pc)
lim.close()
del lim

# Second loop goes through each file in that set and plots it to a figure.
for f in arr:
    level = Plotting2d(f)
    level.plotsilo_2d(lim_min, lim_max, fig, var1, var2)

    level.close()
    del level


# This now saves the figure that has had multiple levels plotted onto it.
plt.savefig(imagefile, bbox_inches='tight', dpi=300)
plt.close()

del fig, imagefile

t1 = time.time()
print("Running Time: {0}".format(t1 - t0))
