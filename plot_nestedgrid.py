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
from argparse_command import InputValues
from colormaps import cmaps
import matplotlib.pyplot as plt
import astropy.units as u
import multiprocessing as mp
from joblib import Parallel, delayed
import dill
import time

# Saves the relevant information entered in the command.
line = InputValues()
time_dicts = line.time_dicts


# ################# This is where you change all the parameters. #######################
# Expected layout: ["variable", max colour value, min colour value, colormap]
var1 = ["Density", -21, -27, "viridis", 'log']
var2 = ["Temperature", 8.2, 4, "gist_heat", 'log']
# Variable choices: "Density", "DivB", "MagneticFieldX", "MagneticFieldY", "Pressure", "Temperature", "VelocityX",
# "VelocityY", "VelocityZ", "glmPSI".
# #######################################################################################


def plotter(files):  # for files in time_dicts:
    # First loop goes through each set of files for each timestep.
    # Function to take all the layer files and combine into one figure
    fig = plt.figure()
    imagefile = "%s%s_%s.png" % (line.img_path, line.img_file, time_dicts[files][0][len(time_dicts[files][0]) - 22:
                                                                                    len(time_dicts[files][0]) - 6])
    # Define the max extents of the figure grid.
    lim_min = Plotting2d(time_dicts[files][0]).xmin().to(u.pc)
    lim_max = Plotting2d(time_dicts[files][0]).xmax().to(u.pc)

    # Second loop goes through each file in that set and plots it to a figure.
    for f in time_dicts[files]:
        Plotting2d(f).plotsilo_2d(lim_min, lim_max, fig, var1, var2)

    # This now saves the figure that has had multiple levels plotted onto it.
    plt.savefig(imagefile, bbox_inches='tight', dpi=300)
    plt.close()


# To run the code with multiple processors:
num_cores = mp.cpu_count()
Parallel(n_jobs=num_cores)(delayed(plotter)(files) for files in time_dicts)
