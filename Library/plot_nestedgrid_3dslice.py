# Author: Sam Green
# Created: 12-08-2019

# Script that plots 3D data from the Nested_Grid_Pion.

from Plotting_Classes import Plotting3d
from argparse_command import InputValues
from colormaps import cmaps
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
import time
import multiprocessing as mp
from joblib import Parallel, delayed
import dill

from pympler.tracker import SummaryTracker
import objgraph

import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

# tracker = SummaryTracker()

# Saves the relevant information entered in the command.
# line = InputValues()
# time_dicts = line.time_dicts

# New method to choose what data gets plotted, everything can be edited in this script, there should no longer be a
# need to edit the Plotting_Classes.py script.

# ################# This is where you change all the parameters. #######################
# Expected layout: ["variable", max colour value, min colour value, colormap]
var1 = ["Density", -21, -27, "viridis", 'log', 'y', 127]
# Variable choices: "Density", "DivB", "MagneticFieldX", "MagneticFieldY", "Pressure", "Temperature", "VelocityX",
# "VelocityY", "VelocityZ", "glmPSI".
# #######################################################################################
# objgraph.show_growth(limit=10)


# def plotter(files):
for files in InputValues().time_dicts:
    t0 = time.time()

    # First loop goes through each set of files for each timestep.
    # Function to take all the layer files and combine into one figure
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.05], width_ratios=[1, 1])
    gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)
    imagefile = "%s%s_%s.png" % (InputValues().img_path, InputValues().img_file,
                                 InputValues().time_dicts[files][0][len(InputValues().time_dicts[files][0]) -
                                                                    22:len(InputValues().time_dicts[files][0]) - 6])
    # Define the max extents of the figure grid.
    lim = Plotting3d(InputValues().time_dicts[files][0])
    lim_min = lim.xmin().to(u.pc)
    lim_max = lim.xmax().to(u.pc)
    lim.close()
    del lim

    # Second loop goes through each file in that set and plots it to a figure.
    for f in InputValues().time_dicts[files]:
        level = Plotting3d(f)
        level.plotsilo_3dslice(lim_min, lim_max, fig, gs, var1)

        level.close()
        del level

    # This now saves the figure that has had multiple levels plotted onto it.
    plt.savefig(imagefile, bbox_inches='tight', dpi=300)
    plt.close()

    del fig, gs, imagefile

    t1 = time.time()
    print("Running Time: {0}".format(t1 - t0))


# objgraph.show_growth()
# roots = objgraph.get_leaking_objects()
# print(len(roots))
# objgraph.show_most_common_types(objects=roots)
# objgraph.show_refs(roots[:3], refcounts=True, filename='roots.png')
# tracker.print_diff()

# To run the code with multiple processors:
# num_cores = mp.cpu_count()
# Parallel(n_jobs=num_cores)(delayed(plotter)(files) for files in time_dicts)

