# Author: Sam Green, Created: 11-10-20

# Script to call all PyPion classes to plot 2d + 3d Silo data.

from .Plotting_Classes import Plotting2d, Plotting3d
from .argparse_command import InputValues
#-------------------------------
import matplotlib
# Using this to stop matplotlib from using a $DISPLAY environment variable.
# i.e. This now works over ssh without the need for a x-server.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#-------------------------------
# The way I plot 3 arrays onto the same grid throws up a warning.
# This just surpresses the warning.
import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
#-------------------------------
line = InputValues()
time_dicts = line.time_dicts
dimen = line.dimen
#-------------------------------

for files in time_dicts:

	arr = time_dicts[files]

	var1 = ["Density", -22, -27, "viridis", 'y', 63]
	# var1 = ["Temperature", 8, 3, "inferno", 'log', 'y', 127]

	fig = plt.figure()

	a# = Plotting2d(arr).plotsilo_2d(var1[0], fig, var1)
	
	a = Plotting3d(arr).XZXYslice(var1[0], fig, var1)

	imagefile = "%s%s_%s.png" % (line.img_path, line.img_file, time_dicts[files][0][len(time_dicts[files][0]) - 13:len(time_dicts[files][0]) - 5])
	plt.savefig(imagefile, bbox_inches='tight', dpi=300)

