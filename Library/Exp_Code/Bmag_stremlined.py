# Author: Sam Green, Created: 12-01-21
# Email: green@cp.dias.ie

#############################################################
# -------------- Set of libraries needed:
from ReadData import ReadData

# Probably don't need most of these, they're just here for debugging purposes.
import matplotlib
# Using this to stop matplotlib from using a $DISPLAY environment variable.
# i.e. This now works over ssh.
matplotlib.use('Agg')
from matplotlib.colorbar import Colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

import numpy as np
from astropy import units as u
import math

import warnings
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

plt.rc('font', **{'size': 12})
# plt.rc('lines', linewidth=2)
plt.rc('font', weight='bold')


class Plotting3d(ReadData):
    def Bmag(self, param, Fig, var1):

        var = var1
        fig = Fig

        #Import the data you want. This function is currently using the B-field data
        # but can easily be changed to use the velocity data.
        magx = self.get_3Darray('MagneticFieldX')['data']
        magy = self.get_3Darray('MagneticFieldY')['data']
        magz = self.get_3Darray('MagneticFieldZ')['data']

        lim_max = (self.get_3Darray(param)['max_extents'] * u.cm).to(u.pc)
        lim_min = (self.get_3Darray(param)['min_extents'] * u.cm).to(u.pc)
        sim_time = self.get_3Darray(param)['sim_time'].to(u.Myr)
        ngrid = self.ngrid()

        # Setting the extents of the plot from the lim_max and lim_min params.
        Xmax = lim_max[0][1].value
        Ymax = lim_max[0][2].value
        Xmin = lim_min[0][1].value
        Ymin = lim_min[0][2].value

        # This variable specifies the speration between streamlines.
        # the bigger the number, the more streamlines on the plot, vice versa.
        sep = (Xmax - Xmin) / 256
        
        # Here you specify what slice from the 3d data you want to use.
        # I.e. if you have a 256^3 sim and you want the middle slice then you use 128.
        grid = 128
        
        # Calculate the magnitude of the 3D data.
        mag_3d = np.sqrt(np.square(magx) + np.square(magy) + np.square(magz))

        # Select the slice to plot from the 3D data and calculate the log.
        # Also enter data into numpy array for ease.
        # The sim I used this code for had 2 nested grids and I needed to use 
        # seperate slice positions for each grid to make the correct image.
        log_magdz = np.log10(mag_3d[0][:, :, grid])
        log_magdz_1 = np.log10(mag_3d[1][:, :, grid])

        # Select slice from data for streamlines.
        mag = mag_3d[0][:, :, grid]
        magy_2d = magy[0][:, :, grid]
        magz_2d = magz[0][:, :, grid]
        magx_2d = magx[0][:, :, grid]
        # Select which axis you want to use.
        Xnew = magy_2d
        Ynew = magz_2d
        Vnew = mag

        # Gives the correct dimensions to the x and y-axis
        # Sets up the correct postion of each vector arrow and it's correct starting point
        Xpos, Ypos = np.meshgrid(np.arange(Xmin, Xmax, sep), np.arange(Ymin, Ymax, sep))

        # Convert the lists to numpy arrays fro the streamplot function
        XPOS = np.array(Xpos)
        YPOS = np.array(Ypos)
        VELx = np.array(Xnew)
        VELy = np.array(Ynew)
        speed = np.array(Vnew)

        # Generate a variable linewidth for the streamlines, using the Velocity-Magnitude data
        lw = lw = 2 * speed / speed.max()
        
        #plot 2d image
        ax1 = fig.add_subplot(1, 1, 1)

        # ax1.set_title('    Time = %5.5f Myr' % sim_time.value)

        ax1.set_xlim(lim_min[0][1].value, lim_max[0][1].value)
        ax1.set_ylim(lim_min[0][2].value, lim_max[0][2].value)
        

        # Need to plot each slice on top of each other.
        im1 = ax1.imshow(log_magdz, interpolation='nearest', cmap='plasma',
                        extent=[lim_min[0][1].value, lim_max[0][1].value, lim_min[0][2].value, lim_max[0][2].value],
                        origin='lower', vmax=-4.5, vmin=-7)

        #plt.xlim(-0.75, 0.75)
        #plt.ylim(-0.75, 0.75)

        im1 = ax1.imshow(log_magdz_1, interpolation='nearest', cmap='plasma',
                        extent=[lim_min[1][1].value, lim_max[1][1].value, lim_min[1][2].value, lim_max[1][2].value],
                        origin='lower', vmax=-4.5, vmin=-7)
        
        # plot streamlines
        Q = plt.streamplot(XPOS, YPOS, VELx, VELy, color = 'white', linewidth = lw, arrowstyle = '->', arrowsize = 1)

        txt1 = ax1.text(0.7, 0.92, r'$log(|B|/G)$', transform=ax1.transAxes)
        txt1 = ax1.text(0.05, 0.92, 'x = -0.3 pc', transform=ax1.transAxes)
        ax1.set_xlabel('y-axis (pc)')
        ax1.set_ylabel('z-axis (pc)')

        divider2 = make_axes_locatable(ax1)
        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
        cbar2 = plt.colorbar(im1, cax=cax2, ticks=MultipleLocator(.5), format="%.2f")

        return fig
