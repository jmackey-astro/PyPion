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
plt.rc('font', weight='normal')


class Plotting3d(ReadData):
    def Vmag(self, param, Fig, var1):
        vlim = 2500  # max velocity on colourbar in km/s
        lev = 2  # level to plot streamlines for

        var = var1
        fig = Fig

        #Import the data you want. This function is currently using the B-field data
        # but can easily be changed to use the velocity data.
        velx = self.get_3Darray('VelocityX')['data']
        vely = self.get_3Darray('VelocityY')['data']
        velz = self.get_3Darray('VelocityZ')['data']

        lim_max = (self.get_3Darray(param)['max_extents'] * u.cm)/1.e13
        lim_min = (self.get_3Darray(param)['min_extents'] * u.cm)/1.e13
        sim_time = self.get_3Darray(param)['sim_time'].to(u.d)
        ngrid = self.ngrid()
        # Here you specify what slice from the 3d data you want to use.
        # I.e. if you have a 256^3 sim and you want the middle slice then you use 128.
        grid = ngrid[0]/2
        
        # Calculate the velnitude of the 3D data.
        vel_3d = np.sqrt(np.square(velx) + np.square(vely) + np.square(velz))

        # Select the slice to plot from the 3D data and calculate the log.
        # Also enter data into numpy array for ease.
        # This is for 4 nested grids.
        veldz_0 = np.array(vel_3d[0][grid, :, :])*1.0e-5
        veldz_1 = np.array(vel_3d[1][grid, :, :])*1.0e-5
        veldz_2 = np.array(vel_3d[2][grid, :, :])*1.0e-5
        veldz_3 = np.array(vel_3d[3][grid, :, :])*1.0e-5

        ##########################################################
        # Setting the extents of the streamline plot and the data
        Xmax = lim_max[lev][0].value
        Ymax = lim_max[lev][1].value
        Xmin = lim_min[lev][0].value
        Ymin = lim_min[lev][1].value
        print(Xmin,Xmax,Ymin,Ymax)

        # This variable specifies the speration between streamlines.
        # the bigger the number, the more streamlines on the plot, vice versa.
        sep = (Xmax - Xmin) / 80
        
        # Select slice in X-Y plane from data for streamlines.
        vel = np.array(vel_3d[lev][grid, :, :])*1.0e-5
        vely_2d = np.array(vely[lev][grid, :, :])*1.0e-5
        velz_2d = np.array(velz[lev][grid, :, :])*1.0e-5
        velx_2d = np.array(velx[lev][grid, :, :])*1.0e-5
        # Select which axis you want to use.
        Xnew = velx_2d
        Ynew = vely_2d
        Vnew = vel

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
        lw = 1.5 * speed / vlim
        ##########################################################
        
        #plot 2d image
        ax1 = fig.add_subplot(1, 1, 1)

        # ax1.set_title('    Time = %5.5f Myr' % sim_time.value)

        ax1.set_xlim(lim_min[0][0].value, lim_max[0][0].value)
        ax1.set_ylim(lim_min[0][1].value, lim_max[0][1].value)
        ax1.set_xlim(-5,5)
        ax1.set_ylim(-5,5)
        

        ##########################################################
        # plot the image with velocity magnitude
        # Need to plot each slice on top of each other.
        im0 = ax1.imshow(veldz_0, interpolation='nearest', cmap='plasma',
                        extent=[lim_min[0][0].value, lim_max[0][0].value, lim_min[0][1].value, lim_max[0][1].value],
                        origin='lower', vmax=vlim, vmin=0)

        im1 = ax1.imshow(veldz_1, interpolation='nearest', cmap='plasma',
                        extent=[lim_min[1][0].value, lim_max[1][0].value, lim_min[1][1].value, lim_max[1][1].value],
                        origin='lower', vmax=vlim, vmin=0)
        im2 = ax1.imshow(veldz_2, interpolation='nearest', cmap='plasma',
                        extent=[lim_min[2][0].value, lim_max[2][0].value, lim_min[2][1].value, lim_max[2][1].value],
                        origin='lower', vmax=vlim, vmin=0)
        
        im3 = ax1.imshow(veldz_3, interpolation='nearest', cmap='plasma',
                        extent=[lim_min[3][0].value, lim_max[3][0].value, lim_min[3][1].value, lim_max[3][1].value],
                        origin='lower', vmax=vlim, vmin=0)
        ##########################################################
        
        
        ##########################################################
        # plot streamlines for level "lev"
        Q = plt.streamplot(XPOS, YPOS, VELx, VELy, color='white', linewidth = lw, arrowstyle = '->', arrowsize = 1)
        ##########################################################

        txt1 = ax1.text(0.65, 0.92, r'$log(|v|/ km/s)$', transform=ax1.transAxes, color='white')
        txt1 = ax1.text(0.05, 0.92, 'z = 0.0', transform=ax1.transAxes, color='white')
        ax1.set_xlabel('x-axis ($10^{13}$ cm)')
        ax1.set_ylabel('y-axis ($10^{13}$ cm)')

        divider2 = make_axes_locatable(ax1)
        cax2 = divider2.append_axes("right", size="5%", pad=0.05)
        cbar2 = plt.colorbar(im1, cax=cax2, ticks=MultipleLocator(500), format="%.2f")

        del  velx, vely, velz, lim_max,  lim_min,  sim_time,  ngrid
        del Xmax, Ymax, Xmin, Ymin, sep, vel_3d, veldz_0, veldz_1, veldz_2,
        del veldz_3, vel, vely_2d, velz_2d, velx_2d, Xnew, Ynew, Vnew
        return fig
