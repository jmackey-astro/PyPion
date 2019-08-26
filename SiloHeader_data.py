# Author: Sam Green, Created: 18-10-17
# Script opens up the header of silo files and saves some variables.

# New comments:
# - 2016-12.14 JM: ReadData has a CLOSE() function that calls the OpenSilo
#   CLOSE() function.  Xmin/Xmax/time etc. all now new directory to /header.
# - 15-02-2019 SG: Added new functions level_max and level_min and nlevels to help with nested grid python.

# -------------- Set of libraries needed:
import Silo
import numpy as np
from astropy import units as u

# -------------- Class to open Silo files and save variables.


class OpenData:
    def __init__(self, files):  # This will open the .SILO file and enters the 'header' directory.
        self.db = Silo.Open(files)
        self.db.SetDir('/header')

    def close(self):  # To close all the variables after use.
        self.db.Close()

    def header_info(self):  # Returns the header information:
        header = self.db.GetToc()
        return header

    # Xmax variable contains the max value for the x- and y-axis of the grid.
    def xmax(self):  # Function that returns the Xmax variable once called:
        self.db.SetDir('/header')
        xmax = self.db.GetVar("Xmax") * u.cm
        # xmax = xmax.to(u.pc)  # Converts the grid to parsec units.
        return xmax

    # Xmin variable contains the min value for the x- and y-axis of the grid.
    def xmin(self):
        self.db.SetDir('/header')
        xmin = self.db.GetVar("Xmin") * u.cm
        # xmin = xmin.to(u.pc)  # Converts the grid to parsec units.
        return xmin

    # level_xmax variable contains the max value for the x- and y-axis of each level.
    def level_max(self):  # Function that returns the level_xmax variable once called:
        self.db.SetDir('/header')
        level_max = self.db.GetVar("level_xmax") * u.cm
        return level_max

    # level_xmin variable contains the min value for the x- and y-axis of each level.
    def level_min(self):
        self.db.SetDir('/header')
        level_min = self.db.GetVar("level_xmin") * u.cm
        return level_min

    # nlevels variable contains the number of levels in the simulation.
    def nlevels(self):  # Function that returns the grid levels variable once called:
        self.db.SetDir('header')
        level = self.db.GetVar("grid_nlevels")
        return level

    # sim_time variable contains the simulation time.
    def sim_time(self):
        self.db.SetDir('/header')
        simtime = self.db.GetVar("t_sim") * u.second
        simtime = simtime.to(u.Myr)  # Converts time to Mega-years
        return simtime

    # ngrid variable contains the size of the grid.
    def ngrid(self):
        self.db.SetDir('/header')
        ngrid = self.db.GetVar("NGrid")  # Gets the size of the image grid
        return ngrid

    # nproc variable contains the number of domains in the silo file.
    def nproc(self):
        self.db.SetDir('/header')
        nproc = self.db.GetVar("MPI_nproc")  # Gets the number of MPI processes (== number of domains)
        return nproc

    # dom_size contains the size of the domain in the grid.
    def dom_size(self):  # Function to set up the domain size of the grid
        ndom = np.array([1, 1, 1])
        cells = np.copy(self.ngrid())
        i = 1
        while i < self.nproc():
            axis = np.argmax(cells)
            i *= 2
            ndom[axis] *= 2
            cells[axis] /= 2
        domsize = self.ngrid() / ndom
        ndom1 = ndom
        return {'DomSize': domsize, 'Ndom': ndom1}
