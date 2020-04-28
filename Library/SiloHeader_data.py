# Author: Sam Green, Created: 18-10-17
# Script opens up the header of silo files and saves some variables.

# New comments:
# - 2016-12.14 JM: ReadData has a CLOSE() function that calls the OpenSilo
#   CLOSE() function.  Xmin/Xmax/time etc. all now new directory to /header.
# - 15-02-2019 SG: Added new functions level_max and level_min and nlevels to help with nested grid python.
# - 22-04-2020 SG: Added all functions that use the init function in the OpenData class.
# - 22-04-2020 SG: Removed hard-coded cm in xmax and xmin.

# -------------- Set of libraries needed:
import Silo
import numpy as np
from astropy import units as u

# -------------- Class to open Silo files and save variables.


class OpenData:
    def __init__(self, files):  # This will open the .SILO file and enters the 'header' directory.
        #print(files)
        self.data = files
        self.db = Silo.Open(files[0])
        self.db.SetDir('/header')

    def open(self,level):
        self.db.Close()
        self.db = Silo.Open(self.data[level])
        self.db.SetDir('/header')
        #print("OPEN",self.level_max())

    def close(self):  # To close all the variables after use.
        self.db.Close()

    def header_info(self):  # Returns the header information:
        header = self.db.GetToc()
        return header

    # Xmax variable contains the max value for the x- and y-axis of the grid.
    def xmax(self):  # Function that returns the Xmax variable once called:
        self.db.SetDir('/header')
        xmax = self.db.GetVar("Xmax")
        return xmax

    # Xmin variable contains the min value for the x- and y-axis of the grid.
    def xmin(self):
        self.db.SetDir('/header')
        xmin = self.db.GetVar("Xmin")
        return xmin

    # level_xmax variable contains the max value for the x- and y-axis of each level.
    def level_max(self):  # Function that returns the level_xmax variable once called:
        self.db.SetDir('/header')
        level_max = self.db.GetVar("level_xmax")
        return level_max

    # level_xmin variable contains the min value for the x- and y-axis of each level.
    def level_min(self):
        self.db.SetDir('/header')
        level_min = self.db.GetVar("level_xmin")
        return level_min

    # nlevels variable contains the number of levels in the simulation.
    def nlevels(self):  # Function that returns the grid levels variable once called:
        self.db.SetDir('/header')
        level = self.db.GetVar("grid_nlevels")
        return level

    # ndim variable: how many spatial dimensions on grid
    def ndim(self):  # Function that returns the grid levels variable once called:
        self.db.SetDir('/header')
        dim = self.db.GetVar("gridndim")
        return dim

    # cycle is the simulation timestep on the finest level.
    def cycle(self):  # Function that returns the simulation cycle:
        self.db.SetDir('/')
        cycle = self.db.GetVar("cycle")
        return cycle

    # sim_time variable contains the simulation time.
    def sim_time(self):
        self.db.SetDir('/header')
        simtime = self.db.GetVar("t_sim") * u.second
        #simtime = simtime.to(u.Myr)  # Converts time to Mega-years
        return simtime

    # ngrid variable contains the size of the grid.
    def ngrid(self):
        self.db.SetDir('/header')
        ngrid = self.db.GetVar("NGrid")  # Gets the size of the image grid
        ngrid = [int(_) for _ in ngrid]
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
        domsize = [int(_) for _ in domsize]
        ndom1 = ndom
        ndom1 = [int(_) for _ in ndom1]
        return {'DomSize': domsize, 'Ndom': ndom1}

    # Retrieves the requested data from the silo file
    def variable(self, par):
        # Saves the selected data as a variable.
        param = self.db.GetVar(par + "_data")
        # Saves the selected data's dimensions.
        param_dims = self.db.GetVar(par + "_dims")
        # Reshapes the dimensions into the correct orientation.
        if (self.ndim()==1):
          param_dims = [param_dims]
        param_dims = param_dims[::-1]
        # Puts the array into the correct format, i.e. (a,b).
        param = np.array(param).reshape(param_dims)
        return param

    def parameter(self, data):
        array_param = []
        for n in range(self.nproc()):
            nn = "%04d" % (n,)  # Sub-domains are called rank_XXXX_domain_XXXX.
            self.db.SetDir("/rank_%s_domain_%s/" % (nn, nn))

            variable_data = self.variable(data)
            array_param.append(variable_data)  # Saves the selected data into the empty array.
        return np.array(array_param)
