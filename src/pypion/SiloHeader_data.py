# Author: Sam Green, Created: 18-10-17
# Script opens up the header of silo files and saves some variables.

# New comments:
# - 2016-12.14 JM: ReadData has a CLOSE() function that calls the OpenSilo
#   CLOSE() function.  Xmin/Xmax/time etc. all now new directory to /header.
# - 15-02-2019 SG: Added new functions level_max and level_min and nlevels to help with nested grid python.
# - 22-04-2020 SG: Added all functions that use the init function in the OpenData class.
# - 22-04-2020 SG: Removed hard-coded cm in xmax and xmin.

# -------------- Set of libraries needed:
import sys

# Point to silo installation folder.
# I recommend installed silo here at ~/.local/silo/
from pathlib import Path
home = str(Path.home())
sys.path.insert(0,home+"/.local/silo/")
import Silo

import numpy as np
from astropy import units as u
from os import path as ospath

# -------------- Class to open Silo files and save variables.


class OpenData:
    # This will open the .SILO file and enters the 'header' directory.
    # If using multiple levels, this will open the level 0 file.
    def __init__(self, files):
        self.data = files
        self.db = Silo.Open(files[0])
        self.db.SetDir('/header')
        self.OpenFile = files[0]

    # Call this function with a for loop to cycle through each level file.
    def open(self, level):
        self.db.Close()
        self.db = Silo.Open(self.data[level])
        self.db.SetDir('/header')
        self.OpenFile = self.data[level]

    # To close all the variables after use.
    def close(self):
        self.db.Close()
        self.OpenFile = ""

    # Returns the header information:
    def header_info(self):
        header = self.db.GetToc()
        return header

    # Xmax variable contains the max value for every axis of the grid.
    def xmax(self):
        self.db.SetDir('/header')
        xmax = self.db.GetVar("Xmax")
        return xmax

    # Xmin variable contains the min value for every axis of the grid.
    def xmin(self):
        self.db.SetDir('/header')
        xmin = self.db.GetVar("Xmin")
        return xmin

    # level_xmax variable contains the max value for every axes of each level.
    def level_max(self):
        self.db.SetDir('/header')
        level_max = self.db.GetVar("level_xmax")
        return level_max

    # level_xmin variable contains the min value for every axis of each level.
    def level_min(self):
        self.db.SetDir('/header')
        level_min = self.db.GetVar("level_xmin")
        return level_min

    # nlevels variable contains the number of levels in the simulation.
    def nlevels(self):
        self.db.SetDir('/header')
        level = self.db.GetVar("grid_nlevels")
        return level

    # ndim variable: how many spatial dimensions on grid
    def ndim(self):
        self.db.SetDir('/header')
        dim = self.db.GetVar("gridndim")
        return dim

    # cycle is the simulation timestep on the finest level.
    def cycle(self):
        self.db.SetDir('/')
        cycle = self.db.GetVar("cycle")
        return cycle

    # sim_time variable contains the simulation time.
    def sim_time(self):
        self.db.SetDir('/header')
        simtime = self.db.GetVar("t_sim") * u.second
        return simtime

    # ngrid variable contains the size of the grid.
    def ngrid(self):
        self.db.SetDir('/header')
        ngrid = self.db.GetVar("NGrid")
        ngrid = [int(_) for _ in ngrid]
        return ngrid

    # nproc variable contains the number of domains in the silo file.
    def nproc(self):
        self.db.SetDir('/header')
        nproc = self.db.GetVar("MPI_nproc")
        return nproc

    # dom_size contains the size of the domain in the grid.
    def dom_size(self):
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
        else:
          param_dims = param_dims[::-1]
        # Puts the array into the correct format, i.e. (a,b).
        param = np.array(param).reshape(param_dims)
        # return the data, and the number of dimensions
        return param, len(param_dims)

    # Opens up the sub-domains and saves the data specified in "variable".
    # "data" is a string corresponding to a scalar variable, e.g. "Density"
    def parameter(self, data):
        array_param = []  # data for each sub-domain
        exten_param = []  # min-extents of each sub-domain
        self.db.SetDir("/")
        extents = self.db.GetVar("MultiMesh_extents")
        i = 0  # index in extents
        pp = self.db.GetVar(data+"_varnames")
        newfile=""
        nowfile=""
        newpar=""
        origfile=""
        #for n in paths:
        for n in pp:
            if n=="":
              continue # 1st and last elements are empty
            # if "n" contains a ":" then the domain is in a different file.
            # We need to close the current file and open another one.
            loc=n.find(":")
            if loc != -1:
              if (newfile==""):
                origfile=self.OpenFile
                nowfile=self.OpenFile
              else:
                nowfile=newfile
              # split newfile into filename + location
              newfile= n[:loc]
              newpar = n[loc+1:]
              # append full path to new filename
              newfile = ospath.dirname(self.OpenFile) + "/" + newfile
              # only close old / open new file if they are not the same file
              if (newfile != nowfile):
                self.db.Close()
                self.db = Silo.Open(newfile)
                self.db.SetDir("/")
                nowfile=newfile
            else:
              # nothing do do because domain is in my file.
              newpar = n
              nowfile= self.OpenFile
            #print(nowfile,newpar)
            d,v = ospath.split(newpar)
            rank = d[6:10]  # save rank so i can get extents.
            self.db.SetDir(d)
            variable_data, idim = self.variable(v)
            array_param.append(variable_data)
            min_extents = extents[i:i+idim]
            i = i + 2*idim
            exten_param.append(min_extents)
        # re-open original file that was being used for the 1st sub-domain.
        self.db.Close()
        self.db = Silo.Open(self.OpenFile)
        self.db.SetDir("/")
        return np.array(array_param), np.array(exten_param)

