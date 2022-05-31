# Author: Sam Green, Created: 18-10-2017
# Script sets up a set of core classes/functions to read data from silo data files.

# Comments from previous script version (Silo_Modules.py):
# -2016-06-10 SG: Opens MultiMesh silo files and can plot 2D structures.
# -2016-09-29 SG: PARAMETER and PARAMETER_MAIN Function added to save density and temperature data depending on
# what variable name is called,
# no longer separate functions to calculate density and temperature data.
# -2016-10-07 SG: PARAMETER Function modified to allow velocity to be chosen as the variable name.
# -2016-12.14 JM: self.Param and self.BigArray_Param no longer class
#   member data.
# -2017-02-23 SG: Added comments.
# -2017-03-10 SG: Added Class to read in VTK data produced by PION 2D projection code and calculate data to
# plot the column density and nebular line emission.

# New comments:
# - 2017-10-28 SG: Rewritten some functions to remove errors, increase efficiency.
# - 2017-10-29 SG: Bugs in reshaped_parameter2d function fixed.
# - 2018-03-26 SG: Added class to plot slices of the 3D data.
# - 22-04-2020 SG: Moved the parameter and variable function to SiloHeader_data.py
# - 22-04-2020 SG: Both functions now open all the files associated with each timestep and saves the data arrays for each of them.

# This is to make the scripts work with nested-grids.
# Works in 1D, 2D, and 3D. Should also work with non-nestedgrid too (i.e. 1 level of data).

# -------------- Set of libraries needed:
import numpy as np
from .SiloHeader_data import OpenData

# --------------Class to access each 2D sub-domain and save density, temperature, etc. data:


class ReadData(OpenData):

    # Function to read 1D data from multiple nested grid levels from
    # a Silo file, and return the data as numpy arrays with extents.
    def get_1Darray(self, param):

        level = self.nlevels()
        arr =  [[None]] * level
        level_max = [[None]] * level
        level_min = [[None]] * level
        sim_time = self.sim_time().value

        i = 0
        for file in self.data:
            self.open(i)

            variable_array = np.zeros((self.ngrid()[0]))

            domain=self.dom_size()
            a = domain['DomSize'][0]
            c = domain['Ndom'][0]
            e, min_ext = self.parameter(param)
            level_max[i] = self.level_max()
            level_min[i] = self.level_min()
            dx = (level_max[i][0]-level_min[i][0])/c

            for idom in range(c):
              # Sets the positions of each process array.
              ix = int((min_ext[idom][0] - level_min[i][0])*1.01/dx)
              x0 = ix * a
              x1 = x0 + a
              # Saves all the values into the 1D image array
              variable_array[x0:x1] = e[idom]

            arr[i] = variable_array
            i += 1

            del a
            del c
            del e
            del variable_array

        return {'data': arr, 'max_extents': level_max, 'min_extents': level_min, 'sim_time': sim_time}

    # Function to read 2D data from multiple nested grid levels from
    # a Silo file, and return the data as numpy arrays with extents.
    def get_2Darray(self, param):

        level = self.nlevels()
        arr =  [[None]] * level
        level_max = [[None]] * level
        level_min = [[None]] * level
        sim_time = self.sim_time().value

        i = 0
        for file in self.data:
            self.open(i)

            variable_array = np.zeros((self.ngrid()[1], self.ngrid()[0]))

            domain=self.dom_size()
            a = domain['DomSize'][0]
            b = domain['DomSize'][1]
            c = domain['Ndom'][0]
            d = domain['Ndom'][1]
            e, min_ext = self.parameter(param)
            level_max[i] = self.level_max()
            level_min[i] = self.level_min()
            dx = (level_max[i][0]-level_min[i][0])/c
            dy = (level_max[i][1]-level_min[i][1])/d
            ndom = c * d

            for idom in range(ndom):
              # get the (ix,iy,iz) location of sub-domain from its extents
              ix = int((min_ext[idom][0] - level_min[i][0])*1.01/dx)
              iy = int((min_ext[idom][1] - level_min[i][1])*1.01/dy)
              x0 = ix * a
              y0 = iy * b
              x1 = x0 + a
              y1 = y0 + b
              variable_array[y0:y1, x0:x1] = e[idom]
              
            arr[i] = variable_array
            i += 1

            del a
            del b
            del c
            del d
            del e
            del variable_array

        return {'data': arr, 'max_extents': level_max, 'min_extents': level_min, 'sim_time': sim_time}

    # Function to read 3D data from multiple nested grid levels from
    # a Silo file, and return the data as numpy arrays with extents.
    def get_3Darray(self, param):

        level = self.nlevels()
        arr =  [[None]] * level
        level_max = [[None]] * level
        level_min = [[None]] * level
        sim_time = self.sim_time()

        i = 0
        for file in self.data:
            self.open(i)

            variable_array = np.zeros((self.ngrid()[2], self.ngrid()[1], self.ngrid()[0]))
            domain=self.dom_size()
            a = domain['DomSize'][0]
            b = domain['DomSize'][1]
            f = domain['DomSize'][2]
            c = domain['Ndom'][0]
            d = domain['Ndom'][1]
            g = domain['Ndom'][2]
            e, min_ext = self.parameter(param)
            level_max[i] = self.level_max()
            level_min[i] = self.level_min()
            dx = (level_max[i][0]-level_min[i][0])/c
            dy = (level_max[i][1]-level_min[i][1])/d
            dz = (level_max[i][2]-level_min[i][2])/g
            ndom = c * d * g

            for idom in range(ndom):
              # get the (ix,iy,iz) location of sub-domain from its extents
              ix = int((min_ext[idom][0] - level_min[i][0])*1.01/dx)
              iy = int((min_ext[idom][1] - level_min[i][1])*1.01/dy)
              iz = int((min_ext[idom][2] - level_min[i][2])*1.01/dz)
              x0 = ix * a
              y0 = iy * b
              z0 = iz * f
              x1 = x0 + a
              y1 = y0 + b
              z1 = z0 + f
              variable_array[z0:z1, y0:y1, x0:x1] = e[idom]
              
            arr[i] = variable_array
            i += 1

            del a
            del b
            del c
            del d
            del e
            del g
            del f
            del variable_array

        return {'data': arr, 'max_extents': level_max, 'min_extents': level_min, 'sim_time': sim_time}
