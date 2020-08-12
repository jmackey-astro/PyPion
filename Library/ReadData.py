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
from SiloHeader_data import OpenData

# --------------Class to access each 2D sub-domain and save density, temperature, etc. data:


class ReadData(OpenData):

    def get_1Darray(self, param):

        level = self.nlevels()
        arr =  [[None]] * level
        level_max = [[None]] * level
        level_min = [[None]] * level
        sim_time = self.sim_time().value

        i = 0
        for file in self.data:
            #print(i,file)
            #opendata = OpenData(file)
            self.open(i)

            variable_array = np.zeros((self.ngrid()[0]))

            a = self.dom_size()['DomSize'][0]
            c = self.dom_size()['Ndom'][0]
            e = self.parameter(param)

            for iD in range(c):
              # Sets the positions of each process array.
              x0 = iD * a
              x1 = x0 + a

              domain = iD
              # Saves all the values into the 1D image array
              variable_array[x0:x1] = e[domain]
              level_max[i] = self.level_max()
              level_min[i] = self.level_min()

            arr[i] = variable_array
            i += 1

            del a
            del c
            del e
            del variable_array

        return {'data': arr, 'max_extents': level_max, 'min_extents': level_min, 'sim_time': sim_time}

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

            a = self.dom_size()['DomSize'][0]
            b = self.dom_size()['DomSize'][1]
            c = self.dom_size()['Ndom'][0]
            d = self.dom_size()['Ndom'][1]
            e = self.parameter(param)

            for jD in range(d):
                for iD in range(c):
                    # Sets the positions of each process array.
                    x0 = iD * a
                    y0 = jD * b
                    x1 = x0 + a
                    y1 = y0 + b

                    domain = jD * c + iD
                    # Saves all the values into the 2D image array
                    variable_array[y0:y1, x0:x1] = e[domain]
                    level_max[i] = self.level_max()
                    level_min[i] = self.level_min()

            arr[i] = variable_array
            i += 1

            del a
            del b
            del c
            del d
            del e
            del variable_array

        return {'data': arr, 'max_extents': level_max, 'min_extents': level_min, 'sim_time': sim_time}

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
            #print("domain",domain)
            e = self.parameter(param)
            #print("parameter",len(e))

            for kD in range(g):
                for jD in range(d):
                    for iD in range(c):
                        # Sets the positions of each process array.
                        x0 = iD * a
                        y0 = jD * b
                        z0 = kD * f
                        x1 = x0 + a
                        y1 = y0 + b
                        z1 = z0 + f

                        domain = kD * d * c + jD * c + iD
                        # Saves all the values into the 3D image array
                        #print(domain,z0,y0,x0)
                        variable_array[z0:z1, y0:y1, x0:x1] = e[domain]

            arr[i] = variable_array
            level_max[i] = self.level_max()
            level_min[i] = self.level_min()
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
