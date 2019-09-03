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

# -------------- Set of libraries needed:
import numpy as np
# import astropy.constants as apc
# from astropy import units as u
# import math
# from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel, CustomKernel, Model2DKernel
# from astropy.modeling.models import Gaussian2D
# from astropy.io import fits
# import scipy
# from scipy import ndimage

from SiloHeader_data import OpenData

# --------------Class to access each 2D sub-domain and save density, temperature, etc. data:


class Read2dSiloData(OpenData):
    def variable(self, data):  # Retrieves the requested data from the silo file
        param = self.db.GetVar(data + "_data")  # Saves the selected data as a variable.
        param_dims = self.db.GetVar(data + "_dims")  # Saves the selected data's dimensions.
        param_dims = param_dims[::-1]  # Reshapes the dimensions into the correct orientation.
        param = np.array(param).reshape(param_dims)  # Puts the array into the correct format, i.e. (a,b).
        return param

    # Function that loops through all the selected data from each sub-domain and saves it to an empty 1D array:
    def parameter2d(self, data):
        array_param = []

        for n in range(self.nproc()):
            nn = "%04d" % (n,)  # Sub-domains are called rank_XXXX_domain_XXXX.
            self.db.SetDir("/rank_%s_domain_%s/" % (nn, nn))

            variable_data = self.variable(data)
            array_param.append(variable_data)  # Saves the selected data into the empty array.
        return array_param

    def reshaped_parameter2d(self, data):
        variable_array = np.zeros((self.ngrid()[1], self.ngrid()[0]))

        a = self.dom_size()['DomSize'][0]
        b = self.dom_size()['DomSize'][1]
        c = self.dom_size()['Ndom'][0]
        d = self.dom_size()['Ndom'][1]
        e = self.parameter2d(data)

        for jD in range(d):
            for iD in range(c):
                x0 = iD * a  # Sets the positions of each process array.
                y0 = jD * b
                x1 = x0 + a
                y1 = y0 + b

                domain = jD * c + iD
                variable_array[y0:y1, x0:x1] = e[domain]  # Saves all the values into the 2D image array

        del a
        del b
        del c
        del d
        del e

        return variable_array


class Read3dSiloData(OpenData):
    def variable(self, data):  # Retrieves the requested data from the silo file
        param = self.db.GetVar(data + "_data")  # Saves the selected data as a variable.
        param_dims = self.db.GetVar(data + "_dims")  # Saves the selected data's dimensions.
        param_dims = param_dims[::-1]  # Reshapes the dimensions into the correct orientation.
        param = np.array(param).reshape(param_dims)  # Puts the array into the correct format, i.e. (a,b).
        return param

    def parameter3d(self, data):
        array_param = []

        for n in range(self.nproc()):
            nn = "%04d" % (n,)  # Sub-domains are called rank_XXXX_domain_XXXX.
            self.db.SetDir("/rank_%s_domain_%s/" % (nn, nn))

            variable_data = self.variable(data)
            array_param.append(variable_data)  # Saves the selected data into the empty array.
        return array_param

    def reshaped_parameter3d(self, data):
        variable_array = np.zeros((self.ngrid()[2], self.ngrid()[1], self.ngrid()[0]))

        a = self.dom_size()['DomSize'][0]
        b = self.dom_size()['DomSize'][1]
        f = self.dom_size()['DomSize'][2]
        c = self.dom_size()['Ndom'][0]
        d = self.dom_size()['Ndom'][1]
        g = self.dom_size()['Ndom'][2]
        e = self.parameter3d(data)

        for kD in range(g):
            for jD in range(d):
                for iD in range(c):
                    x0 = iD * a  # Sets the positions of each process array.
                    y0 = jD * b
                    z0 = kD * f
                    x1 = x0 + a
                    y1 = y0 + b
                    z1 = z0 + f

                    domain = kD * d * d + jD * c + iD
                    variable_array[z0:z1, y0:y1, x0:x1] = e[domain]  # Saves all the values into the 2D image array

        del a
        del b
        del c
        del d
        del e
        del g
        del f

        return variable_array

