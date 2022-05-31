# Author: Sam Green, Created: 18-10-2017
# Script to allow the user to use commandline arguments when running scripts.

# Comments from previous script version (Silo_Modules.py):
# - 2016.06.15 JM: sorting of list, adding path to filenames,
#   relative position for annotation.  Also added command-line args.
# - 2016.06.16 JM: improved removal of non-requested files.

# New Comments:
# - 15-02-2019 SG: Added a new argument to choose the type of file being used. Should've added it over a year ago.
# - 15-02-2019 SG: Added method to group files by timestep (useful for nested grid sims).

# --------------set of libraries needed:
import argparse
from os import listdir
import sys
sys.path.append('.')


class InputValues:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Process some values.',
            usage='script.py <path-to-files> <base-filename> <image-path> <image-filename> <dimension>')
        parser.add_argument('Path', help='path to Silo files')
        parser.add_argument('file_base', help='base filename of silo files')
        parser.add_argument('img_path', help='path to save images into')
        parser.add_argument('img_base', help='base filename of image files')
        parser.add_argument('dimension', help='dimension of data (1d, 2d, 3d)')
        args = parser.parse_args()

        self.dimen = args.dimension

        # make sure that path has a trailing /
        file_path = args.Path
        if not file_path.endswith('/'):
            file_path += "/"

        self.img_path = args.img_path
        if not self.img_path.endswith('/'):
            self.img_path += "/"

        self.img_file = args.img_base

        # Reads all the files in a directory and then puts their names into a 1D array.
        data_files = [f for f in listdir(file_path)]  # if isfile(join(file_path, f))]

        # Sort list by timestep
        data_files.sort()

        search = ".silo"

        # Remove non-silo files, and non-requested files.
        # search = None
        # if args.file_type == 'silo':
        #    search = ".silo"
        # elif args.file_type == 'vtk':
        #    search = ".vtk"
        # elif args.file_type == 'fits':
        #    search = ".fits"

        data_files = [f for f in data_files if search in f]

        search = args.file_base
        data_files = [f for f in data_files if search in f]

        # Want to include only the primary data files, with "_0000.*.silo"
        #print(len(data_files))
        search = "_0001."
        data_files = [f for f in data_files if search not in f]
        search = "_0002."
        data_files = [f for f in data_files if search not in f]
        search = "_0003."
        data_files = [f for f in data_files if search not in f]
        search = "_0004."
        data_files = [f for f in data_files if search not in f]
        #print(len(data_files))

        # Add path to files
        self.data_files_with_path = [f.replace(f, file_path+f) for f in data_files]

        #  Establish dictionaries to group files by timestep (i.e. such that all layers in a given timestep are in
        #  a list assigned to a key with that timestep)
        self.time_dicts = {}
        for index in range(len(self.data_files_with_path)):

            filename = self.data_files_with_path[index]
            time = filename[len(filename) - 13: len(filename) - 5]
            try:
                self.time_dicts[time].append(filename)
            except:
                self.time_dicts[time] = [filename]
