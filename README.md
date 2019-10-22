# Python Library for 2D and 3D Nested Grid PION
Created - 26|08|2019,
Modified - 04|09|2019,
Owner - Samuel Green,
Email - green@cp.dias.ie

Welcome to the the Python Library that does post-processing on the .Silo simulation data files outputted from PION. The libray works on nested-grid and uniform-grid, and 3D and 2D Silo files.


- - -
1. [How to Run](#run)
2. [Project Aims](#aims)
3. [List of python packages that are needed](#package)
4. [How to use the Library](#how)

- - -
<a name="run"></a>

## 1\. How to Run

The best way to run this library at the moment is to use the bash script 'list-files'. To run:
```sh
$ bash list-files 'file_location' 'image_destination'
```

04-09-19: This script has only been tested on Linux (working on other OS)

<a name="aims"></a>

## 2\. Project Aims

* We want a library of python scripts that can read in Silo files and use the data for post-processing.
* Make the code easy to use and understand for anyone that needs it.
* To have all the modules needed to analyse all the variables saved in the Silo file.
* For it to be modular, so you can pick and choose what to run and not have to run code that is unnecessary.
* For it to work!

<a name="package"></a>

## 3\. List of python packages that are needed:
To be able to use all the features of this library you will need to have the following python 
modules installed on your system. Obviously you don't need all of these if you only need parts of this library.

* Silo
* numpy
* astropy
* matplotlib
* tk (python-tk)

All of these modules can be installed through _pip_ on Linux (might have to install pip first using 'apt install'):

```sh
$ pip install 'python-module'
```
Edit 13-09-19 -> You have to install the Silo package using:

```sh
$ sudo apt install python-silo
```
All of these modules should be available on windows through whatever python ide you use (i.e. Anaconda). Not sure how it's done on a Mac, need to ask a Mac user...

<a name="how"></a>

## 4\. How to use the Library

At the moment the main scripts in this library are:

* argparse_command.py - Saves the options entered into the command line when the python script is run. (Currently being taken out, 04-09-19)
* SiloHeader\_data.py - Which opens the silo file and saves all of the important header variables (eg. sim_time, xmax, xmin, etc.).
* ReadData.py - Opens the directory in the silo (or vtk, or fits) file and saves the requested variable data (eg. density, temp, etc.).
* Plotting_Classes.py - Sets up the plotting function and the figure. 
* plot\_nestedgrid\_3dslice.py -  Uses functions from the previous scripts to create 2d slice figures from 3D nested\_grid_pion data.
* plot\_nestedgrid.py - Uses functions from the previous scripts to create figures from 2D nested\_grid_pion data. 


This description is a work-in-progress...
