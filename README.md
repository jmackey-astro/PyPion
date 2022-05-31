# Python Library for PION
Created - 26|08|2019,
Modified - 31|05|2022

Welcome to the Python Library that does post-processing on the .Silo simulation data files outputted from PION. The library works on nested-grid and uniform-grid, and 3D and 2D Silo files.

The relevant files can be found in src/pypion. You can use PlotData.py as a test script, it can call a 2D and a 3D function from PlottingClasses.py. See comments insted the python script on how to use it and change for your needs.

# Install via pip:

You can now install PyPion via PyPi:

* Just run `python3 -m pip install pypion`

In some cases you may need to use `sudo` or install pip with: `sudo apt install python3-pip`

*Note: This is still in testing phase so errors may occur when installing. Please report problems to green@cp.dias.ie

# Install Silo package for python3:

The install script for Silo can be found at /src/silo/ in this repository. Please download it to ~/.local/silo on your computer and install it there using `bash install_silo.sh` 

# Info on how to get started

* Example usage can be found at: [https://www.pion.ie/docs/python.html](https://www.pion.ie/docs/python.html)

* Check out the [Wiki](https://git.dias.ie/compastro/pion_python/-/wikis/home) for more detailed usage.

* How to use the PyPion docker image at: [https://www.pion.ie/docs/python.html](https://www.pion.ie/docs/python.html)


# Developers/Maintainers:

* Samuel Green - green@cp.dias.ie
* Jonathan Mackey

