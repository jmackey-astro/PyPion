__author__ = "Thomas Jones, Jonathan Mackey"

from ReadData import ReadData
import yt
import os
import glob
import numpy as np
import re
from derived_fields_yt import derived_fields

##########################################################################
def make_snapshots(data_path, file_base):
    """
    Function to make snapshots of silo files.

    Input: path to silo files, base of filename
    Output: list of snapshots - each snapshot is a list of silo files of different levels at the same time instant.

    Example:
    evolution = make_snapshots(data_path, file_base)
    
    """
    ########## Cataloging silo files ###############################
    # Get the list of silo files
    silo_files = os.listdir(data_path)
    #silo_files = sorted([f for f in silo_files if f.contains(file_base)]
    silo_files = sorted([f for f in silo_files if f.endswith('.silo')])
    filename=(silo_files[0]).replace('_level00_0000.00000000.silo','')
    filename=file_base
    
    print("Info from silo files:")
    print(f"Basename of silo files: {filename}")

    file_list = glob.glob(os.path.join(data_path, '*.silo'), recursive=False)

    level_list = []
    files = []

    for f in file_list:
        level = re.search('_level(.*)_', f)
        if level == None:
            pass
        else:
            level = level.group(1)
            if not level in level_list:
                level_list.append(level)
    level_list.sort()

    ########## Categorizing data files into levels #################
    if len(level_list) == 1: 
        print('Simulation Info: Single level')
        catalog = []
        files = sorted(glob.glob(filename + '_0000.*.silo'))
        catalog.append(files)
    else:
        print(f'Simulation Info: {len(level_list)} levels')
        catalog = []
        for i in range(len(level_list)):
            files = sorted(glob.glob(os.path.join(data_path, f"{filename}_level{level_list[i]}_0000.*.silo")))
            catalog.append(files)
            

    # Bundle silo files of different levels of same time instant into a snapshot.
    evolution = np.array(catalog).T
    print(f"Number of snapshots: {evolution.shape[0]}")

    return evolution
##########################################################################


##########################################################################
def get_ds(filename, quantities=["density"], **kwargs) -> yt.data_objects.static_output.Dataset:

    """
    Function to create a yt dataset from a silo file.
    Input: snapshot of silo file (list of filepaths for a single snapshot), list of quantities to be loaded
    Output: yt dataset

    Example:
    evolution = make_snapshots(data_path)
    file = evolution[0]
    ds = get_ds(file, quantities=["density", "temperature", "velocity", "magnetic_field"])

    """
    data = ReadData(filename)
    ndim = data.ndim()

    if (ndim==3):
      ds = get_ds3D(data, quantities, **kwargs)
    elif (ndim==2):
      ds = get_ds2D(data, quantities, **kwargs)
    elif (ndim==1):
      ds = get_ds1D(data, quantities, **kwargs)
    else:
      print("Bad ndim in silo file:", ndim)
      quit()

    return ds
##########################################################################


##########################################################################
def get_ds3D(data, quantities=["density"], **kwargs) -> yt.data_objects.static_output.Dataset:
    sim_time = data.sim_time()
    N_levels = data.nlevels()
    Dom_size = max(np.array(data.level_max()) - np.array(data.level_min()))
    rdata =  data.get_3Darray("Density")
    data_den = rdata['data']
    # get extents normalised to length 1, and with negative corner at the origin:
    lmin = np.array(rdata['min_extents']) #/Dom_size
    lmax = np.array(rdata['max_extents']) #/Dom_size
    # switch ordering or array because pypion returns ordering (z,y,x)
    N_grids = data.ngrid()[::-1]
    lmin = np.fliplr(lmin)
    lmax = np.fliplr(lmax)
    print(f"Simulation Info: {N_levels} levels, {N_grids} grids, {Dom_size} cm")
    

    # Arrays for quantities
    if quantities.__contains__("density"):
      data_den = rdata['data']
    if quantities.__contains__("logdensity"):
      data_lden = np.log10(rdata['data'])
    if quantities.__contains__("temperature"):
      data_temp = data.get_3Darray("Temperature")['data']
    if quantities.__contains__("pressure"):
      data_pressure = data.get_3Darray("Pressure")['data']
    if quantities.__contains__("NG_Mask"):
      data_ngmask = data.get_3Darray("NG_Mask")['data']
    if quantities.__contains__("windtracer"):
      data_windtr = data.get_3Darray("Tr000_WIND")['data']
    if quantities.__contains__("windtracer0"):
      data_windtr0 = data.get_3Darray("Tr000_WIND0")['data']
    if quantities.__contains__("windtracer1"):
      data_windtr1 = data.get_3Darray("Tr000_WIND1")['data']
    if quantities.__contains__("velocity"):
      data_velx = data.get_3Darray("VelocityX")['data']
      data_vely = data.get_3Darray("VelocityY")['data']
      data_velz = data.get_3Darray("VelocityZ")['data']
    if quantities.__contains__("magnetic_field"):
      data_Bx = data.get_3Darray("MagneticFieldX")['data']
      data_By = data.get_3Darray("MagneticFieldY")['data']
      data_Bz = data.get_3Darray("MagneticFieldZ")['data']
    if quantities.__contains__("xray_emission"):
      df = derived_fields()
      out = df.xray_em_interp3D(data)
      data_emx001 = out['emx_001']
      data_emx002 = out['emx_002']
      data_emx003 = out['emx_003']
      data_emx005 = out['emx_005']
      data_emx010 = out['emx_010']
      data_emx020 = out['emx_020']
      data_emx025 = out['emx_025']
      data_emx030 = out['emx_030']
      data_emx050 = out['emx_050']
      data_emx100 = out['emx_100']

    grid_data = np.array([dict(left_edge=lmin[n], 
                            right_edge=lmax[n], 
                            level=n, dimensions=N_grids) for n in range(N_levels)])

    bbox = np.array([lmin[0],lmax[0]]).T

    i = 0
    for g in grid_data:
        if quantities.__contains__("density"):
            g[("gas", "density")] = (data_den[i], "g/cm**3")
        if quantities.__contains__("logdensity"):
            g[("gas", "logdensity")] = (data_lden[i], "")
        if quantities.__contains__("temperature"):
            g[("gas", "temperature")] = (data_temp[i], "K")
        if quantities.__contains__("pressure"):
            g[("gas", "pressure")] = (data_pressure[i], "dyne/cm**2")
        if quantities.__contains__("NG_Mask"):
            g[("gas", "NG_Mask")] = (data_ngmask[i], "K")
        if quantities.__contains__("windtracer"):
            g[("gas", "windtracer")] = (data_windtr[i], "")
        if quantities.__contains__("velocity"):
            g[("gas", "velocity_x")] = (data_velx[i], "cm/s")
            g[("gas", "velocity_y")] = (data_vely[i], "cm/s")
            g[("gas", "velocity_z")] = (data_velz[i], "cm/s")
        if quantities.__contains__("magnetic_field"):
            g[("gas", "magnetic_field_x")] = (data_Bx[i], "G")
            g[("gas", "magnetic_field_y")] = (data_By[i], "G")
            g[("gas", "magnetic_field_z")] = (data_Bz[i], "G")
        if quantities.__contains__("xray_emission"):
            g[("gas", "xray_0.1")]  = (data_emx001[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_0.2")]  = (data_emx002[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_0.3")]  = (data_emx003[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_0.5")]  = (data_emx005[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_1.0")]  = (data_emx010[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_2.0")]  = (data_emx020[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_2.5")]  = (data_emx025[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_3.0")]  = (data_emx030[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_5.0")]  = (data_emx050[i], "erg/cm**3/s/arcmin**2")
            g[("gas", "xray_10.0")] = (data_emx100[i], "erg/cm**3/s/arcmin**2")
            print("adding emissivity to g,",i)

        i += 1

    data.close()

    ds = yt.load_amr_grids(grid_data, N_grids, length_unit=f"cm",
                           geometry=("cartesian"), axis_order=("z","y","x"),
                           sim_time=sim_time.value, bbox=bbox, time_unit="s",
                           unit_system="cgs",periodicity=(False, False, False))
    return ds
##########################################################################

##########################################################################
def get_ds2D(data, quantities=["density"], **kwargs) -> yt.data_objects.static_output.Dataset:
    sim_time = data.sim_time()
    N_levels = data.nlevels()
    Dom_size = max(np.array(data.level_max()) - np.array(data.level_min()))
    rdata =  data.get_2Darray("Density")
    data_den = rdata['data']
    # get extents in cgs units
    lmin = np.array(rdata['min_extents'])
    lmax = np.array(rdata['max_extents'])
    # switch ordering or array because pypion returns ordering (z,y,x)
    lmin[:, [0, 1]] = lmin[:, [1,0]] # swap z,r
    lmax[:, [0, 1]] = lmax[:, [1,0]] # swap z,r
    N_grids = np.array(data.ngrid())
    N_grids[[0, 1]]= N_grids[[1, 0]] # swap z,r
    lmin[:,2] = 0         # theta dir should have finite extent
    lmax[:,2] = 2.0*np.pi # theta dir should have finite extent
    #print("lmin",lmin)
    #print("lmax",lmax)
    print(f"Simulation Info: {N_levels} levels, {N_grids} grids, {Dom_size} cm")

    # Arrays for quantities
    if quantities.__contains__("density"):
      data_den = rdata['data']
    if quantities.__contains__("logdensity"):
      data_lden = np.log10(rdata['data'])
    if quantities.__contains__("temperature"):
      data_temp = data.get_2Darray("Temperature")['data']
    if quantities.__contains__("pressure"):
      data_pressure = data.get_2Darray("Pressure")['data']
    if quantities.__contains__("NG_Mask"):
      data_ngmask = data.get_2Darray("NG_Mask")['data']
    if quantities.__contains__("windtracer"):
      data_windtr = data.get_2Darray("Tr000_WIND")['data']
    if quantities.__contains__("windtracer0"):
      data_windtr0 = data.get_2Darray("Tr000_WIND0")['data']
    if quantities.__contains__("windtracer1"):
      data_windtr1 = data.get_2Darray("Tr000_WIND1")['data']
    if quantities.__contains__("velocity"):
      data_velx = data.get_2Darray("VelocityX")['data']
      data_vely = data.get_2Darray("VelocityY")['data']
      data_velz = data.get_2Darray("VelocityZ")['data']
    if quantities.__contains__("magnetic_field"):
      data_Bx = data.get_2Darray("MagneticFieldX")['data']
      data_By = data.get_2Darray("MagneticFieldY")['data']
      data_Bz = data.get_2Darray("MagneticFieldZ")['data']
    if quantities.__contains__("xray_emission"):
      df = derived_fields()
      out = df.xray_em_interp2D(data)
      data_emx001 = out['emx_001']
      data_emx002 = out['emx_002']
      data_emx003 = out['emx_003']
      data_emx005 = out['emx_005']
      data_emx010 = out['emx_010']
      data_emx020 = out['emx_020']
      data_emx025 = out['emx_025']
      data_emx030 = out['emx_030']
      data_emx050 = out['emx_050']
      data_emx100 = out['emx_100']

    grid_data = np.array([dict(left_edge=lmin[n], 
                        right_edge=lmax[n], 
                        level=n, dimensions=N_grids) for n in range(N_levels)])
    bbox = np.array([lmin[0],lmax[0]]).T

    i = 0
    for g in grid_data:
      if quantities.__contains__("density"):
        g[("gas", "density")] = (data_den[i].reshape((data_den[i].shape[0], data_den[i].shape[1], 1)), "g/cm**3")
        #  g[("gas", "density")] = (data_den[i], "g/cm**3")
      if quantities.__contains__("logdensity"):
        g[("gas", "logdensity")] = (data_lden[i].reshape((data_lden[i].shape[0], data_lden[i].shape[1], 1)), "")
      if quantities.__contains__("temperature"):
        g[("gas", "temperature")] = (data_temp[i].reshape((data_temp[i].shape[0], data_temp[i].shape[1], 1)), "K")
      if quantities.__contains__("pressure"):
        g[("gas", "pressure")] = (data_pressure[i].reshape((data_pressure[i].shape[0], data_pressure[i].shape[1], 1)), "dyne/cm**2")
      if quantities.__contains__("NG_Mask"):
        g[("gas", "NG_Mask")] = (data_ngmask[i].reshape((data_ngmask[i].shape[0], data_ngmask[i].shape[1], 1)), "K")
      if quantities.__contains__("windtracer"):
        g[("gas", "windtracer")] = (data_windtr[i].reshape((data_windtr[i].shape[0], data_windtr[i].shape[1], 1)), "")
      if quantities.__contains__("velocity"):
        g[("gas", "velocity_x")] = (data_velx[i].reshape((data_velx[i].shape[0], data_velx[i].shape[1], 1)), "cm/s")
        g[("gas", "velocity_y")] = (data_vely[i].reshape((data_vely[i].shape[0], data_vely[i].shape[1], 1)), "cm/s")
        g[("gas", "velocity_z")] = (data_velz[i].reshape((data_velz[i].shape[0], data_velz[i].shape[1], 1)), "cm/s")
      if quantities.__contains__("magnetic_field"):
        g[("gas", "magnetic_field_x")] = (data_Bx[i].reshape((data_Bx[i].shape[0], data_Bx[i].shape[1], 1)), "G")
        g[("gas", "magnetic_field_y")] = (data_By[i].reshape((data_By[i].shape[0], data_By[i].shape[1], 1)), "G")
        g[("gas", "magnetic_field_z")] = (data_Bz[i].reshape((data_Bz[i].shape[0], data_Bz[i].shape[1], 1)), "G")
      if quantities.__contains__("xray_emission"):
        g[("gas", "xray_0.1")]  = (data_emx001[i].reshape((data_emx001[i].shape[0], data_emx001[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_0.2")]  = (data_emx002[i].reshape((data_emx002[i].shape[0], data_emx002[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_0.3")]  = (data_emx003[i].reshape((data_emx003[i].shape[0], data_emx003[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_0.5")]  = (data_emx005[i].reshape((data_emx005[i].shape[0], data_emx005[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_1.0")]  = (data_emx010[i].reshape((data_emx010[i].shape[0], data_emx010[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_2.0")]  = (data_emx020[i].reshape((data_emx020[i].shape[0], data_emx020[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_2.5")]  = (data_emx025[i].reshape((data_emx025[i].shape[0], data_emx025[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_3.0")]  = (data_emx030[i].reshape((data_emx030[i].shape[0], data_emx030[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_5.0")]  = (data_emx050[i].reshape((data_emx050[i].shape[0], data_emx050[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        g[("gas", "xray_10.0")] = (data_emx100[i].reshape((data_emx100[i].shape[0], data_emx100[i].shape[1], 1)), "erg/cm**3/s/arcmin**2")
        print("adding emissivity to g,",i)

      i += 1

    data.close()

    ds = yt.load_amr_grids(grid_data, N_grids, length_unit=f"cm", geometry=("cylindrical"), axis_order=("r","z","theta"), sim_time=sim_time.value, bbox=bbox, time_unit="s", unit_system="cgs",periodicity=(False, False,True))
    return ds
##########################################################################


##########################################################################
# Create list of yt datasets representing a time series of simulation
# snapshots.
def get_ts(evolution, **kwargs):

    ds_list = []

    start = kwargs.get('start', 1)
    end = kwargs.get('end', 10)
    step = kwargs.get('step', 1)

    # Load the desired snapshots
    for i in range(start, end, step):
        ds_list.append(get_ds(evolution[i], quantities=kwargs.get("quantities", ["density"]), start_time=kwargs.get("start_time", 0)))

    print("Number of datasets: ", len(ds_list))
    # Create time series object
    return yt.DatasetSeries(ds_list)
    del ds_list
##########################################################################


