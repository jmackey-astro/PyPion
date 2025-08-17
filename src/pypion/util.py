from .SiloHeader_data import OpenData
from .ReadData import ReadData
import astropy.units as unit
import os
import glob
import re

######################################################################################
# batch silo files
######################################################################################
def batch_silos(dir, filebase, start_time=None, finish_time=None, time_unit=None,
                out_frequency=None):
    '''
    Organizes silo files into groups based on their corresponding time instants.
    This function scans a specified directory for silo files with a given base name,
    and if the simulation involves multiple grid levels, it groups the files by time instant.
    For each time instant, it creates a sublist containing silo files from different grid levels.

    The function allows running from a particular start time to a finish time, enabling users to
    process specific intervals of simulation data. Additionally, one can restart the processing
    by adjusting the start time, making it flexible for iterative analysis or simulations.

    :param dir: Directory path where the silo files are located.
    :param filebase: Base name of the silo files to search for.
    :param start_time: The start time for processing, in the given time_unit.
    :param finish_time: The finish time for processing, in the given time_unit.
    :param time_unit: The time unit for start_time and finish_time ('sec', 'yr', or 'kyr').
    :return: A list of lists, where each sublist contains silo files for a specific time instant,
             with each file in the sublist representing a different grid level.
    :raises: RuntimeError if no silo files are found or if expected files for certain levels are missing.
    '''
    conversion_factors = {
        'kyr': 3.154e+10,
        'yr': 3.154e+7,
        'sec': 1
    }

    inverse_conversion_factors = {
        'kyr': 1 / 3.154e+10,
        'yr': 1 / 3.154e+7,
        'sec': 1
    }

    # Convert start_time and finish_time to seconds based on time_unit
    factor = conversion_factors.get(time_unit, 1)  # Default to 1 if time_unit is None or unrecognized
    start_time_sec = start_time * factor if start_time is not None else None
    finish_time_sec = finish_time * factor if finish_time is not None else None

    print(" ---------------------------")
    print(" batching silo files into time instances:")
    print(f" starting time: {start_time} {time_unit}")
    print(f" finishing time: {finish_time} {time_unit}")
    print(f" output frequency: {out_frequency}")

    # Construct search pattern for silo files
    search_pattern = os.path.join(dir, f"{filebase}_*.silo")
    all_silos = glob.glob(search_pattern)

    if not all_silos:
        nebula_exit_with_error(f"no '{filebase}' silo files found in '{dir}'")

    # Assume the number of grid levels; this may need to be adjusted based on actual usage
    # Open header data
    header_data = OpenData(all_silos)
    # Set the directory to '/header'
    header_data.db.SetDir('/header')
    # Retrieve what coordinate system is used
    coord_sys = header_data.db.GetVar("coord_sys")
    # Retrieve no of nested grid levels
    Nlevels = header_data.db.GetVar("grid_nlevels")
    # close the object
    header_data.close()
    # set number of time instance
    Ninstances = 0

    # for uniform grid *********************************************************************
    if Nlevels == 1:
        print(f" grid: uniform")
        silos = sorted(all_silos)
        batched_silos = [[silo] for silo in silos]
        selected_silos = []
        for i, silo in enumerate(batched_silos):
            data = ReadData(silo)
            basic = data.get_1Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()

            if (start_time_sec is not None and sim_time < start_time_sec) or \
                    (finish_time_sec is not None and sim_time > finish_time_sec):
                continue  # Skip this silo if it doesn't meet the criteria
            selected_silos.append(silo)

        batched_silos = selected_silos
        # Keep files based on the output frequency if specified
        if out_frequency is not None:
            Ninstances = len(batched_silos)
            # indices to keep: multiples of out_freq, plus the first and last index
            #indices_to_keep = sorted(set(range(0, Ninstances, out_frequency)) | {0, Ninstances - 1})
            indices_to_keep = sorted(set(range(0, Ninstances, out_frequency)))
            batched_silos = [batched_silos[i] for i in indices_to_keep]

        Ninstances = len(batched_silos)
        print(f" {Ninstances} time instances between {start_time} {time_unit} and {finish_time} {time_unit}")
        print(" batching completed")
        return batched_silos

    # for nested grid *********************************************************************
    else:
        print(f" grid: nested with {Nlevels} levels")
        # Pattern to match level 00 files
        pattern = re.compile(r'_level00_0000\.\d+\.silo$')
        # Find and sort level 00 files, one per time instant
        level00_instants = [file for file in all_silos if pattern.search(file)]
        batched_silos = [[file] for file in sorted(level00_instants)]

        # check if level 00 batched silos are within the asked time range
        selected_silos = []
        for i, silo in enumerate(batched_silos):
            data = ReadData(silo)
            if coord_sys == 3:
                basic = data.get_1Darray('Density')
            elif coord_sys == 2:
                basic = data.get_2Darray('Density')
            elif coord_sys == 1:
                basic = data.get_3Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()

            if (start_time_sec is not None and sim_time < start_time_sec) or \
                    (finish_time_sec is not None and sim_time > finish_time_sec):
                continue  # Skip this silo if it doesn't meet the criteria
            selected_silos.append(silo)

        batched_silos = selected_silos

        # Keep files based on the output frequency if specified
        if out_frequency is not None:
            Ninstances = len(batched_silos)
            # indices to keep: multiples of out_freq, plus the first and last index
            indices_to_keep = sorted(set(range(0, Ninstances, out_frequency)) | {0, Ninstances - 1})
            batched_silos = [batched_silos[i] for i in indices_to_keep]

        # appending other level silos to corresponding time instant
        for i, instant in enumerate(batched_silos):
            # Extract the time instant from the level 00 file
            file_name = instant[0].split('/')[-1]
            instant_extension = file_name.split('_')[-1].replace('.silo', '')

            # Find and append silo files for the same time instant across other levels
            for level in range(1, Nlevels):
                level_pattern = re.compile(f'level{str(level).zfill(2)}_{instant_extension}.silo')

                # Find the file for the current level
                level_instant_file = next((file for file in all_silos if level_pattern.search(file)), None)

                if level_instant_file is None:
                    (f"missing silo file for level {level} instant {instant_extension}")
                else:
                    # Append the found file to the corresponding time instant group
                    batched_silos[i].append(level_instant_file)

        # common to both uniform and nested grid
        if start_time is None:
            data = ReadData(batched_silos[0])
            if coord_sys == 3:
                basic = data.get_1Darray('Density')
            elif coord_sys == 2:
                basic = data.get_2Darray('Density')
            elif coord_sys == 1:
                basic = data.get_3Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()
            start_time = sim_time * inverse_conversion_factors.get(time_unit, 1)

        if finish_time is None:
            data = ReadData(batched_silos[-1])
            if coord_sys == 3:
                basic = data.get_1Darray('Density')
            elif coord_sys == 2:
                basic = data.get_2Darray('Density')
            elif coord_sys == 1:
                basic = data.get_3Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()
            finish_time = sim_time * inverse_conversion_factors.get(time_unit, 1)


        if finish_time is not None:
            data = ReadData(batched_silos[-1])
            if coord_sys == 3:
                basic = data.get_1Darray('Density')
            elif coord_sys == 2:
                basic = data.get_2Darray('Density')
            elif coord_sys == 1:
                basic = data.get_3Darray('Density')
            sim_walltime_sec = (basic['sim_time'] * unit.s).value
            data.close()
            #finish_time_sec = finish_time * conversion_factors.get(time_unit, 1)
            if finish_time_sec > sim_walltime_sec:
                sim_walltime = sim_walltime_sec * inverse_conversion_factors.get(time_unit, 1)
                nebula_warning(
                    f"specified finish time {float(finish_time):.3f} {time_unit} exceeds the simulation"
                    f" walltime {float(sim_walltime):.3f} {time_unit}"
                )

        Ninstances = len(batched_silos)
        print(f" {Ninstances} time instances between {start_time:.3f} {time_unit} and {finish_time:.3f} {time_unit}")
        print(" silo batching completed")
        if not batched_silos:
            nebula_exit_with_error('no silo files found in the specified time range, check your selection criteria')
        return batched_silos