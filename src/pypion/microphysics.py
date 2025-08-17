from . import util
from .SiloHeader_data import OpenData
import re
import numpy as np

# PION constants
pion_elements = {'H': "Hydrogen", 'He': "Helium", 'C': "Carbon", 'N': "Nitrogen",
                   'O': "Oxygen", 'Ne': "Neon", 'Si': "Silicon", 'S': "Sulfur",
                   'Fe': "Iron"}

mass = {'H': 1.6738e-24, 'He': 6.6464768e-24, 'C': 1.994374e-23, 'N': 2.325892e-23,
        'O': 2.6567628e-23, 'Ne': 3.3509177e-23, 'Si': 4.6637066e-23, 'S': 5.3245181e-23,
        'Fe': 9.2732796e-23}

atomic_number = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10, 'Si': 14, 'S': 16,
                 'Fe': 26}

top_level_ions = {'H1+', 'He2+', 'C6+', 'N7+', 'O8+', 'Ne10+', 'Si14+', 'S16+', 'Fe26+'}

class microphysics():

    ######################################################################################
    # initialize microphysics
    ######################################################################################
    def __init__(self, silo_set):
        self.silo_set = silo_set
        self.chemistry_container = {}

    ######################################################################################
    # load chemistry
    ######################################################################################
    def load_chemistry(self):
        '''
        This method extracts information related to the chemistry and chemical tracers,
        transforming the chemical tracer names to a format that PyPion can directly
        read from the Silo file. This method can be included in the next version of
        PyPion.

        Parameters
        ----------
        instant_silo_set : The instance for which chemical data is to be extracted

        Returns
        -------
        Generates and stores the following in self.chemistry_container:
        - 'dynamics': Dynamics data retrieved from the Silo file
        - 'chemistry': Chemistry flag indicating if chemistry data is available
        - 'E_update': Energy update information (if chemistry flag is true)
        - 'chemistry_code': The code indicating the type of chemistry (if chemistry flag is true)
        - 'microphysics': List of microphysics processes (if chemistry flag is true)
        - 'Ntracers': Number of chemical tracers
        - 'mpv10_elements': List of elements identified for MPv10 chemistry code
        - 'mpv10_tracers': List of tracers corresponding to each element for MPv10 chemistry code
        '''

        # Open the data for the first silo instant silo
        header_data = OpenData(self.silo_set[0])
        # Set the directory to '/header'
        header_data.db.SetDir('/header')
        # print(header_data.header_info())
        # Retrieve the value of "EP_chemistry" from the header data
        chemistry_flag = header_data.db.GetVar("EP_chemistry")
        self.chemistry_container['chemistry'] = chemistry_flag

        # Define the list of process variable names
        processes = ['EP_coll_ionisation', 'EP_rad_recombination',
                     'EP_cooling', 'EP_raytracing', 'EP_phot_ionisation',
                     'EP_charge_exchange']

        # Define the list of process names corresponding to the process variable names
        processes_name = ['coll_ionisation', 'rad_recombination', 'cooling',
                          'raytracing', 'phot_ionisation', 'charge_exchange']

        # Check if chemistry_flag is true
        if chemistry_flag:
            # Retrieve the value of "EP_update_erg"
            energy_update = header_data.db.GetVar("EP_update_erg")
            # save the energy_update value in the chemistry_container dictionary
            self.chemistry_container['E_update'] = energy_update
            # Retrieve the value of "chem_code"
            chemistry_code = header_data.db.GetVar("chem_code")[0]
            # save the chemistry_code value in the chemistry_container dictionary
            self.chemistry_container['chemistry_code'] = chemistry_code

            # Initialize an empty list to store microphysics processes
            microphysics = []
            # Check if the chemistry_code is not 'MPv10'
            if not chemistry_code == 'MPv10':
                # Exit with an error if the chemistry_code is not 'MPv10'
                util.pypion_exit_with_error(" PION is not running NEMO v1.0; NelubaPy functionality is limited.")
            else:
                # print the chemistry code
                print(f" Microphysics module loaded")
                print(" Module : NEMO")

                # Loop through each process
                for index, process in enumerate(processes):
                    # Check if the process variable exists in the header data
                    if header_data.db.GetVar(process):
                        # Append the corresponding process name to the microphysics list
                        microphysics.append(processes_name[index])

                # save the microphysics list in the chemistry_container dictionary
                self.chemistry_container['microphysics'] = microphysics
                # Retrieve the number of tracers
                Ntracers = header_data.db.GetVar('num_tracer')
                # elements in the tracer list
                tracer_elements = []
                # mass_fraction
                mass_fractions = {}
                # list of element wise tracer list
                elementWiseTracers = [[] for _ in range(len(pion_elements))]
                # list of element names from the nebula_elements dictionary keys
                element_list = list(pion_elements.keys())
                # save the number of tracers in the chemistry_container dictionary
                self.chemistry_container['Ntracers'] = Ntracers
                # If verbose is enabled, print the number of chemical tracers


                # Loop through each tracer index
                for i in range(Ntracers):
                    # create a tracer index string with leading zeros
                    tracer_index = f'Tracer{i:03}'
                    # retrieve the tracer value
                    chem_tracer = header_data.db.GetVar(tracer_index)[0]

                    # check if the tracer is an element ('X' denoting elemental mass fraction)
                    if 'X' in chem_tracer and chem_tracer.replace("_", "").replace("X", "") in pion_elements:
                        # extract the element name
                        element = chem_tracer.replace("_", "").replace("X", "")
                        tracer_elements.append(element)
                        # get the full element name from the nebula_elements dictionary
                        mass_fractions[element] = f'Tr{i:03}_' + chem_tracer
                        # get the index of the element in the element_list
                        element_index = element_list.index(element)
                        # append the tracer with the corresponding element to the mpv10tracers list
                        if 0 <= element_index < len(elementWiseTracers):
                            elementWiseTracers[element_index].append(f'Tr{i:03}_' + chem_tracer)

                    # check if the tracer is a corresponding ion
                    if re.sub(r'\d{1,2}\+', '', chem_tracer) in pion_elements:
                        self.chemistry_container[chem_tracer] = f'Tr{i:03}_' + chem_tracer.replace('+', 'p')
                        # extract the element name
                        element = re.sub(r'\d{1,2}\+', '', chem_tracer)
                        # get the index of the element in the element_list
                        element_index = element_list.index(element)
                        # gppend the tracer with the corresponding ion to the mpv10tracers list
                        elementWiseTracers[element_index].append(f'Tr{i:03}_' + chem_tracer.replace('+', 'p'))

                print(f" Tracers ({Ntracers}): {', '.join(tracer_elements)}")
                # save mass fraction to chemistry_container dictionary
                # self.chemistry_container['mass_fractions'] = mass_fractions
                self.element_list = tracer_elements
                self.chemistry_container['elemental_mass_fraction'] = mass_fractions
                self.chemistry_container['tracer_elements'] = tracer_elements
                self.chemistry_container['element_wise_tracer_list'] = elementWiseTracers
        else:
            util.pypion_exit_with_error('Chemistry flag is inactive in the loaded simulation file.')
        header_data.close()

    ######################################################################################
    # get elements
    ######################################################################################
    def get_elements(self):
        return np.array(self.chemistry_container['tracer_elements'])

    ######################################################################################
    # get chemical tracers
    ######################################################################################
    def get_chemical_tracers(self):
        """
        Retrieve the list of chemical tracer strings for each tracer in the chemistry
        container dictionary, processed element by element. Each sublist starts with the
        mass fraction of the element followed by the tracers.

        Returns:
            list of lists: Each sublist contains the mass fraction followed by the values of
            the tracers for a specific element.
        """
        elements = self.get_elements()
        tracers = []

        for element in elements:
            # Retrieve tracers for the element
            element_tracers = [self.chemistry_container[f"{element}{q}+" if q > 0 else element]
                               for q in range(atomic_number[element])]

            tracers.append(element_tracers)

        return tracers
