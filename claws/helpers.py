#!/usr/bin/env python
"""Helper functions.
"""

# Import modules
import itertools
import numpy as np

from claws.custom_exceptions import InputError

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Variables 
# ---------------------------------------------------------------------------- #

# Available concentration units (1. corresponds to program units)
available_conc_units = {'ng/L': 1., 'ug/L': 1.e3, 'mg/L': 1.e6, 'g/L': 1.e9}
                        
# Available mole units (1. corresponds to program units)
available_mol_units = {'umol': 1.e-6, 'mmol': 1.e-3, 'mol': 1., 'kmol': 1.e3}

# Available length, area, and volume units (1. corresponds to program units)
available_len_units = {'in': 0.0254, 'ft': 0.3048, 'm': 1., 'km': 1.e3,
                       'mi': 1.60934e3}
                       
available_area_units = {k + '^2': v**2 for (k,v) in available_len_units.items()}   
available_vol_units = {k + '^3': v**3 for (k,v) in available_len_units.items()}                                          
available_vol_units = {**available_vol_units, 'L': 1.e-3, 'M m^3': 1.e6}

# Available time units (1. corresponds to program units)
available_time_units = {'s': 1./3600., 'min': 1./60., 'h': 1., 'day': 24.,
                        'yr': 8760.}
                        
# Available mass units (1. corresponds to program units)
available_mass_units = {'ng': 1.e-12, 'ug': 1.e-9, 'mg': 1.e-6, 'g': 1.e-3,
                        'kg': 1., 'tonne': 1.e3}

# Available flow rate units                        
available_flowrate_units = {k1 + '/' + k2: v1/v2 for ((k1,v1), (k2,v2)) in
    itertools.product(available_vol_units.items(), available_time_units.items())} 


# ---------------------------------------------------------------------------- #
# Functions 
# ---------------------------------------------------------------------------- #

def indent(indent_level):
    """Indent a string according to an indentation level
    
    Arguments:
        indent_level: integer; indent level
    """
    return ''.join('\t' for _ in range(indent_level))
    
def flatten(inputlist):
    """Flatten a python list
    
    Arguments:
        inputlist: list; a python list of lists
    """
    return np.array([item for sublist in inputlist for item in sublist])
    
def convert_conc_to_prog_units(input_value, input_conc_units, class_name):
    """Convert an input concentration in prog. units
    
    Arguments:
        input_value: float; Concentration value to convert in prog. units
        
        input_conc_units: string; Input concentration units. Available
            concentration units are listed in available_conc_units
            
        class_name: string; name of the class calling this function
    """
    if input_conc_units not in available_conc_units.keys():
        raise InputError(input_conc_units,
            "Concentration units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_conc_units])))
                
    return input_value*available_conc_units.get(input_conc_units)
    
def convert_mol_to_prog_units(input_value, input_mol_units, class_name):
    """Convert an input molar quantity in prog. units
    
    Arguments:
        input_value: float; Molar quantity value to convert in prog. units
        
        input_mol_units: string; Input molar quantity units. Available
            molar quantity units are listed in available_mol_units
            
        class_name: string; name of the class calling this function
    """
    if input_mol_units not in available_mol_units.keys():
        raise InputError(input_mol_units,
            "Molar quantity units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_mol_units])))
                
    return input_value*available_mol_units.get(input_mol_units)
    
def convert_len_to_prog_units(input_value, input_len_units, class_name):
    """Convert an input length in prog. units
    
    Arguments:
        input_value: float; Length value to convert in prog. units
        
        input_len_units: string; Input length units. Available length units
            are listed in available_len_units
            
        class_name: string; name of the class calling this function
    """
    if input_len_units not in available_len_units.keys():
        raise InputError(input_len_units,
            "Length units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_len_units])))
    
    return input_value*available_len_units.get(input_len_units)
    
def convert_area_to_prog_units(input_value, input_area_units, class_name):
    """Convert an input area in prog. units
    
    Arguments:
        input_value: float; Area value to convert in prog. units
        
        input_area_units: string; Input area units. Available area units
            are listed in available_area_units
            
        class_name: string; name of the class calling this function
    """
    if input_area_units not in available_area_units.keys():
        raise InputError(input_len_units,
            "Area units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_area_units])))
    
    return input_value*available_area_units.get(input_area_units)
    
def convert_vol_to_prog_units(input_value, input_vol_units, class_name):
    """Convert an input volume in prog. units
    
    Arguments:
        input_value: float; Volume value to convert in prog. units
        
        input_vol_units: string; Input volume units. Available volume units
            are listed in available_vol_units
            
        class_name: string; name of the class calling this function
    """
    if input_vol_units not in available_vol_units.keys():
        raise InputError(input_len_units,
            "Volume units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_vol_units])))
    
    return input_value*available_vol_units.get(input_vol_units)
    
def convert_time_to_prog_units(input_value, input_time_units, class_name):
    """Convert an input time in prog. units
    
    Arguments:
        input_value: float; Time value to convert in prog. units
        
        input_time_units: string; Input time units. Available time units
            are listed in available_time_units
            
        class_name: string; name of the class calling this function
    """
    if input_time_units not in available_time_units.keys():
        raise InputError(input_time_units,
            "Time units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_time_units])))
    
    return input_value*available_time_units.get(input_time_units)
    
def convert_mass_to_prog_units(input_value, input_mass_units, class_name):
    """Convert an input mass in prog. units
    
    Arguments:
        input_value: float; Mass value to convert in prog. units
        
        input_mass_units: string; Input mass units. Available
            mass units are listed in available_mass_units
            
        class_name: string; name of the class calling this function
    """
    if input_mass_units not in available_mass_units.keys():
        raise InputError(input_mass_units,
            "Concentration units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_mass_units])))
    
    return input_value*available_mass_units.get(input_mass_units)
    
def convert_conc(input_value, output_conc_units, class_name):
    """Convert an input concentration from prog units to output_conc_units
    
    Arguments:
        input_value: float; Concentration value to convert in output_conc_units
        
        output_conc_units: string; Output concentration units. Available
            concentration units are listed in available_conc_units
            
        class_name: string; name of the class calling this function
    """
    if output_conc_units not in available_conc_units.keys():
        raise InputError(output_conc_units,
            "Concentration units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_conc_units])))
    
    return input_value/available_conc_units.get(output_conc_units)
    
def convert_mol(input_value, output_mol_units, class_name):
    """Convert an input molar quantity from prog units to output_mol_units
    
    Arguments:
        input_value: float; Molar quantity value to convert in output_mol_units
        
        output_mol_units: string; Output molar quantity units. Available
            molar quantity units are listed in available_mol_units
            
        class_name: string; name of the class calling this function
    """
    if output_mol_units not in available_mol_units.keys():
        raise InputError(output_mol_units,
            "Molar quantity units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_mol_units])))
    
    return input_value/available_mol_units.get(output_mol_units)
    
def convert_len(input_value, output_len_units, class_name):
    """Convert an input length from prog units to output_len_units
    
    Arguments:
        input_value: float; Length value to convert in output_len_units
        
        output_len_units: string; Input length units. Available length units
            are listed in available_len_units
            
        class_name: string; name of the class calling this function
    """
    if output_len_units not in available_len_units.keys():
        raise InputError(output_len_units,
            "Length units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_len_units])))
    
    return input_value/available_len_units.get(output_len_units)
    
def convert_area(input_value, output_area_units, class_name):
    """Convert an input area from prog units to output_area_units
    
    Arguments:
        input_value: float; Area value to convert in output_area_units
        
        output_area_units: string; Input area units. Available area units
            are listed in available_area_units
            
        class_name: string; name of the class calling this function
    """
    if output_area_units not in available_area_units.keys():
        raise InputError(output_area_units,
            "Area units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_area_units])))
    
    return input_value/available_area_units.get(output_area_units)
    
def convert_vol(input_value, output_vol_units, class_name):
    """Convert an input volume from prog units to output_vol_units
    
    Arguments:
        input_value: float; Volume value to convert in output_vol_units
        
        output_vol_units: string; Input volume units. Available volume units
            are listed in available_vol_units
            
        class_name: string; name of the class calling this function
    """
    if output_vol_units not in available_vol_units.keys():
        raise InputError(output_vol_units,
            "Volume units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_vol_units])))
    
    return input_value/available_vol_units.get(output_vol_units)
    
def convert_time(input_value, output_time_units, class_name):
    """Convert an input time from prog units to output_time_units
    
    Arguments:
        input_value: float; Time value to convert in output_time_units
        
        output_time_units: string; Output time units. Available time units
            are listed in available_time_units
            
        class_name: string; name of the class calling this function
    """
    if output_time_units not in available_time_units.keys():
        raise InputError(output_time_units,
            "Time units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_time_units])))
    
    return input_value/available_time_units.get(output_time_units)
    
def convert_mass(input_value, output_mass_units, class_name):
    """Convert an input mass from prog units to output_mass_units
    
    Arguments:
        input_value: float; Mass value to convert in output_mass_units
        
        output_mass_units: string; Output mass units. Available mass units
            are listed in available_mass_units
            
        class_name: string; name of the class calling this function
    """
    if output_mass_units not in available_mass_units.keys():
        raise InputError(output_mass_units,
            "Mass units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_mass_units])))
    
    return input_value/available_mass_units.get(output_mass_units)
    
def convert_flowrate(input_value, output_fr_units, class_name):
    """Convert an input flow rate from prog units to output_fr_units
    
    Arguments:
        input_value: float; Flow rate value to convert in output_fr_units
        
        output_fr_units: string; Input flow rate units. Available flow rate
            units are listed in available_flowrate_units
            
        class_name: string; name of the class calling this function
    """
    if output_fr_units not in available_flowrate_units.keys():
        raise InputError(output_fr_units,
            "Flow rate units not recognised in {}. "\
            "Options are: {}".format(class_name, ', '.join([units for units in
                available_flowrate_units])))
    
    return input_value/available_flowrate_units.get(output_fr_units)
