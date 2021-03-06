#!/usr/bin/env python
"""Main module defining dictionaries to record user inputs and functions to
sanitise these inputs
"""

# Import modules
import os
import numpy as np
from pyproj import Proj
import opendrift

from claws.custom_exceptions import InputError
from claws.helpers import *

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

output_options = {
    # Period at which to create sub-histograms and plot concentration maps.
    # default is None, in which case the time_bins entry is considered
    # If both are None, only the last time step is considered
    "output_period": None,
    # time_bins: integer or list; output time bins, as defined by the user.
    # default is None. A value of None means that the variable 'output_period'
    # is considered instead 
    "time_bins": None,
    # Concentration units: string; ng/L, ug/L, mg/L or g/L. default is ng/L
    "concentration_units": "ng/L",
    # Length units: string; m, km, mi, ft. default is km
    "length_units": "km",
    # Time units: string; s, min, h or day. default is h
    "time_units": "h",
    # write_time_series_to_file: bool; whether or not to print time series
    # (peak concentration and area greater than EQS/MAC) to a .dat file
    "write_time_series_to_file": True,
}

animation = {
    # animate: bool; whether or not to animate concentration fields and vertical
    # distribution histograms. default is True 
    "animate": True,
    # frame_rate: integer; number of frames per second. default is 5
    "frame_rate": 5,
    # image_rate: integer; number of images per second, if different from
    # 'frame_rate', images get duplicated. default is 5
    # Example: an 'image_rate' of 1 means that the same image is displayed
    # during a full second and will be duplicated X times with a frame rate of X
    "image_rate": 5,
}

# Concentration, length, and time conversion factors between the program units
# and the user-defined output units set in the output_options dictionary
_unit_factors = [1., 1., 1.]


# ---------------------------------------------------------------------------- #
# Functions 
# ---------------------------------------------------------------------------- #

def sanitise_working_folder(working_folder):
    """Sanitise working folder and create media directory
    
    Arguments:
        working_folder: string; working folder
    """
    if working_folder[-1] != '/':
        working_folder += '/'
        
    media_folder = create_media_directory(working_folder)
    return working_folder, media_folder
    
def create_media_directory(working_folder):
    """Create a folder to store media: images, videos created during
    post-processing
    """
    media_folder = os.path.join(working_folder, 'media/')
    try:
        os.mkdir(media_folder)
    except FileExistsError:
        pass
    return media_folder
    
def sanitise_output_file(working_folder, of, script_args):
    """Sanitise output file
    
    Arguments:
        working_folder: string; working folder
        
        of: string; output file
        
        script_args: list of strings or integer value; run script command-line
            arguments or number of independent simulations
    """
    if type(script_args) is list:
        if len(script_args) > 1:
            simulation_index = '_' + str(script_args[1])
        else:
            simulation_index = ''
        outfile = working_folder + of + simulation_index + '.nc'
    elif type(script_args) is int:
        nsimulations = script_args
        if nsimulations > 1:
            outfile = working_folder + of + '_1.nc'
        else:
            outfile = working_folder + of + '.nc'
    else:
        raise InputError(script_args,
            'Wrong argument type passed to claws.claws.sanitise_output_file')
    return outfile
    
def sanitise_output_options(ndt):
    """Sanitise the output_options dictionary
    
    Arguments:
        ndt: int; number of time steps stored in the .nc file
    """
    
    _sanitise_time_bins(ndt)
    _sanitise_output_units()
    
def _sanitise_time_bins(ndt):
    """Sanitise time bins
    
    Arguments:
        ndt: int; number of time steps stored in the .nc file
    """    
    global output_options
    
    # Shorthand
    time_bins = output_options["time_bins"]
    
    # If no time bins are specified by the user, the global variable
    # 'output_period' is used
    if time_bins is None:
        if output_options["output_period"] is None:
            time_bins = np.arange(ndt-1, ndt)
        else:    
            time_bins = np.arange(0, ndt, output_options["output_period"])
    elif type(time_bins) == int and time_bins == -1:
        time_bins = np.arange(ndt-1, ndt)
    else:
        if len(np.shape(time_bins)) > 1:
            time_bins = flatten(time_bins)
        else:
            time_bins = np.array(time_bins)
        
    output_options["time_bins"] = time_bins
    
def _sanitise_output_units():
    """Sanitise the units defined by the user in the output_options dictionary
    """
    global _unit_factors
    global output_options
    
    # Output concentration unit
    output_conc_units = output_options["concentration_units"]
    conc_factor = convert_conc(1.0, output_conc_units, "_sanitise_output_units")
    
    # Output length unit
    output_len_units = output_options["length_units"]
    len_factor = convert_len(1.0, output_len_units, "_sanitise_output_units")
    
    # Output time unit
    output_time_units = output_options["time_units"]
    available_time_units = {'s': 1., 'min': 1./60., 'h': 1./3600.,
                            'day': 1./86400.}
        
    if output_time_units not in available_time_units.keys():
        raise InputError(output_time_units,
            "Time units not recognised in output_options "\
            "dictionary. Options are: {}".format(', '.join([unit for unit in
            available_time_units])))
            
    time_factor = available_time_units.get(output_time_units, None)
        
    # Program to output conversion factors
    _unit_factors = [conc_factor, len_factor, time_factor]
    
def get_unit_factors():
    """Return the unit conversion factors between the program units and the
    user-defined output units
    """
    return _unit_factors
