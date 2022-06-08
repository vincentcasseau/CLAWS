#!/usr/bin/env python
"""Abstract class defining a sea lice treatment.
"""

# Import modules
import re
import numpy as np
import claws.helpers as helpers
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
# Classes 
# ---------------------------------------------------------------------------- #

class Treatment(object):
    def __init__(self, tarpaulin_height, tarpaulin_radius, seeding_times,
                 nparticles, Chemicals, input_len_units='m',
                 input_time_units='h', reference=""):
        # __tarpaulin_height: float; Tarpaulin height (meters)
        # A random number in the range [-__tarpaulin_height, 0] will the drawn
        # to decide on the initial depth position of a particle
        self.__tarpaulin_height = tarpaulin_height
        # __tarpaulin_radius: float; Tarpaulin radius (meters)
        # Particles will be released within this radius around the prescribed
        # site location
        self.__tarpaulin_radius = tarpaulin_radius
        # __seeding_times: int/float/list/np.ndarray; Seeding times (hours)
        self.__seeding_times = np.array(seeding_times, dtype=float)
        # __nparticles: int/float; Number of particles to seed at each seeding
        # time
        self.__nparticles = nparticles
        # __chemicals_obj: ChemicalSubstance; Chemicals used at each seeding
        # time
        self.__chemicals_obj = Chemicals
        # __reference: string; Reference. default: empty string
        self.__reference = reference
        
        self._sanitize(input_len_units, input_time_units)
        
    def __str__(self, indent_lvl=1):
        return """{1}\n{0}Tarpaulin height (m) = {2}\
            \n{0}Tarpaulin radius (m) = {3}\n{0}Seeding times (h) = {4}\
            \n{0}Number of particles = {5}\n{0}Chemicals = {6}\
            \n{0}Marker index = {7}\n{0}Reference = {8}""".format(
            helpers.indent(indent_lvl), self.name(),self.__tarpaulin_height,
            self.__tarpaulin_radius, self.__seeding_times, self.__nparticles,
            [chem.name() for chem in self.__chemicals_obj], self.__marker_index,
            self.__reference)
    
    def tarpaulin_height(self):
        return self.__tarpaulin_height
        
    def tarpaulin_radius(self):
        return self.__tarpaulin_radius
        
    def nparticles(self):
        return self.__nparticles
        
    def seeding_times(self):
        return self.__seeding_times
        
    def Chemicals(self):
        return self.__chemicals_obj
        
    def marker_index(self):
        return self.__marker_index
        
    def set_marker_index(self, value, indices=None):
        if indices is None:
            self.__marker_index = [value for mi in self.__marker_index]
        else:
            self.__marker_index[indices] = value
      
    def _sanitize(self, input_len_units, input_time_units):
        ntreatments = len(self.__seeding_times)
    
        if type(self.__tarpaulin_height) in [float, int]:
            if self.__tarpaulin_height <= 0.0:
                raise InputError(self.__tarpaulin_height, 
                    "Tarpaulin height must be a positive number")
            self.__tarpaulin_height = np.full(ntreatments,
                                              self.__tarpaulin_height)
        elif (type(self.__tarpaulin_height) in [list, np.ndarray] and
              len(self.__tarpaulin_height) != ntreatments):
            raise IndexError("Incorrect size for tarpaulin_height list:", 
                             self.__tarpaulin_height, " in class",
                             type(self).__name__)
                             
        if type(self.__tarpaulin_radius) in [float, int]:
            if self.__tarpaulin_radius <= 0.0:
                raise InputError(self.__tarpaulin_radius, 
                    "Tarpaulin radius must be a positive number")
            self.__tarpaulin_radius = np.full(ntreatments,
                                              self.__tarpaulin_radius)
        elif (type(self.__tarpaulin_radius) in [list, np.ndarray] and
              len(self.__tarpaulin_radius) != ntreatments):
            raise IndexError("Incorrect size for tarpaulin_radius list:", 
                             self.__tarpaulin_radius, " in class",
                             type(self).__name__)
        
        if type(self.__nparticles) in [float, int]:
            if self.__nparticles <= 0.0:
                raise InputError(self.__nparticles, 
                    "Number of particles to seed must be a positive number")
            self.__nparticles = np.full(ntreatments, self.__nparticles)
        elif (type(self.__nparticles) in [list, np.ndarray] and
              len(self.__nparticles) != ntreatments):
            raise IndexError("Incorrect size for nparticles list:", 
                             self.__nparticles, " in class",
                             type(self).__name__)
                             
        self.__nparticles = np.array(self.__nparticles)
        nparticles_rounded = np.array(self.__nparticles, dtype=int)
        nparticles_diff = self.__nparticles - nparticles_rounded
        for i in range(len(self.__nparticles)):
            if np.random.rand() < nparticles_diff[i]:
                nparticles_rounded[i] += 1
        self.__nparticles = nparticles_rounded

        if (type(self.__chemicals_obj) in [list, np.ndarray] and
                len(self.__chemicals_obj) != ntreatments):
            raise IndexError("Incorrect size for chemicals list:", 
                             [c.name() for c in self.__chemicals_obj],
                             " in class", type(self).__name__)
        else:
            self.__chemicals_obj = np.full(ntreatments, self.__chemicals_obj)

        # Conversions to program units
        self.__tarpaulin_height = helpers.convert_len_to_prog_units(
            self.__tarpaulin_height, input_len_units, self.name())
        self.__tarpaulin_radius = helpers.convert_len_to_prog_units(
            self.__tarpaulin_radius, input_len_units, self.name())
        self.__seeding_times = helpers.convert_time_to_prog_units(
            self.__seeding_times, input_time_units, self.name())
        
        self.__marker_index = np.full(ntreatments, -1)
                
    def name(self):
        return ' '.join(re.findall('([A-Z][a-z]+)', type(self).__name__))
        
    def tarpaulin_area(self):
        return np.pi*self.__tarpaulin_radius**2
    
    def tarpaulin_volume(self):
        return self.__tarpaulin_height*self.tarpaulin_area()
        

# ---------------------------------------------------------------------------- #
# Functions 
# ---------------------------------------------------------------------------- #        
        
def compute_marker_indices(seeding_locations, chemicals):
    """Compute marker indices for all seeding locations
    
    Arguments:
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
        
        chemicals: list of ChemicalSubstance objects as implemented in
            claws.ChemicalSubstance; All chemicals used in 
            any farms
    """
    nchemicals = len(chemicals)
    if nchemicals > 1:
        for loc in seeding_locations:
            for c in range(nchemicals):
                indices = np.where(loc.Treatment().Chemicals() == chemicals[c])[0]
                loc.Treatment().set_marker_index(c, indices)
    else:
        for loc in seeding_locations:
            loc.Treatment().set_marker_index(0)
        
    return seeding_locations
