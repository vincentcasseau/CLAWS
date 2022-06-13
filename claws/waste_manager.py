#!/usr/bin/env python
"""Root class to manage waste.
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

class WasteManager(object):
    def __init__(self, feed_requirement, feed_water_percentage=9.0,
                 feed_waste_percentage=3.0, feed_absorbed_percentage=85.0,
                 feed_carbon_percentage=49.0, faeces_carbon_percentage=30.0,
                 reference="", name=""):
        # self.__name: string; Treatment name. default is class name
        self.__name = name
        # __feed_requirement: float; Feed requirement
        self.__feed_requirement = feed_requirement
        # __feed_water_content: float; Feed water percentage. default is 0.09
        self.__feed_water_content = feed_water_percentage/100.
        # __feed_wastage_rate: float; Feed wastage rate. default is 0.03
        self.__feed_wastage_rate = feed_waste_percentage/100.
        # __feed_absorbed_rate: float; Feed absorbed percentage. default is 0.85
        self.__feed_absorbed_rate = feed_absorbed_percentage/100.
        # __feed_carbon_content: float; Feed carbon content. default is 0.49
        self.__feed_carbon_content = feed_carbon_percentage/100.
        # __faeces_carbon_rate: float; Faeces carbon rate. default is 0.3
        self.__faeces_carbon_rate = faeces_carbon_percentage/100.
        # __reference: string; Reference. default: empty string
        self.__reference = reference
        
        self._sanitize()
        
    def __str__(self, indent_lvl=1):
        return """{1}\
            \n{0}Reference = {2}""".format(
            helpers.indent(indent_lvl), self.name(), self.__reference)
    
    def _sanitize(self):
        if not self.__name:
            self.__name = ' '.join(re.findall('([A-Z][a-z]+)',
                                   type(self).__name__))
                                   
        if type(self.__feed_requirement) is float:
            if self.__feed_requirement < 0.0:
                raise InputError(self.__feed_requirement, 
                    "Feed requirement must be a positive number")

    def name(self):
        return self.__name
        
    def convert_mass_io(mass, input_mass_units, output_mass_units):
        mass = convert_mass_to_prog_units(mass, input_mass_units, self.name())
        return convert_mass(mass, output_mass_units, self.name())
        
    def waste_solids_ratio(self):
        return (1.0 - self.__feed_water_content)*self.__feed_wastage_rate
        
    def discharged_carbon_ratio(self):
        return self.waste_solids_ratio()*self.__feed_carbon_content
            
    def excreted_solids_ratio(self):
        return (1.0 - self.__feed_water_content)*\
            (1.0 - self.__feed_wastage_rate)*(1.0 - self.__feed_absorbed_rate)
    
    def excreted_carbon_ratio(self):
        return self.excreted_solids_ratio()*self.__faeces_carbon_rate
    
    def waste_solids(self, feed_load, input_mass_units='kg',
                     output_mass_units='kg'):
        return convert_mass_io(feed_load*self.waste_solids_ratio(),
                               input_mass_units, output_mass_units)
        
    def discharged_carbon(self, feed_load, input_mass_units='kg',
                          output_mass_units='kg'):
        return convert_mass_io(feed_load*self.discharged_carbon_ratio(),
                               input_mass_units, output_mass_units)
                               
    def excreted_solids(self, feed_load, input_mass_units='kg',
                        output_mass_units='kg'):
        return convert_mass_io(feed_load*self.excreted_solids_ratio(),
                               input_mass_units, output_mass_units)
                               
    def excreted_carbon(self, feed_load, input_mass_units='kg',
                        output_mass_units='kg'):
        return convert_mass_io(feed_load*self.excreted_carbon_ratio(),
                               input_mass_units, output_mass_units)
