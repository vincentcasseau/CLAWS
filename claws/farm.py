#!/usr/bin/env python
"""Abstract class defining a farm.
"""

# Import modules
import re
import numpy as np
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

# Total nutrient discharge (kg N per tonne of fish produced)
total_nutrient_discharge = {"salmon": 48.2, "halibut": 67.1, "turbot": 86.9,
                            "cod": 72.3, "haddock": 72.3}

# Total nutrient discharge expressed as a 'species factor' relative to salmon
species_factor = {"salmon": 1.0, "halibut": 1.4, "turbot": 1.8, "cod": 1.5,
                  "haddock": 1.5}


# ---------------------------------------------------------------------------- #
# Functions 
# ---------------------------------------------------------------------------- #

def nutrient_enhancement_index(ECE):
    """Return the nutrient enhancement index using an input ECE in umol/L
    
    References:
        https://consultation.sepa.org.uk/permits/loch-long-salmon-beinn-
              reithe-car-application/user_uploads/210113-bnrt-par-nutrient-
              modelling_redacted.pdf, p. 13
    """
    if ECE > 10.:
        return 5
    elif ECE > 3.:
        return 4
    elif ECE > 1.:
        return 3
    elif ECE > 0.3:
        return 2
    elif ECE > 0:
        return 1
    return 0


# ---------------------------------------------------------------------------- #
# Classes 
# ---------------------------------------------------------------------------- #

class Farm(object):
    def __init__(self, Site, Treatment=None, yearly_fish_production=0.0,
                 total_nutrient_discharge=0.0, input_mass_units='kg',
                 reference="", name=""):
        # __name: string; Farm name. default is class name
        self.__name = ' '.join(re.findall('([A-Z][a-z]+)', type(self).__name__))
        # __site_obj: Site object; Farm location
        self.__site_obj = Site
        # __treatment_obj: Treatment object; Treatment being applied at the
        # farm. default: None
        self.__treatment_obj = Treatment
        # self.__yearly_fish_production: float; Yearly fish production (kg)
        self.__yearly_fish_production = yearly_fish_production
        # Combined source of nitrogen from dissolved ammonia and particulate
        # waste emissions (units: kg of nitrogen per tonne of fish produced)
        self.__total_nutrient_discharge = total_nutrient_discharge
        # __reference: string; Reference. default: empty string
        self.__reference = reference
        
        self._sanitize(input_mass_units)
        
    def __str__(self):
        return """{}\n\tSite = {}\n\tTreatment = {}\
            \n\tYearly fish production (tonnes) = {}\
            \n\tTotal nutrient discharge (kg N per tonne of fish produced) = {}\
            \n\tReference = {}""".format(self.name(),
            self.__site_obj.__str__(2), self.__treatment_obj.__str__(2),
            self.yearly_fish_production('tonne'),
            self.__total_nutrient_discharge, self.__reference)
    
    def Site(self):
        return self.__site_obj
        
    def Treatment(self):
        return self.__treatment_obj
        
    def yearly_fish_production(self, units='kg'):
        return convert_mass(self.__yearly_fish_production, units, self.name())
    
    def total_nutrient_discharge(self):
        return self.__total_nutrient_discharge
        
    def set_total_nutrient_discharge(self, value):
        self.__total_nutrient_discharge = value
        self._sanitize_total_nutrient_discharge()
        
    def _sanitize_total_nutrient_discharge(self):
        assert(type(self.__total_nutrient_discharge) is float)
        if self.__total_nutrient_discharge < 0.0:
            raise InputError(self.__total_nutrient_discharge,
                "Farm's total nutrient discharge should be a positive number")
    
    def _sanitize(self, input_mass_units):
        assert(type(self.__yearly_fish_production) is float)
        if self.__yearly_fish_production < 0.0:
            raise InputError(self.__yearly_fish_production,
                "Farm's yearly fish production should be a positive number")
                
        self._sanitize_total_nutrient_discharge()
        
        # Yearly fish production stored in kg
        self.__yearly_fish_production = convert_mass_to_prog_units(
            self.__yearly_fish_production, input_mass_units, self.name())
                
    def name(self):
        return self.__name
        
    def ece(self, total_max_consented_biomass, loch_obj,
            input_mass_units='tonne', units='umol/L'):
        """Return the equilibrium concentration enhancement (ECE) for nitrogen
        in kg/m^3
        
            ECE = total_nutrient_discharge * total_max_consented_biomass
                / flushing_rate
                
        Arguments:
            total_max_consented_biomass: float; Total maximum consented biomass
            
            loch_obj: Loch object, from which the flushing rate in m^3/yr is
                obtained
            
            input_mass_units: string; input mass units for the total maximum
                consented biomass. default is tonne
            
            units: string; output ECE units. default is umol/L
                
        References:
          - https://consultation.sepa.org.uk/permits/loch-long-salmon-beinn-
                reithe-car-application/user_uploads/210113-bnrt-par-nutrient-
                modelling_redacted.pdf, pp. 12-13
                
          - Gillibrand, P.A., Gubbins, G.J., Greathead, C., Davies, I.M. (2002),
            "Scottish executive locational guidelines for fish farming:
            predicted levels of nutrient enhancement and benthic impact",
            Scottish Fisheries Research Report Number 63 / 2002
        """
        tmcb_kg = convert_mass_to_prog_units(total_max_consented_biomass,
                                             input_mass_units, self.name())
        tmcb_tonne = convert_mass(tmcb_kg, 'tonne', self.name())
        
        flushing_rate = loch_obj.flushing_rate('m^3/yr')
        
        # ECE in kg/m^3
        ece = self.__total_nutrient_discharge*tmcb_tonne/flushing_rate
        
        # Conversion to output units
        u1, u2 = units.split('/')
        f2 = convert_vol(1.0, u2, self.name())
        if "mol" in u1:
            molweight_nitrogen = 0.014
            f1 = convert_mol(1.0, u1, self.name())/molweight_nitrogen
        elif "g" in u1:
            f1 = convert_mass(1.0, u1, self.name())
        else:
            raise InputError(units, "ECE units not accepted in {}".format(
                             self.name()))
        return ece*f1/f2
