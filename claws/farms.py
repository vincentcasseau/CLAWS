#!/usr/bin/env python
"""Derived classes defining Farms."""

# Import modules
import numpy as np
from claws.farm import Farm, total_nutrient_discharge

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

class SalmonFarm(Farm):
    def __init__(self, Site, Treatment, yearly_fish_production=0.0,
                 input_mass_units='kg', reference=""):
        super().__init__(Site=Site, Treatment=Treatment,
                         yearly_fish_production=yearly_fish_production,
                         total_nutrient_discharge=
                             total_nutrient_discharge["salmon"],
                         input_mass_units=input_mass_units, reference=reference)
                         
class HalibutFarm(Farm):
    def __init__(self, Site, Treatment, yearly_fish_production=0.0,
                 input_mass_units='kg', reference=""):
        super().__init__(Site=Site, Treatment=Treatment,
                         yearly_fish_production=yearly_fish_production,
                         total_nutrient_discharge=
                             total_nutrient_discharge["halibut"],
                         input_mass_units=input_mass_units, reference=reference)
                         
class TurbotFarm(Farm):
    def __init__(self, Site, Treatment, yearly_fish_production=0.0,
                 input_mass_units='kg', reference=""):
        super().__init__(Site=Site, Treatment=Treatment,
                         yearly_fish_production=yearly_fish_production,
                         total_nutrient_discharge=
                             total_nutrient_discharge["turbot"],
                         input_mass_units=input_mass_units, reference=reference)
                         
class CodFarm(Farm):
    def __init__(self, Site, Treatment, yearly_fish_production=0.0,
                 input_mass_units='kg', reference=""):
        super().__init__(Site=Site, Treatment=Treatment,
                         yearly_fish_production=yearly_fish_production,
                         total_nutrient_discharge=
                             total_nutrient_discharge["cod"],
                         input_mass_units=input_mass_units, reference=reference)
                         
class HaddockFarm(Farm):
    def __init__(self, Site, Treatment, yearly_fish_production=0.0,
                 input_mass_units='kg', reference=""):
        super().__init__(Site=Site, Treatment=Treatment,
                         yearly_fish_production=yearly_fish_production,
                         total_nutrient_discharge=
                             total_nutrient_discharge["haddock"],
                         input_mass_units=input_mass_units, reference=reference)
