#!/usr/bin/env python
"""Test script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "hystrath@gmail.com"
__status__ = "Production"

from claws.sites import *
from claws.lochs import *
from claws.chemicals import *
from claws.treatments import *
from claws.farm import nutrient_enhancement_index
from claws.farms import *

def print_biomass_info(biomass_tonne, header_string):
    ece1 = farm.ece(biomass_tonne, loch, units='kg/m^3')
    ece2 = farm.ece(biomass_tonne, loch, units='ug/L')
    ece3 = farm.ece(biomass_tonne, loch)
    nei = nutrient_enhancement_index(ece3)                       
    
    print("\n" + header_string)
    print("{}'s ECE (kg/m^3) = {:.6e}".format(farm.name(), ece1))
    print("{}'s ECE (ug/L) = {:.6f}".format(farm.name(), ece2))
    print("{}'s ECE (umol/L) = {:.6f}".format(farm.name(), ece3))
    print("{}'s nutrient enhancement index = {:d}".format(farm.name(), nei))

if __name__ == "__main__":
    loch = LochLong()
    
    farm = SalmonFarm(Ardentinny(),
                      SeaLiceTreatment(tarpaulin_height=3.,
                                       tarpaulin_radius=19.,
                                       seeding_times=[0., 3., 6.],
                                       nparticles=10000,
                                       Chemicals=Azamethiphos()),
                      yearly_fish_production=50.,
                      input_mass_units='tonne')
    
    print(loch)
    print(farm)
    
    existing_biomass_tonne = 500.
    print_biomass_info(existing_biomass_tonne, "Existing biomass:")
    
    option1_biomass_tonne = 2000.
    print_biomass_info(option1_biomass_tonne, "Released biomass - Option 1:")
    
    total_max_consented_biomass_tonne = existing_biomass_tonne + \
        option1_biomass_tonne
    print_biomass_info(total_max_consented_biomass_tonne,
                       "Released biomass - Total:")
