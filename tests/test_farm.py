#!/usr/bin/env python
"""Test script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"

from claws.sites import *
from claws.chemicals import *
from claws.treatments import *
from claws.farms import *

if __name__ == "__main__":
    farm = HaddockFarm(SouthBute(),
                       BathMedicine(tarpaulin_height=3.,
                                    tarpaulin_radius=19.,
                                    seeding_times=[0., 3., 6.],
                                    nparticles=10000,
                                    Chemicals=Cypermethrin()),
                       yearly_fish_production=50.,
                       input_mass_units='tonne')
    print(farm)
