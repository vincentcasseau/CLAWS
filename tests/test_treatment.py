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

from claws.chemicals import *
from claws.treatments import *

if __name__ == "__main__":
    treatment = SeaLiceTreatment(tarpaulin_height=3., tarpaulin_radius=19.,
                                 seeding_times=[0., 3., 6.], nparticles=10000,
                                 Chemicals=Deltamethrin())
    print(treatment)
