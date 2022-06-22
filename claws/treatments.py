#!/usr/bin/env python
"""Derived classes defining Treatments."""

# Import modules
import numpy as np
from claws.treatment import Treatment

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

class NoTreatment(Treatment):
    def __init__(self):
        super().__init__(tarpaulin_height=1.0, tarpaulin_radius=1.0,
                         seeding_times=[], nparticles=[], Chemicals=[])
                         
class Dyeing(Treatment):
    def __init__(self, nparticles, Chemicals, seeding_times=[0.0],
                 input_time_units='h', reference="", name=""):
        super().__init__(tarpaulin_height=1.0, tarpaulin_radius=1.0,
                         seeding_times=seeding_times, nparticles=nparticles,
                         Chemicals=Chemicals, input_time_units=input_time_units,
                         reference=reference, name=name)
                         
class BathMedicine(Treatment):
    def __init__(self, tarpaulin_height, tarpaulin_radius, seeding_times,
                 nparticles, Chemicals, input_len_units='m',
                 input_time_units='h', reference="", name=""):
        super().__init__(tarpaulin_height, tarpaulin_radius, seeding_times,
                         nparticles, Chemicals, input_len_units,
                         input_time_units, reference, name)
