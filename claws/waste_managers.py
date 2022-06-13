#!/usr/bin/env python
"""Derived classes defining WasteManagers."""

# Import modules
import numpy as np
from claws.waste_manager import WasteManager

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

class NoWasteManager(WasteManager):
    def __init__(self):
        super().__init__(feed_requirement=0.0,
                         feed_water_percentage=0.0,
                         feed_waste_percentage=0.0,
                         feed_absorbed_percentage=0.0,
                         feed_carbon_percentage=0.0,
                         faeces_carbon_percentage=0.0)
                         
class WasteFeedManager(WasteManager):
    def __init__(self, feed_requirement=0.0, feed_water_percentage=9.0,
                 feed_waste_percentage=3.0, feed_absorbed_percentage=85.0,
                 feed_carbon_percentage=49.0, faeces_carbon_percentage=30.0,
                 reference="", name=""):
        super().__init__(feed_requirement, feed_water_percentage,
                         feed_waste_percentage, feed_absorbed_percentage,
                         feed_carbon_percentage, faeces_carbon_percentage,
                         reference, name)
