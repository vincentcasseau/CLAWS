#!/usr/bin/env python
"""Derived classes defining Sites."""

# Import modules
import numpy as np
from claws.site import Site

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

class GreatCumbrae(Site):
    def __init__(self):
        super().__init__(longitude=-4.898109, latitude=55.75061)
        
class LittleCumbrae(Site):
    def __init__(self):
        super().__init__(longitude=-4.95608, latitude=55.73220)        
        
class SouthBute(Site):
    def __init__(self):
        super().__init__(longitude=-5.001, latitude=55.738552)   

class Ardentinny(Site):
    def __init__(self):
        super().__init__(longitude=-4.8941, latitude=56.0191)   

class NorthKilbrannan(Site):
    def __init__(self):
        super().__init__(longitude=-5.43519, latitude=55.693)
