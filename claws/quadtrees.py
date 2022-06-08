#!/usr/bin/env python
"""Derived classes defining different quadtree structures."""

# Import modules
import numpy as np
from claws.quadtree import Quadtree

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

class NoQuadtree(Quadtree):
    def __init__(self):
        super().__init__(is_active=False)
        
class M1Quadtree(Quadtree):
    def __init__(self, max_quadtree_depth=1, min_particles_per_bin=1000):
        super().__init__(is_active=True, max_quadtree_depth=max_quadtree_depth,
                         min_particles_per_bin=min_particles_per_bin)
                         
class M2Quadtree(Quadtree):
    def __init__(self, min_particles_per_bin=1000, leaf_bin_width=np.nan):
        super().__init__(is_active=True,
                         min_particles_per_bin=min_particles_per_bin,
                         leaf_bin_width=leaf_bin_width)
                         
class M3Quadtree(Quadtree):
    def __init__(self, max_quadtree_depth=1, min_particles_per_bin=1000,
                 concentration_target=1.0, input_conc_units='ng/L'):
        super().__init__(is_active=True, max_quadtree_depth=max_quadtree_depth,
                         min_particles_per_bin=min_particles_per_bin,
                         variable_root_bin_width=True,
                         concentration_target=concentration_target,
                         input_conc_units=input_conc_units)
                         
class M4Quadtree(Quadtree):
    def __init__(self, min_particles_per_bin=1000, concentration_target=1.0,
                 leaf_bin_width=np.nan, input_conc_units='ng/L'):
        super().__init__(is_active=True,
                         min_particles_per_bin=min_particles_per_bin,
                         variable_root_bin_width=True,
                         concentration_target=concentration_target,
                         leaf_bin_width=leaf_bin_width,
                         input_conc_units=input_conc_units)
                         
class M5Quadtree(Quadtree):
    def __init__(self, min_particles_per_bin=1000, concentration_target=1.0,
                 leaf_to_seeding_area_ratio=5.0, input_conc_units='ng/L'):
        super().__init__(is_active=True,
                         min_particles_per_bin=min_particles_per_bin,
                         variable_root_bin_width=True,
                         concentration_target=concentration_target,
                         variable_leaf_bin_width=True,
                         leaf_to_seeding_area_ratio=leaf_to_seeding_area_ratio,
                         input_conc_units=input_conc_units)
