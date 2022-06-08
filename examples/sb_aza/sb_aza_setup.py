#!/usr/bin/env python
"""Setup script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"

# Import modules
import numpy as np
import datetime
from pyproj import Proj

import claws.claws as claws
from claws.sites import *
from claws.lochs import *
from claws.chemicals import *
from claws.treatments import *
from claws.farms import *
from claws.probe import Probe
from claws.quadtrees import *


# ------------------------------ FILE OPERATIONS ----------------------------- #
# Working folder (enter './' if current folder)
working_folder = './'
# Selafin inpout file
selafin_file = "r3d_tide-ES_real_gen_new_3.slf"
# Name of the Opendrift outfile minus its .nc extension
outfile = 'sb_aza_full'
# Number of simulations (>1: ensemble averaging)
nsimulations = 1


# ------------------------------ TIME CONTROLS ------------------------------- #
# Simulation start time (date)
start_time = datetime.datetime(2018,4,29,7,0)
# Time-step (seconds)
timestep_seconds = 600
# Total simulation time (hours)
run_duration_hours = 240. # 10 days; 8 days dosage (0,1,2,3,4,5,6,7) then 3 days free flow.
#run_duration_hours = 6.
# Period at which to create sub-histograms and plot concentration maps.
# default is None, in which case the time_bins entry is considered
# If both are None, only the last time step is considered
claws.output_options["output_period"] = 40 # every 4*600 s
# time_bins: integer or list; output time bins, as defined by the user.
# default is None. A value of None means that the variable 'output_period'
# is considered instead 
# ("time_bins" units: "time_units")
claws.output_options["time_bins"] = None


# ------------------------------ DOMAIN & PROJ ------------------------------- #
# Projection and transformation
proj4_params = '+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

# Extent of the domain [lon_min, lon_max, lat_min, lat_max]
domain_extent = [-5.2, -4.8, 55.6725, 55.85]


# ------------------------------- FARM LOCATIONS ----------------------------- #
# Loch name
loch = LochCreran()

# Chemicals used
chemicals = [Azamethiphos(half_life=5.6, Loch=loch, input_time_units='day')]

# Seeding times (hours)
seeding_times = [0., 24., 48., 72., 96., 120., 144., 168.]
#seeding_times = [0.]

# Weight of a single particle (grams)
particle_weight_grams = 2.06e-2

# Farms
farms = [SalmonFarm(GreatCumbrae(),
                    SeaLiceTreatment(tarpaulin_height=3.,
                                     tarpaulin_radius=7.,
                                     seeding_times=seeding_times,
                                     nparticles=10000,
                                     Chemicals=chemicals[0]))]


# ------------------------------ POST-PROCESSING ----------------------------- #
# Probe locations
probes = [Probe(Station(longitude=-4.911, latitude=55.742)),
          Probe(Station(longitude=-4.922, latitude=55.749)),
          Probe(Station(longitude=-4.932, latitude=55.748))]
                      
# Longitude/latitude bin size in meters
pixelsize_meters = 35
# Depth bin range in meters: (max_depth, min_depth)
depth_bin_range_meters = (-3, 0)

# Quadtree structure
quadtree = NoQuadtree()


# ------------------------------ MISCELLANEOUS ------------------------------- #
# Random seed number (None for random or integer >=0 for reproducibility)
seed = None
