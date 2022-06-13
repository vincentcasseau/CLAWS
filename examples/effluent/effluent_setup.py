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
from claws.site import *
from claws.lochs import *
from claws.chemical import *
from claws.treatments import *
from claws.farms import *
from claws.probe import Probe
from claws.quadtrees import *


# ------------------------------ FILE OPERATIONS ----------------------------- #
# Working folder (enter './' if current folder)
working_folder = './'
# Selafin inpout file
selafin_file = "r3d_effluent_k_epsilon.slf"
# Name of the Opendrift outfile minus its .nc extension
outfile = 'effluent'
# Number of simulations (>1: ensemble averaging)
nsimulations = 1


# ------------------------------ TIME CONTROLS ------------------------------- #
# Simulation start time (date)
start_time = datetime.datetime(2018,4,22,7,0,0)
# Time-step (seconds)
timestep_seconds = 1
# Total simulation time (hours)
run_duration_hours = 0.2
# Time units to show on time series plots
claws.output_options["time_units"] = "min"
# Period at which to create sub-histograms and plot concentration maps.
# default is None, in which case the time_bins entry is considered
# If both are None, only the last time step is considered
claws.output_options["output_period"] = None
# time_bins: integer or list; output time bins, as defined by the user.
# default is None. A value of None means that the variable 'output_period'
# is considered instead 
# ("time_bins" units: "time_units")
claws.output_options["time_bins"] = [[30 + i*60, 45 + i*60] for i in range(12)]


# ------------------------------ DOMAIN & PROJ ------------------------------- #
# Projection and transformation
proj4_params = '+proj=ortho +ellps=WGS84'

# Extent of the domain [lon_min, lon_max, lat_min, lat_max]
xy2lonlat = Proj(proj4_params, preserve_units=False)
lon_min, lat_min = xy2lonlat(0., 0., inverse=True)
lon_max, lat_max = xy2lonlat(20.8, 40., inverse=True)
domain_extent = [lon_min, lon_max, lat_min, lat_max]


# ------------------------------- FARM LOCATIONS ----------------------------- #
# Chemicals used
chemicals = [ChemicalSubstance(name='Effluent', eqs_72hr=0.2, mac_72hr=0.5,
                               input_conc_units='ug/L', input_area_units='m^2')]

# Weight of a particle (grams)
particle_weight_grams = 1e-6

# Seeding times (min)
seeding_times = [float(i) for i in range(12)]

# Farms
farms = [SalmonFarm(Site(x=20.8, y=10.2, proj_params=proj4_params),
                    BathMedicine(tarpaulin_height=1.,
                                 tarpaulin_radius=0.01,
                                 seeding_times=seeding_times,
                                 nparticles=500,
                                 Chemicals=chemicals[0],
                                 input_time_units='min'))]


# ---------------------------------- PHYSICS --------------------------------- #
# Horizontal diffusivity (m^2/s)
horizontal_diffusivity = 0.001
# Vertical diffusivity (m^2/s)
vertical_diffusivity = 0.001


# ------------------------------ POST-PROCESSING ----------------------------- #
# Probe locations
probes = [Probe(x=19.0, y=21.0, proj_params=proj4_params),
          Probe(x=19.0, y=25.0, proj_params=proj4_params),
          Probe(x=19.0, y=29.0, proj_params=proj4_params),
          Probe(x=19.0, y=35.0, proj_params=proj4_params),
          Probe(x=19.0, y=39.9, proj_params=proj4_params),
          Probe(x=20.15, y=10.15, proj_params=proj4_params),
          Probe(x=0.01, y=39.9, proj_params=proj4_params)]
                      
# Longitude/latitude bin size in meters
pixelsize_meters = 0.4
# Depth bin range in meters: (max_depth, min_depth)
depth_bin_range_meters = (-0.5, 0)

# Quadtree structure
quadtree = M1Quadtree(min_particles_per_bin=20)

# Length and concentration units to be used in plots & animations
claws.output_options["length_units"] = "m"
claws.output_options["concentration_units"] = "ug/L"


# ------------------------------ MISCELLANEOUS ------------------------------- #
# Random seed number (None for random or integer >=0 for reproducibility)
seed = None
