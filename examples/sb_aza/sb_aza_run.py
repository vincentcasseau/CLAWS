#!/usr/bin/env python
"""Run script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "hystrath@gmail.com"
__status__ = "Production"

# Import modules
import os
import sys
import numpy as np
from datetime import timedelta

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_telemac_selafin

import claws.claws as claws
from claws.treatment import compute_marker_indices
from sb_aza_setup import *


# ---------------------------------------------------------------------------- #
# OceanDrift parameters 
# ---------------------------------------------------------------------------- #

# Instantiate OpenDrift object
o = OceanDrift(loglevel=20, seed=seed)

# Set OpenDrift configuration
o.set_config('general:coastline_action', 'previous')
o.set_config('drift:vertical_mixing', True)
o.set_config('drift:advection_scheme', 'runge-kutta4')
o.set_config('drift:horizontal_diffusivity', 0.1)
o.set_config('drift:half_life', chemicals[0].half_life())
o.set_config('vertical_mixing:diffusivitymodel', 'environment')
o.set_config('vertical_mixing:timestep', timestep_seconds)
o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.001)

# List properties which can be configured
o.list_configspec()


# ---------------------------------------------------------------------------- #
# Computations
# ---------------------------------------------------------------------------- #

# Set the random seed
np.random.seed(seed)

# Sanitise working folder and output file
working_folder = claws.sanitise_working_folder(working_folder)
of = claws.sanitise_output_file(working_folder, outfile, sys.argv)
# Compute marker indices
compute_marker_indices(farms, chemicals)

# Create selafin reader
selafin = reader_telemac_selafin.Reader(filename=working_folder + selafin_file,
                                        proj4=proj4_params,
                                        start_time=start_time)

print(selafin)

o.add_reader(selafin)

# Seed particles
for farm in farms:
    # Shorthands
    lon = farm.Site().lon()
    lat = farm.Site().lat()
    max_depth = -farm.Treatment().tarpaulin_height()
    rad = farm.Treatment().tarpaulin_radius()
    npart = farm.Treatment().nparticles()
    seeding_times = farm.Treatment().seeding_times()
    ntreatments = len(seeding_times)
    marker = farm.Treatment().marker_index()
    
    # Seed for each treatment
    for t in range(ntreatments):
        zt = max_depth[t]*np.random.rand(npart[t])
        tseed = timedelta(hours=seeding_times[t])
        o.seed_elements(lon, lat, z=zt, radius=rad[t], number=npart[t],
                        origin_marker=marker[t], time=start_time + tseed)

# Run model
o.run(duration=timedelta(hours=run_duration_hours), time_step=timestep_seconds,
      time_step_output=timestep_seconds, outfile=of)
