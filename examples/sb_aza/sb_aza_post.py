#!/usr/bin/env python
"""Post-processing script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"

# Import modules
import os
import numpy as np
import datetime

import opendrift
from pyproj import Proj

import claws.helpers as helpers
import claws.claws as claws
import claws.opendrift_wrapper as od_wrapper
import claws.postpro as postpro
from claws.treatment import compute_marker_indices
from claws.probes import *
from sb_aza_setup import *


# ---------------------------------------------------------------------------- #
# Computations
# ---------------------------------------------------------------------------- #

# Set the random seed
np.random.seed(seed)

# Sanitise working folder and output file
working_folder = claws.sanitise_working_folder(working_folder)
of = claws.sanitise_output_file(working_folder, outfile, nsimulations)
# Compute marker indices
compute_marker_indices(farms, chemicals)

# Useful variables
z0 = depth_bin_range_meters[1]
dz = abs(depth_bin_range_meters[1] - depth_bin_range_meters[0])
particle_weight_ng = 1.e9*particle_weight_grams

pixelsize_meters = quadtree.update(pixelsize_m=pixelsize_meters, bin_depth=dz,
                                   mass_particle=particle_weight_ng,
                                   seeding_radius=farms[0].Treatment()
                                      .tarpaulin_radius())

root_voxelsize_dm3 = 1.e3*dz*pixelsize_meters**2

# Open the output file with Xarray
oa = opendrift.open_xarray(of)

# Loop over chemical substances
for sp in range(len(chemicals)):
    species_name = chemicals[sp].name()
    if len(chemicals) > 1:
        file_prefix = species_name + "_"
    else:
        file_prefix = ''
    
    # Create array with the right shape using data from the first simulation
    concentration = od_wrapper.get_histogram(outfile=of,
                                             pixelsize_m=pixelsize_meters,
                                             corners=domain_extent, status=0,
                                             z0=z0, dz=dz, marker_index=sp)

    # Compute bounds of the concentration domain                  
    conc_domain_extent = od_wrapper.compute_concentration_domain(
        concentration[0])

    # Sanitise output_options dictionary entries
    claws.sanitise_output_options(concentration)
    cf, lf, tf = claws.get_unit_factors()

    # Extract time data from histogram
    ndt = np.shape(concentration)[0]
    run_duration_sec = float((concentration.time[-1]-concentration.time[0])/1e9)
    time_last_treatment = postpro.find_last_treatment(farms)
    len_time_bins = len(claws.output_options["time_bins"])

    # Create empty arrays
    quadtree_area_over_eqs = np.zeros(shape=(len_time_bins))
    quadtree_area_over_mac = np.zeros(shape=(len_time_bins))

    # Average concentration over all ensembles and seeding locations
    if nsimulations > 1:
        concentration = od_wrapper.average_ensembles(
            h=concentration, nsimulations=nsimulations,
            outfile=working_folder+outfile, pixelsize_m=pixelsize_meters,
            corners=domain_extent, status=0, z0=z0, dz=dz, marker_index=sp)

    # Create sub-histograms to refine regions of high particle density
    # Note that the sign of the histogram bins that present sub-histograms is
    # flipped to keep track of them without storing additional data
    quadtree_conc_lvl, concentration = \
        quadtree.get_histograms(nsimulations=nsimulations,
                                outfile=working_folder+outfile,
                                root_histo=concentration,
                                root_pixelsize_m=pixelsize_meters,
                                status=0, z0=z0, dz=dz, marker_index=sp)

    # Convert concentration of all seedings in ng/L
    pixelsize_level = [pixelsize_meters*(0.5**lvl) for lvl in
        helpers.flatten(quadtree_conc_lvl)]
    voxelsize_dm3 = [root_voxelsize_dm3*(0.25**lvl) for lvl in
        helpers.flatten(quadtree_conc_lvl)]
    particle_to_ngPerL = [particle_weight_ng/vol for vol in voxelsize_dm3]
    concentration = [concentration[i]*particle_to_ngPerL[i] for i in
        range(len(concentration))]

    # Record concentration time series for all probes
    check_probes_validity(probes, conc_domain_extent, concentration)
    record_probes_history(probes, conc_domain_extent, concentration, quadtree,
                          quadtree_conc_lvl)

    # Compute area (in sq meters) where concentration exceeds EQS
    eqs = chemicals[sp].eqs_72hour()
    area_over_eqs = od_wrapper.get_series_area_over_limit(
        np.abs(concentration[0]), pixelsize_m=pixelsize_meters, limit=eqs)

    # Compute area (in sq meters) where concentration exceeds MAC
    mac = chemicals[sp].mac_72hour()
    area_over_mac = od_wrapper.get_series_area_over_limit(
        np.abs(concentration[0]), pixelsize_m=pixelsize_meters, limit=mac)

    # Compute area (in sq meters) where concentration exceeds EQS/MAC for all
    # sub-histograms. Note that because the sign of bins presenting a
    # sub-histogram was flipped, the contribution of depth 'd' is not considered
    # when a depth 'd+1' exists
    if quadtree.is_active():
        offset = 1
        for i, tbin in enumerate(claws.output_options["time_bins"]):
            quadtree_area_over_eqs[i] = od_wrapper.get_series_area_over_limit(
                concentration[0][tbin], pixelsize_m=pixelsize_level[0],
                limit=eqs)
            quadtree_area_over_mac[i] = od_wrapper.get_series_area_over_limit(
                concentration[0][tbin], pixelsize_m=pixelsize_level[0],
                limit=mac)
                    
            for c in np.arange(offset, offset + len(quadtree_conc_lvl[1+i])):
                quadtree_area_over_eqs[i] += od_wrapper.get_series_area_over_limit(
                    concentration[c], pixelsize_m=pixelsize_level[c],
                    limit=eqs)
                quadtree_area_over_mac[i] += od_wrapper.get_series_area_over_limit(
                    concentration[c], pixelsize_m=pixelsize_level[c],
                    limit=mac)
        
            # Increment histogram offset
            offset += len(quadtree_conc_lvl[1+i])

    # Restore the positive sign of the histogram bins that present
    # sub-histograms
    concentration = [np.abs(conc) for conc in concentration]

    # Compute peak concentration
    peakconc = od_wrapper.get_series_peak_concentration(concentration[0])

    # Compute particle vertical distribution
    bar_height, nmax, vdist = od_wrapper.get_series_vertical_distribution(
        outfile=working_folder + outfile, nsimulations=nsimulations, ndt=ndt,
        marker_index=sp, maxdepth=-20.0)

    # Conversion to user-defined output units  
    time_last_treatment *= tf*3600.
    concentration = [cf*conc for conc in concentration]
    peakconc *= cf
    area_over_eqs *= lf**2
    quadtree_area_over_eqs *= lf**2
    area_over_mac *= lf**2 
    quadtree_area_over_mac *= lf**2


# ---------------------------------------------------------------------------- #
# Plots 
# ---------------------------------------------------------------------------- #

    # Plot concentration on terrain map for all output times
    quadtree_peakconc = postpro.plot_concentration_map(
        working_folder + file_prefix, domain_extent, concentration,
        quadtree_conc_lvl, chemicals[sp], quadtree, farms, probes)

    # Plot concentration time series at probe stations
    t = np.linspace(start=0.0, stop=run_duration_sec*tf, num=ndt)
    plot_probes_concentration(probes, working_folder, t, time_last_treatment,
                              chemicals[sp], pixelsize_meters, dz, quadtree,
                              ylabel='{} concentration'.format(species_name),
                              filename=file_prefix + 'probe')

    # Plot time series of peak concentration
    postpro.plot_series_peak_concentration(working_folder + file_prefix, t,
        peakconc, quadtree_peakconc, chemicals[sp], quadtree, pixelsize_meters,
        dz, time_last_treatment, last_treatment_label='Last treatment')

    # Plot time series of area greater than EQS
    postpro.plot_series_area_greater_EQS(working_folder + file_prefix, t,
        area_over_eqs, quadtree_area_over_eqs, chemicals[sp], quadtree,
        pixelsize_meters, dz, time_last_treatment,
        last_treatment_label='Last treatment')
        
    # Plot time series of area greater than MAC
    postpro.plot_series_area_greater_MAC(working_folder + file_prefix, t,
        area_over_mac, quadtree_area_over_mac, chemicals[sp], quadtree,
        pixelsize_meters, dz, time_last_treatment,
        last_treatment_label='Last treatment')

    # Plot vertical distribution bar graph for all output times
    postpro.plot_vertical_distribution_bar_graph(working_folder + file_prefix,
                                                 vdist, bar_height, nmax)

    # Animate concentration on terrain map
    postpro.animate_concentration_terrain(working_folder + file_prefix)

    # Animate particle vertical distribution
    postpro.animate_vertical_distribution(working_folder + file_prefix)
            
    # Print time series to file
    postpro.print_time_series_to_file(working_folder + file_prefix,
                                      [("time", t),
                                       ("peak_concentration", peakconc),
                                       ("area_over_eqs", area_over_eqs)])
