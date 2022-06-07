#!/usr/bin/env python
"""Post-processing script"""

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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datetime

import cartopy.crs as ccrs
import opendrift
from pyproj import Proj

import claws.helpers as helpers
import claws.claws as claws
import claws.opendrift_wrapper as od_wrapper
import claws.postpro as postpro
from claws.treatment import compute_marker_indices
from claws.farm import nutrient_enhancement_index
from claws.probes import *
from effluent_setup import *


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
    eqs_aze = chemicals[sp].eqs_aze()
    area_over_eqs = od_wrapper.get_series_area_over_limit(
        np.abs(concentration[0]), pixelsize_m=pixelsize_meters, limit=eqs)

    # Compute area (in sq meters) where concentration exceeds MAC
    mac = chemicals[sp].mac_72hour()
    mac_aze = chemicals[sp].mac_aze()
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
    bar_height, zmax, nmax, vdist = od_wrapper.get_series_vertical_distribution(
        outfile=working_folder+outfile, nsimulations=nsimulations, ndt=ndt,
        marker_index=sp, nvbinsmax=20, maxdepth=-3.0)

    # Conversion to user-defined output units  
    time_last_treatment *= tf*3600.
    concentration = [cf*conc for conc in concentration]
    peakconc *= cf
    area_over_eqs *= lf**2
    quadtree_area_over_eqs *= lf**2
    area_over_mac *= lf**2 
    quadtree_area_over_mac *= lf**2

    if np.isfinite(chemicals[sp].eqs_3hour()):
        eqs_3hour = chemicals[sp].eqs_3hour()*cf
    if np.isfinite(chemicals[sp].eqs_72hour()):    
        eqs_72hour = chemicals[sp].eqs_72hour()*cf
    if np.isfinite(chemicals[sp].mac_72hour()):    
        mac_72hour = chemicals[sp].mac_72hour()*cf
    eqs_aze *= lf**2
    mac_aze *= lf**2


# ---------------------------------------------------------------------------- #
# Plots 
# ---------------------------------------------------------------------------- #

    globe = ccrs.Mercator().globe

    # Pre-compute all-time min and max concentrations to use the same scale
    # throughout
    cminmin, cmaxmax, quadtree_peakconc = \
        od_wrapper.get_alltime_min_max_concentrations(
            concentration, quadtree_conc_lvl, quadtree)

    # Plot concentration on terrain map for all output times
    offset = 1
    for i, tbin in enumerate(claws.output_options["time_bins"]):
        fig, ax = postpro.create_terrain(corners=domain_extent,
                                         seeding_locations=farms,
                                         display_probe_markers=True,
                                         probes=probes)

        # Fix the scale for all plots
        axmin, axmax = postpro.round_logscale_axis([cminmin, cmaxmax])
        
        # Superimpose histograms on top of the root histogram
        conc = concentration[0][tbin].transpose()
        c = ax.pcolormesh(concentration[0].lon_bin, concentration[0].lat_bin,
                          conc.where(conc>0), cmap='jet',
                          norm=colors.LogNorm(vmin=axmin, vmax=axmax),
                          transform=ccrs.PlateCarree(globe=globe))
                          
        if quadtree.is_active():
            for hi in np.arange(offset, offset + len(quadtree_conc_lvl[1+i])):
                conc = concentration[hi].transpose()
                c = ax.pcolormesh(concentration[hi].lon_bin,
                                  concentration[hi].lat_bin,
                                  conc.where(conc>0), cmap='jet',
                                  norm=colors.LogNorm(vmin=axmin, vmax=axmax),
                                  transform=ccrs.PlateCarree(globe=globe))
                              
            # Increment histogram offset
            offset += len(quadtree_conc_lvl[1+i])

        # Set colorbar and save
        cbar = fig.colorbar(c, ax=ax, orientation="horizontal", shrink=0.65,
                            pad=0.075, aspect=25, extend='min')
        cbar.set_label('{} concentration {}'.format(species_name,
            claws.output_options["concentration_units"]))
        plt.title(np.datetime64(concentration[0].time[tbin].item(),
            'ns').astype('datetime64[us]').astype(datetime.datetime).strftime(
            '%d %B %Y, %I:%M:%S %p'))   
        plt.savefig(working_folder + file_prefix +
            'concentration_{:04d}.png'.format(tbin))
        plt.close()
            
    # Plot concentration time series at probe stations
    t = np.linspace(start=0.0, stop=run_duration_sec*tf, num=ndt)
    plot_probes_concentration(probes, working_folder, t, time_last_treatment,
                              chemicals[sp], pixelsize_meters, dz, quadtree,
                              last_treatment_label='Last spillage',
                              ylabel='{} concentration'.format(species_name),
                              filename=file_prefix+'probe')

    # Plot time series of peak concentration
    fig, ax = plt.subplots(1)
    psize_str = ('%f'%pixelsize_meters).rstrip('0').rstrip('.')
    dz_str = ('%f'%dz).rstrip('0').rstrip('.')
    plt.plot(t, peakconc, color='black', linestyle='-', lw=1,
             label='{0} m x {0} m x {1} m'.format(psize_str, dz_str))
    if quadtree.is_active():
        plt.scatter(claws.output_options["time_bins"]*tf,
                    quadtree_peakconc, color='blue', marker='+',
                    label='quadtree (max depth: {})'.format(
                        quadtree.max_quadtree_depth()))         
    if np.isfinite(chemicals[sp].eqs_3hour()):
        plt.axhline(y=eqs_3hour, color='green', linestyle='dotted',
                    label='3-hour EQS')
    if np.isfinite(chemicals[sp].eqs_72hour()):
        plt.axhline(y=eqs_72hour, color='green', linestyle='--',
                    label='72-hour EQS')
    if np.isfinite(chemicals[sp].mac_72hour()):
        plt.axhline(y=mac_72hour, color='red', linestyle='--',
                    label='72-hour MAC')
    ax.set_yscale('log')
    postpro.round_logscale_yaxis(ax)
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    postpro.print_last_treatment(plt, ax, time_last_treatment,
                                 label='Last spillage')
    plt.xlabel('Time {}'.format(claws.output_options["time_units"]))
    plt.ylabel('Peak {} concentration {}'.format(species_name,
        claws.output_options["concentration_units"]))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=2)
    plt.savefig(working_folder + file_prefix + 'peak_concentration.png')
    plt.close()

    # Plot time series of area greater than EQS
    fig, ax = plt.subplots(1)
    plt.plot(t, area_over_eqs, color='black', linestyle='-', lw=1,
             label='{0} m x {0} m x {1} m'.format(psize_str, dz_str))
    if quadtree.is_active():
        plt.scatter(claws.output_options["time_bins"]*tf,
                    quadtree_area_over_eqs,
                    color='blue', marker='+',
                    label='quadtree (max depth: {})'.format(
                        quadtree.max_quadtree_depth()))  
    plt.axhline(y=eqs_aze, color='g', linestyle='--', label='72-hour EQS AZE')
    ax.set_yscale('log')    
    postpro.round_logscale_yaxis(ax)  
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    postpro.print_last_treatment(plt, ax, time_last_treatment,
                                 label='Last spillage')
    plt.xlabel('Time {}'.format(claws.output_options["time_units"]))
    plt.ylabel('Area > EQS ({}$^2$)'.format(
        claws.output_options["length_units"]))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=2)
    plt.savefig(working_folder + file_prefix + 'area_over_EQS.png')
    plt.close()

    # Plot time series of area greater than MAC
    fig, ax = plt.subplots(1)
    plt.plot(t, area_over_mac, color='black', linestyle='-', lw=1,
             label='{0} m x {0} m x {1} m'.format(psize_str, dz_str))
    if quadtree.is_active():
        plt.scatter(claws.output_options["time_bins"]*tf,
                    quadtree_area_over_mac,
                    color='blue', marker='+',
                    label='quadtree (max depth: {})'.format(
                        quadtree.max_quadtree_depth())) 
    plt.axhline(y=mac_aze, color='r', linestyle='--', label='72-hour MAC AZE')
    ax.set_yscale('log')    
    postpro.round_logscale_yaxis(ax)  
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    postpro.print_last_treatment(plt, ax, time_last_treatment,
                                 label='Last spillage')
    plt.xlabel('Time {}'.format(claws.output_options["time_units"]))
    plt.ylabel('Area > MAC ({}$^2$)'.format(
        claws.output_options["length_units"]))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=2)
    plt.savefig(working_folder + file_prefix + 'area_over_MAC.png')
    plt.close()

    # Plot vertical distribution bar graph for all output times
    for tbin in claws.output_options["time_bins"]:
        arr = np.array([100.*v[tbin].item() for v in vdist])
        len_arr = len(arr)
        zbin = np.linspace(start=0, stop=-bar_height*(len_arr-1), num=len_arr)
        
        fig, ax = plt.subplots(1)
        plt.barh(zbin, arr, height=-bar_height, facecolor='navy', alpha=0.35,
                 lw=1, edgecolor='k', align='edge')
        nyticks, mltp = postpro.compute_nticks(len_arr)
        plt.yticks(np.linspace(start=0, stop=-bar_height*mltp*(nyticks-1),
                               num=nyticks))
        plt.xlabel('Number of particles (%)')
        plt.ylabel('Depth (m)')
        plt.xlim([0, min(100, round(100.*nmax,0)//5*5.+5.)])
        plt.title(np.datetime64(vdist[0].time[tbin].item(), 'ns').astype(
                  'datetime64[us]').astype(datetime.datetime).strftime(
                  '%d %B %Y, %I:%M:%S %p')) 
        plt.savefig(working_folder + file_prefix + 'vertical_distribution'\
                        '_{:04d}.png'.format(tbin))
        plt.close()

    # Animate concentration on terrain map
    postpro.animate_concentration_terrain(working_folder + file_prefix)

    # Animate particle vertical distribution
    postpro.animate_vertical_distribution(working_folder + file_prefix)
            
    # Print time series to file
    postpro.print_time_series_to_file(working_folder + file_prefix,
                                      [("time", t),
                                       ("peak_concentration", peakconc),
                                       ("area_over_eqs", area_over_eqs),
                                       ("area_over_mac", area_over_mac),
                                       ("time_bins", claws.output_options["time_bins"]*tf),
                                       ("quadtree_area_over_eqs", quadtree_area_over_eqs),
                                       ("quadtree_peak_concentration", quadtree_peakconc)])
