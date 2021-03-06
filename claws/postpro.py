#!/usr/bin/env python
"""Post-processing functions.
"""

# Import modules
import os
import numpy as np
import xarray as xr
import datetime
import pyproj
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle, Polygon, Circle
import matplotlib.transforms as mtransforms

from opendrift.readers import reader_global_landmask
import cartopy.crs as ccrs
from cartopy.crs import CRS, Globe
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature

from claws.helpers import *
from claws.claws import output_options, animation, get_unit_factors
from claws.opendrift_wrapper import (get_alltime_min_max_concentrations,
                                     get_nparticles_in_polygon)

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Functions 
# ---------------------------------------------------------------------------- #

def find_last_treatment(seeding_locations):
    """Find the time of latest treatment over all seeding locations.
    
    Arguments:
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
    """
    time_last_treatment = -np.inf
    for loc in seeding_locations:
        if len(loc.Treatment().seeding_times()) > 0:
            loc_last_treatment = np.max(loc.Treatment().seeding_times())
            time_last_treatment = np.max([time_last_treatment,
                                          loc_last_treatment])
    return time_last_treatment
    
def plot_selafin_map(working_folder, corners, loch_obj, seeding_locations,
                        probes_obj, proj4_params, selafin_folder=''):
    """Plot selafin data on terrain map
    
    Arguments:
        working_folder: string; working folder
        
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the terrain map
        
        loch_obj: Loch object as implemented in claws.Loch
        
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
        
        probes_obj: list of Probe objects as implemented in claws.Probe
        
        proj4_params: projection parameters
        
        selafin_folder: selafin data input folder, containing xarrays
            (engine: zarr). default is empty string
    """
   
    if not (selafin_folder and os.path.isdir(selafin_folder)):
        return
    
    globe = ccrs.Mercator().globe
    proj = pyproj.Proj(proj4_params)

    # Load selafin data
    arr = xr.open_dataset(selafin_folder, engine='zarr')
    x, y = arr.node_x, arr.node_y
    for i in range(len(x)):
        x[i], y[i] = proj(x[i], y[i], inverse=True)
    
    selafin_vars = {
      'altitude': ['m', 'Bathymetry', 'bathymetry', 'ocean'],
      'sea_water_salinity': ['PSU', 'Sea water salinity', 'salinity', 'summer'],
      'sea_water_temperature': ['??C', 'Sea water temperature', 'temperature', 'jet'],
      'turbulent_kinetic_energy': ['J/kg', 'Turbulent kinetic energy', 'tke', 'jet'],
      'x_sea_water_velocity': ['m/s', 'Velocity, Vx', 'vx', 'jet'],
      'y_sea_water_velocity': ['m/s', 'Velocity, Vy', 'vy', 'jet'],
      'upward_sea_water_velocity': ['m/s', 'Velocity, Vz', 'vz', 'jet'],}
      
    for var in arr.keys():
        # Skip variables that aren't listed in the dictionary above
        if var not in selafin_vars.keys():
            continue
        
        var_units = selafin_vars[var][0]
        var_legend = '{} ({})'.format(selafin_vars[var][1], var_units) 
        var_cmap = selafin_vars[var][3]
        
        for tbin in range(len(arr[var].time)):
            var_time = arr[var].time[tbin].item()
            if var != "altitude":
                var_title = np.datetime64(var_time, 'ns').astype(
                    'datetime64[us]').astype(datetime.datetime).strftime(
                    '%d %B %Y, %I:%M:%S %p')
                var_filename = '{}_{:04d}.png'.format(selafin_vars[var][2],
                                                      tbin)
            else:
                var_title = ' '
                var_filename = '{}.png'.format(selafin_vars[var][2])
            z = arr[var][tbin][0] # first of the 8 layers TODO?
            
            # A contour plot of irregularly spaced data coordinates via
            # interpolation on a grid. Create grid values first
            ngridx = 1000
            ngridy = 1000
            xi = np.linspace(corners[0], corners[1], ngridx)
            yi = np.linspace(corners[2], corners[3], ngridy)

            # Linearly interpolate the data (x, y) on a grid defined by (xi, yi)
            triang = tri.Triangulation(x, y)
            interpolator = tri.LinearTriInterpolator(triang, z)
            Xi, Yi = np.meshgrid(xi, yi)
            zi = interpolator(Xi, Yi)
            if var == 'altitude':
                zi = np.minimum(zi, 0.0)

            # Plot concentration on terrain map for all output times
            fig, ax = create_terrain(corners=corners, loch_obj=loch_obj,
                                     seeding_locations=seeding_locations,
                                     display_probe_markers=True,
                                     probes=probes_obj)

            # Draw contours
            cntr = ax.contourf(xi, yi, zi, levels=20, cmap=var_cmap,
                               transform=ccrs.PlateCarree(globe=globe))
                              
            # Set colorbar and save
            cbar = fig.colorbar(cntr, ax=ax, orientation="horizontal",
                                shrink=0.65, pad=0.075, aspect=25)
            cbar.set_label(var_legend)
            plt.suptitle(var_title, y=0.97) 
            plt.savefig(working_folder + var_filename)
            plt.close()
            break
    
def plot_bathymetry_map(working_folder, corners, loch_obj, seeding_locations,
                        probes_obj, proj4_params, bathymetry_file=''):
    """Plot bathymetry on terrain map
    
    Arguments:
        working_folder: string; working folder
        
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the terrain map
        
        loch_obj: Loch object as implemented in claws.Loch
        
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
        
        probes_obj: list of Probe objects as implemented in claws.Probe
        
        proj4_params: projection parameters
        
        bathymetry_file: bathymetry input file in XYZ format.
            default is empty string
    """
   
    if not (bathymetry_file and os.path.isfile(bathymetry_file)):
        return
    
    globe = ccrs.Mercator().globe
    proj = pyproj.Proj(proj4_params)

    # Load bathymetry data
    a = np.loadtxt(fname=bathymetry_file, comments=['#', ':'])
    for p in range(np.shape(a)[0]):
         a[p][0], a[p][1] = proj(a[p][0], a[p][1], inverse=True)
         a[p][2] = min(a[p][2], 0.0)
    x, y, z = a.transpose()
    
    # A contour plot of irregularly spaced data coordinates via interpolation on
    # a grid. Create grid values first
    ngridx = 1000
    ngridy = 1000
    xi = np.linspace(corners[0], corners[1], ngridx)
    yi = np.linspace(corners[2], corners[3], ngridy)

    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    # Plot concentration on terrain map for all output times
    fig, ax = create_terrain(corners=corners, loch_obj=loch_obj,
                             seeding_locations=seeding_locations,
                             display_probe_markers=True,
                             probes=probes_obj)

    # Draw contours
    cntr = ax.contourf(xi, yi, zi, levels=20, cmap="ocean",
                       transform=ccrs.PlateCarree(globe=globe))
                      
    # Set colorbar and save
    cbar = fig.colorbar(cntr, ax=ax, orientation="horizontal", shrink=0.65,
                        pad=0.075, aspect=25)
    cbar.set_label('Bathymetry (m)')
    plt.suptitle(' ', y=0.97) 
    plt.savefig(working_folder + 'bathymetry_raw.png')
    plt.close()
    
def plot_concentration_map(working_folder, corners, concentration,
                           quadtree_conc_lvl, loch_obj, chemical_obj,
                           quadtree_obj, seeding_locations, probes_obj,
                           normalised=False):
    """Plot concentration on terrain map for all output times
    
    Arguments:
        working_folder: string; working folder
        
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the terrain map
        
        concentration: list of Xarray histograms; as generated by Opendrift
            get_histogram function from the data stored in the .nc outfile
            and the get_quadtree_histograms function
        
        quadtree_conc_lvl: nested python list; output of the
            get_quadtree_histograms function. Indicates the number of
            sub-histograms per output time bins and their depth level.
        
        loch_obj: Loch object as implemented in claws.Loch
        
        chemical_obj: ChemicalSubstance object as implemented in
            claws.ChemicalSubstance
        
        quadtree_obj: Quadtree object as implemented in
            claws.quadtree
        
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
        
        probes_obj: list of Probe objects as implemented in claws.Probe
        
        normalised: bool; whether the concentration field is normalised or not.
            If so, no there won't be units on the contour map. default is False
    """
   
    species_name = chemical_obj.name()
    globe = ccrs.Mercator().globe

    # Pre-compute all-time min and max concentrations to use the same scale
    # throughout
    cminmin, cmaxmax, quadtree_peakconc = \
        get_alltime_min_max_concentrations(
            concentration, quadtree_conc_lvl, quadtree_obj)

    # Plot concentration on terrain map for all output times
    offset = 1
    for i, tbin in enumerate(output_options["time_bins"]):
        fig, ax = create_terrain(corners=corners, loch_obj=loch_obj,
                                 seeding_locations=seeding_locations,
                                 display_probe_markers=True,
                                 probes=probes_obj)

        # Fix the scale for all plots
        axmin, axmax = round_logscale_axis([cminmin, cmaxmax])
        
        # Superimpose histograms on top of the root histogram
        conc = concentration[0][tbin].transpose()
        c = ax.pcolormesh(concentration[0].lon_bin, concentration[0].lat_bin,
                          conc.where(conc>0), cmap='jet',
                          norm=colors.LogNorm(vmin=axmin, vmax=axmax),
                          transform=ccrs.PlateCarree(globe=globe))
                          
        if quadtree_obj.is_active():
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
        if normalised:
            prefix_str = 'Normalised '
            suffix_str = 'density'
        else:
            prefix_str = ''
            suffix_str = 'concentration ({})'.format(
                output_options["concentration_units"])
        cbar.set_label('{}{} {}'.format(prefix_str, species_name, suffix_str))
        plt.suptitle(np.datetime64(concentration[0].time[tbin].item(),
            'ns').astype('datetime64[us]').astype(datetime.datetime).strftime(
            '%d %B %Y, %I:%M:%S %p'), y=0.97)   
        plt.savefig(working_folder + 'concentration_{:04d}.png'.format(tbin))
        plt.close()
    return quadtree_peakconc
    
def create_terrain(corners, loch_obj, seeding_locations, tile_style='',
                   display_locations=False, abbreviate_locations=True,
                   display_probe_markers=False, probes=None,
                   custom_crs=ccrs.Geodetic(),
                   globe=None):
    """Create a terrain map using a Stamen tile of size domain_extent. Seeding
    sources are labelled on the map.
    
    Arguments:
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the terrain map
            
        loch_obj: Loch object as implemented in claws.Loch
        
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
            
        tile_style: string; terrain rendering style. default is Opendrift's
            style (no tile). For further details about the available options,
            see: http://maps.stamen.com/#watercolor/12/37.7706/-122.3782
            Enter '' (blank string) to not use Stamen tiles  
            
        display_locations: bool; whether or not to print the locations' names
            on the map. default is False
            
        abbreviate_locations: bool; whether or not to abbreviate seeding
            locations on the map. If so, the first letter of each word will be
            printed. default is True
            
        display_probe_markers: bool; whether or not to print probe markers on
            the map. default is False
        
        probes: list of Probe objects as implemented in
            claws.Probe
        
        custom_crs: CRS; default is Geodetic
        
        globe: a CRS globe; default is None
    """
    # Store a backup of the original lon/lat min/max values
    lon_min_ini, lon_max_ini, lat_min_ini, lat_max_ini = corners
    # Add small margins (default 2%) around the domain's corners
    lon_min, lon_max, lat_min, lat_max = corners
    ratio = dlon_to_dlat_ratio(lat_min, lat_max)
    dlon = lon_max - lon_min
    dlat = lat_max - lat_min
    
    offset_pct = 0.02
    lon_offset = dlon*offset_pct
    lat_offset = dlat*offset_pct*ratio
    domain_w_margins = [lon_min - lon_offset, lon_max + lon_offset,
                        lat_min - lat_offset, lat_max + lat_offset]
    lon_min, lon_max, lat_min, lat_max = domain_w_margins
    
    # Force an aspect ratio of at least 1.28 for best render
    target_aspect_ratio = 1.28
    dlon = lon_max - lon_min
    dlat = lat_max - lat_min
    lon_mid = lon_min + dlon/2.
    lat_mid = lat_min + dlat/2.
    ratio = dlon_to_dlat_ratio(lat_min, lat_max)
    aspect_ratio = dlon/(dlat*ratio)
    if aspect_ratio < target_aspect_ratio:
        # Increase longitude range (lon_mid and ratio are unchanged)
        dlon = dlat*ratio*target_aspect_ratio
        lon_min = lon_mid - dlon/2.
        lon_max = lon_mid + dlon/2.
        domain_w_margins[0] = lon_min
        domain_w_margins[1] = lon_max
    else:
        # Increase latitude range (lat_mid and ratio are unchanged)
        dlat = dlon/(ratio*target_aspect_ratio)
        lat_min = lat_mid - dlat/2.
        lat_max = lat_mid + dlat/2.
        domain_w_margins[2] = lat_min
        domain_w_margins[3] = lat_max
    # Create figure
    figsize = [9., 6.8]
    fig = plt.figure(figsize=figsize)
    
    # Define axes location
    width_ratios = [5, 1]
    gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios)
    gsmap = gridspec.GridSpec(5, 2, width_ratios=width_ratios)
            
    # Create a Stamen instance
    available_tile_styles = ['terrain', 'terrain-background', 'watercolor',
                             'toner-lite', 'toner-background']
    if tile_style in available_tile_styles:
        tiler = cimgt.Stamen(tile_style)
        # Create a GeoAxes in the tile's projection
        ax = fig.add_subplot(gs[0], projection=tiler.crs)
        globe = tiler.crs.globe
    else:
        ax = fig.add_subplot(gs[0], projection=ccrs.Mercator())
        globe = ccrs.Mercator().globe
    ax.set_extent(domain_w_margins, crs=custom_crs)
    
    # Create a second subplot to contain compass and legend
    ax2 = fig.add_subplot(gs[1])
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, width_ratios[0]])
    ax2.axis('off')
    
    # Draw gridlines on the main map
    gl = ax.gridlines(ccrs.PlateCarree(globe=globe), draw_labels=True,
                      linewidth=0.5, color='grey', alpha=0.75, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    
    # Shift ticks labels away from the axes
    gl.xpadding = 10
    gl.ypadding = 10
    
    # Draw canvas to extract x- and y-ticks locations
    fig.canvas.draw()
    xsegs = gl.xline_artists[0].get_segments()
    ysegs = gl.yline_artists[0].get_segments()
    xticks = [xseg[0,0] for xseg in xsegs]
    yticks = [yseg[0,1] for yseg in ysegs]
    
    # Remove non-visible ticks
    xticks = [t for t in xticks if t>=lon_min and t<=lon_max]
    yticks = [t for t in yticks if t>=lat_min and t<=lat_max]
    
    if tile_style in available_tile_styles:
        # Add the Stamen data at the maximum zoom level of 11
        ax.add_image(tiler, 11)
    else:
        # Add land and shorelines as in OpenDrift
        reader_global_landmask.plot_land(ax, domain_w_margins[0],
                                         domain_w_margins[2],
                                         domain_w_margins[1],
                                         domain_w_margins[3], False,
                                         lscale='full', globe=globe)
    
    # Add frame to delineate the domain extent
    ax.plot([lon_min_ini, lon_min_ini], [lat_min_ini, lat_max_ini],
             color='grey', lw=0.5, alpha=0.5, transform=custom_crs, zorder=3)
    ax.plot([lon_max_ini, lon_max_ini], [lat_min_ini, lat_max_ini],
             color='grey', lw=0.5, alpha=0.5, transform=custom_crs, zorder=3)
    ax.plot([lon_min_ini, lon_max_ini], [lat_min_ini, lat_min_ini],
             color='grey', lw=0.5, alpha=0.5, transform=custom_crs, zorder=3)
    ax.plot([lon_min_ini, lon_max_ini], [lat_max_ini, lat_max_ini],
             color='grey', lw=0.5, alpha=0.5, transform=custom_crs, zorder=3)
    
    # Blur domain area that is outside the original domain
    rectangle = Rectangle((lon_min, lat_min),
                          width=lon_max - lon_min,
                          height=lat_min_ini - lat_min,
                          facecolor='white', alpha=0.65, edgecolor='none',
                          clip_on=False, transform=custom_crs, zorder=2)
    ax.add_patch(rectangle)
    rectangle = Rectangle((lon_min, lat_max_ini),
                          width=lon_max - lon_min,
                          height=lat_max - lat_max_ini,
                          facecolor='white', alpha=0.65, edgecolor='none',
                          clip_on=False, transform=custom_crs, zorder=2)
    ax.add_patch(rectangle)
    rectangle = Rectangle((lon_min, lat_min_ini),
                          width=lon_min_ini - lon_min,
                          height=lat_max_ini - lat_min_ini,
                          facecolor='white', alpha=0.65, edgecolor='none',
                          clip_on=False, transform=custom_crs, zorder=2)
    ax.add_patch(rectangle)
    rectangle = Rectangle((lon_max_ini, lat_min_ini),
                          width=lon_max - lon_max_ini,
                          height=lat_max_ini - lat_min_ini,
                          facecolor='white', alpha=0.65, edgecolor='none',
                          clip_on=False, transform=custom_crs, zorder=2)
    ax.add_patch(rectangle)
             
    # Use the cartopy interface to create a matplotlib transform object for the
    # Geodetic coordinate system. 
    geodetic_transform = custom_crs._as_mpl_transform(ax)
    
    if display_probe_markers:
        if probes is not None:
            for P in probes:    
                # Print probe marker: '+'
                ax.plot(P.lon(), P.lat(), marker='+', color='black',
                        markersize=4, alpha=0.7, transform=custom_crs)
    
    # Add markers and text at seeding locations
    for i, loc in enumerate(seeding_locations):
        # Skip locations where no particles are seeded
        if len(loc.Treatment().seeding_times()) == 0: continue
        # Shorthands
        lon = loc.Site().lon()
        lat = loc.Site().lat()
        # Skip locations that fall outside the domain bounds
        if lon > lon_max or lon < lon_min:
            continue
        if lat > lat_max or lat < lat_min:
            continue
        # Seeding location marker
        ax.plot(lon, lat, marker='o', color='black', markersize=4, alpha=0.7,
                transform=custom_crs)
        # Seeding location name
        if display_locations:
            if abbreviate_locations:
                loc_text = ''.join(word[0] for word in loc.split())
            else:
                loc_text = loc
            # Use matplotlib's offset_copy function to define a coordinate
            # system which translates the text by x pixels to the left and y
            # pixels up.
            if i%4 == 0:
                va = 'center'
                ha = 'left'
                nx = 15
                ny = 0
            elif i%4 == 1:
                va = 'bottom'
                ha = 'center'
                nx = 0
                ny = -25
            elif i%4 == 2:
                va = 'center'
                ha = 'right'
                nx = -15
                ny = 0
            else:
                va = 'top'
                ha = 'center'
                nx = 0
                ny = 25
                
            text_transform = mtransforms.offset_copy(geodetic_transform,
                                                     units='dots', x=nx, y=ny)
            ax.text(lon, lat, loc_text,
                     verticalalignment=va, horizontalalignment=ha,
                     transform=text_transform, fontsize=9,
                     bbox=dict(facecolor='sandybrown', alpha=0.5,
                               boxstyle='round'))

    # Add frame to delineate the domain extent
    draw_frame(ax, domain_w_margins, xticks, yticks, custom_crs)
    
    # Add scale
    draw_scale(ax, domain_w_margins)
    
    # Add compass
    draw_compass(ax2)
    
    # Add legend
    draw_legend(ax2, seeding_locations, probes)
    
    # Add minimap
    draw_minimap(fig, gsmap[1,1], seeding_locations, loch_obj.name())
    
    # Activate figure zooming and maximise figure window size
    fig.canvas.draw()
    fig.set_tight_layout(True)
    try:  
        mng = plt.get_current_fig_manager()
        mng.toolbar.zoom()
        mng.resize(*mng.window.maxsize())
    except:
        pass
    
    return [fig, ax]
    
def draw_frame(ax, corners, lon_breaks, lat_breaks, custom_crs=ccrs.Geodetic(),
               width_pct=0.01):
    """Draw a thick white and black frame around the terrain map
    
    Arguments:
        ax: matplotlib axis; axis object as returned by plt.subplots
        
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the terrain map
        
        lon_breaks: numpy array; list of xticks, ie, longitudinal positions
            where to break (switch colors) the frame
        
        lat_breaks: numpy array; list of yticks, ie, latitudinal positions
            where to break (switch colors) the frame
        
        custom_crs: CRS; default is Geodetic
        
        width_pct: float; frame thickness expressed as a percentage of the
            terrain map width. default is 1%
    """ 
    # Compute frame width and height from total longitude/latitude range
    # and multiply by the width percentage
    ratio = dlon_to_dlat_ratio(corners[2], corners[3])
    frame_width = max((corners[1] - corners[0]),
                      (corners[3] - corners[2])*ratio)*width_pct
    frame_height = frame_width/ratio                  
    
    # Insert min/max lon/lat values to the list of ticks if absent
    if abs(lon_breaks[0] - corners[0]) > 1e-6:
        lon_breaks = np.insert(lon_breaks, 0, corners[0])
    if abs(lon_breaks[-1] - corners[1]) > 1e-6:
        lon_breaks = np.append(lon_breaks, corners[1])
    if abs(lat_breaks[0] - corners[2]) > 1e-6:
        lat_breaks = np.insert(lat_breaks, 0, corners[2])
    if abs(lat_breaks[-1] - corners[3]) > 1e-6:
        lat_breaks = np.append(lat_breaks, corners[3])
    
    # Draw frame (rectangles: bar, polygons: corners)
    colour = 'black'
    polygon = Polygon([(lon_breaks[0] - frame_width,
                        lat_breaks[0] - frame_height),
                       (lon_breaks[1], lat_breaks[0] - frame_height),
                       (lon_breaks[1], lat_breaks[0]),
                       (lon_breaks[0], lat_breaks[0])],
                       facecolor=colour, edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
    
    polygon = Polygon([(lon_breaks[0] - frame_width,
                        lat_breaks[0] - frame_height),
                       (lon_breaks[0], lat_breaks[0]),
                       (lon_breaks[0], lat_breaks[1]),
                       (lon_breaks[0] - frame_width, lat_breaks[1])],
                       facecolor='white', edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
    
    # Bottom bar
    for i, (from_lon, to_lon) in enumerate(list(zip(lon_breaks[1:-2],
                                                    lon_breaks[2:-1]))):
        colour = 'white' if colour == 'black' else 'black'
        rectangle = Rectangle((from_lon, lat_breaks[0] - frame_height),
                              width=to_lon - from_lon,
                              height=frame_height, facecolor=colour,
                              edgecolor='black', clip_on=False, ls='solid',
                              lw=1, transform=custom_crs, zorder=4)
        ax.add_patch(rectangle)
        
    colour = 'white' if colour == 'black' else 'black'
    polygon = Polygon([(lon_breaks[-2], lat_breaks[0] - frame_height),
                       (lon_breaks[-1] + frame_width,
                        lat_breaks[0] - frame_height),
                       (lon_breaks[-1], lat_breaks[0]),
                       (lon_breaks[-2], lat_breaks[0])],
                       facecolor=colour, edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
    colour = 'white' if colour == 'black' else 'black'
    polygon = Polygon([(lon_breaks[-1], lat_breaks[0]),
                       (lon_breaks[-1] + frame_width,
                        lat_breaks[0] - frame_height),
                       (lon_breaks[-1] + frame_width, lat_breaks[1]),
                       (lon_breaks[-1], lat_breaks[1])],
                        facecolor=colour, edgecolor='black', clip_on=False,
                        ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
        
    # Right bar
    for i, (from_lat, to_lat) in enumerate(list(zip(lat_breaks[1:-2],
                                                    lat_breaks[2:-1]))):    
        colour = 'white' if colour == 'black' else 'black'
        rectangle = Rectangle((lon_breaks[-1], from_lat),
                              width=frame_width,
                              height=to_lat - from_lat, facecolor=colour,
                              edgecolor='black', clip_on=False, ls='solid',
                              lw=1, transform=custom_crs, zorder=4)
        ax.add_patch(rectangle)
        
    colour = 'white' if colour == 'black' else 'black'
    polygon = Polygon([(lon_breaks[-1], lat_breaks[-2]),
                       (lon_breaks[-1] + frame_width, lat_breaks[-2]),
                       (lon_breaks[-1] + frame_width,
                        lat_breaks[-1] + frame_height),
                       (lon_breaks[-1], lat_breaks[-1])],
                       facecolor=colour, edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
    colour = 'white' if colour == 'black' else 'black'
    polygon = Polygon([(lon_breaks[-2], lat_breaks[-1]),
                       (lon_breaks[-1], lat_breaks[-1]),
                       (lon_breaks[-1] + frame_width,
                        lat_breaks[-1] + frame_height),
                       (lon_breaks[-2], lat_breaks[-1] + frame_height)],
                       facecolor=colour, edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
        
    # Top bar, from right to left
    for i, (from_lon, to_lon) in enumerate(list(zip(np.flip(lon_breaks[1:-2], 0),
                                                    np.flip(lon_breaks[2:-1], 0)))):
        colour = 'white' if colour == 'black' else 'black'
        rectangle = Rectangle((from_lon, lat_breaks[-1]),
                              width=to_lon - from_lon,
                              height=frame_height, facecolor=colour,
                              edgecolor='black', clip_on=False, ls='solid',
                              lw=1, transform=custom_crs, zorder=4)
        ax.add_patch(rectangle)
        
    colour = 'white' if colour == 'black' else 'black'
    polygon = Polygon([(lon_breaks[0], lat_breaks[-1]),
                       (lon_breaks[1], lat_breaks[-1]),
                       (lon_breaks[1], lat_breaks[-1] + frame_height),
                       (lon_breaks[0] - frame_width,
                        lat_breaks[-1] + frame_height)],
                       facecolor=colour, edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
    colour = 'white' if colour == 'black' else 'black'
    polygon = Polygon([(lon_breaks[0] - frame_width, lat_breaks[-2]),
                       (lon_breaks[0], lat_breaks[-2]),
                       (lon_breaks[0], lat_breaks[-1]),
                       (lon_breaks[0] - frame_width,
                        lat_breaks[-1] + frame_height)],
                       facecolor=colour, edgecolor='black', clip_on=False,
                       ls='solid', lw=1, transform=custom_crs, zorder=4)
    ax.add_patch(polygon)
        
    # Left bar, from top to bottom
    for i, (from_lat, to_lat) in enumerate(list(zip(np.flip(lat_breaks[1:-2], 0),
                                                    np.flip(lat_breaks[2:-1], 0)))):
        colour = 'white' if colour == 'black' else 'black'
        rectangle = Rectangle((lon_breaks[0] - frame_width, from_lat),
                              width=frame_width,
                              height=to_lat - from_lat, facecolor=colour,
                              edgecolor='black', clip_on=False, ls='solid',
                              lw=1, transform=custom_crs, zorder=4)
        ax.add_patch(rectangle)
        
def draw_compass(ax, size=0.8):
    """Draw a compass on the terrain map
    
    Arguments:
        ax: matplotlib axis; axis object as returned by plt.subplots
        
        size: float; compass size factor to adjust its size. default is 0.8
    """ 
    # Get axis extent
    lon_min, lon_max = ax.get_xlim()
    lat_min, lat_max = ax.get_ylim()
    dlon = (lon_max - lon_min)
    dlat = (lat_max - lat_min)
    
    # Compass parameters
    compass_size = size*min(dlon, dlat)
    compass_center = np.array([lon_min+0.5*dlon, lat_min+0.91*dlat])
    compass_radius = (compass_size*0.75)/2.
    
    compass_N = np.array([0., compass_radius])
    compass_E = np.array([compass_radius, 0.])
    compass_S = np.array([0., -compass_radius])
    compass_W = np.array([-compass_radius, 0.])
    
    circle = Circle(compass_center, compass_radius, facecolor='w',
                    edgecolor='k', zorder=2)
    ax.add_patch(circle)
    
    # Main black needles
    colour = 'black'
    polygon = Polygon([compass_center,
                       compass_center + 0.2*(compass_N + compass_W),
                       compass_center + compass_N],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center,
                       compass_center + 0.2*(compass_E + compass_N),
                       compass_center + compass_E],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center,
                       compass_center + 0.2*(compass_S + compass_E),
                       compass_center + compass_S],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center,
                       compass_center + 0.2*(compass_W + compass_S),
                       compass_center + compass_W],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    
    # Secondary black needles
    polygon = Polygon([compass_center + 0.2*(compass_N + compass_W),
                       compass_center + 0.57*(compass_N + compass_W),
                       compass_center + 0.1*compass_N + 0.37*compass_W],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center + 0.2*(compass_N + compass_E),
                       compass_center + 0.57*(compass_N + compass_E),
                       compass_center + 0.1*compass_E + 0.37*compass_N],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center + 0.2*(compass_S + compass_E),
                       compass_center + 0.57*(compass_S + compass_E),
                       compass_center + 0.1*compass_S + 0.37*compass_E],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center + 0.2*(compass_S + compass_W),
                       compass_center + 0.57*(compass_S + compass_W),
                       compass_center + 0.1*compass_W + 0.37*compass_S],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    
    # Main white needles
    colour = 'white'
    polygon = Polygon([compass_center,
                       compass_center + compass_N,
                       compass_center + 0.2*(compass_N + compass_E)],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center,
                       compass_center + compass_E,
                       compass_center + 0.2*(compass_E + compass_S)],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center,
                       compass_center + compass_S,
                       compass_center + 0.2*(compass_S + compass_W)],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center,
                       compass_center + compass_W,
                       compass_center + 0.2*(compass_W + compass_N)],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=3)
    ax.add_patch(polygon)
    
    # Secondary white needles
    polygon = Polygon([compass_center + 0.2*(compass_N + compass_W),
                       compass_center + 0.57*(compass_N + compass_W),
                       compass_center + 0.1*compass_W + 0.37*compass_N],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center + 0.2*(compass_N + compass_E),
                       compass_center + 0.57*(compass_N + compass_E),
                       compass_center + 0.1*compass_N + 0.37*compass_E],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center + 0.2*(compass_S + compass_E),
                       compass_center + 0.57*(compass_S + compass_E),
                       compass_center + 0.1*compass_E + 0.37*compass_S],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    polygon = Polygon([compass_center + 0.2*(compass_S + compass_W),
                       compass_center + 0.57*(compass_S + compass_W),
                       compass_center + 0.1*compass_S + 0.37*compass_W],
                       facecolor=colour, edgecolor='k', clip_on=False,
                       ls='solid', lw=1, zorder=2)
    ax.add_patch(polygon)
    
    # Write directions
    ax.text(compass_center[0] + 1.05*compass_N[0],
            compass_center[1] + 1.05*compass_N[1],
            "N", verticalalignment='bottom', horizontalalignment='center',
             fontsize=8, weight="bold", zorder=4)
    ax.text(compass_center[0] + 1.05*compass_E[0],
            compass_center[1] + 1.05*compass_E[1],
            "E", verticalalignment='center', horizontalalignment='left',
             fontsize=8, weight="bold", zorder=4)
    ax.text(compass_center[0] + 1.05*compass_S[0],
            compass_center[1] + 1.05*compass_S[1],
            "S", verticalalignment='top', horizontalalignment='center',
             fontsize=8, weight="bold", zorder=4)
    ax.text(compass_center[0] + 1.05*compass_W[0],
            compass_center[1] + 1.05*compass_W[1],
            "W", verticalalignment='center', horizontalalignment='right',
             fontsize=8, weight="bold", zorder=4)
             
def draw_legend(ax, seeding_locations, probes_obj):
    """Draw a legend on the terrain map listing the different elements
    present (land, coastline, farms, probes), and the name of the application
    
    Arguments:
        ax: matplotlib axis; axis object as returned by plt.subplots
        
        seeding_locations: list of Farm objects as implemented in claws.Farm;
            Seeding locations
        
        probes_obj: list of Probe objects as implemented in claws.Probe
    """ 
    lon_min, lon_max = ax.get_xlim()
    lat_min, lat_max = ax.get_ylim()
    dlon = (lon_max - lon_min)
    dlat = (lat_max - lat_min)
    
    symb_lon = lon_min + dlon*0.2
    text_lon = lon_min + dlon*0.4
    
    # Legend frame
    rectangle = Rectangle((0., lat_min + 0.221*dlat),
                           width=dlon, height=(1.-0.221)*dlat,
                           facecolor='white', edgecolor='black', clip_on=False,
                           ls='solid', lw=1, zorder=1)
    ax.add_patch(rectangle)
    
    # Print application name
    ax.text(lon_min + dlon*0.5, lat_min + dlat*0.24, "CLAWS",
            verticalalignment='center', horizontalalignment='center',
            fontsize=15, weight="bold", zorder=4)
    
    # Print land, coastline, farms and probes legend
    text_spacing = 0.04*dlat
    text_pos = lat_min + dlat*0.5
    
    rectangle = Rectangle((lon_min + dlon*0.1, text_pos-0.01*dlat),
                           width=0.2*dlon, height=0.02*dlat,
                           facecolor='none', edgecolor='grey', clip_on=False,
                           ls='solid', lw=0.5, zorder=5)
    ax.add_patch(rectangle)
    ax.text(text_lon, text_pos, "Domain",
            verticalalignment='center', horizontalalignment='left',
            fontsize=10, zorder=4)
            
    text_pos -= text_spacing
    rectangle = Rectangle((lon_min + dlon*0.1, text_pos-0.01*dlat),
                           width=0.2*dlon, height=0.02*dlat,
                           facecolor=cfeature.COLORS['land'],
                           edgecolor=cfeature.COLORS['land'], clip_on=False,
                           ls='solid', lw=1, zorder=5)
    ax.add_patch(rectangle)
    ax.text(text_lon, text_pos, "Land",
            verticalalignment='center', horizontalalignment='left',
            fontsize=10, zorder=4)
            
    text_pos -= text_spacing 
    ax.plot([lon_min + dlon*0.1, lon_min + dlon*0.3],
             [text_pos, text_pos],
             linestyle = '-', lw=1, color='black', zorder=4)
    ax.text(text_lon, text_pos, "Shoreline",
            verticalalignment='center', horizontalalignment='left',
            fontsize=10, zorder=4)
            
    text_pos -= text_spacing 
    ax.plot(symb_lon, text_pos, marker='o', color='black',
            markersize=4, alpha=0.7)
    ax.text(text_lon, text_pos, "Farm ({})".format(len(seeding_locations)),
            verticalalignment='center', horizontalalignment='left',
            fontsize=10, zorder=4)
    
    text_pos -= text_spacing 
    ax.plot(symb_lon, text_pos, marker='+', color='black',
            markersize=4, alpha=0.7)        
    ax.text(text_lon, text_pos, "Probe ({})".format(len(probes_obj)),
            verticalalignment='center', horizontalalignment='left',
            fontsize=10, zorder=4)
        
def draw_scale(ax, corners, length_pct=0.125):
    """Draw a thick white and black scale on the terrain map
    
    Arguments:
        ax: matplotlib axis; axis object as returned by plt.subplots
        
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the terrain map
        
        length_pct: float; scale length expressed as a percentage of the
            terrain map width. default is 12.5%
    """ 
    lon_min, lon_max = ax.get_xlim()
    lat_min, lat_max = ax.get_ylim()
    dlon = lon_max - lon_min
    dlat = lat_max - lat_min
    lon_mid = lon_min + dlon/2.
    
    aspect_ratio = 0.07
    # List of distances that can be printed on the scale
    accepted_distances = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1., 2.5, 5., 10.,
                          25., 50., 100., 250.]
    
    r = abs(lon_max - lon_min)/abs(corners[1] - corners[0])
    
    # Target scale length in meters
    target_scale_length_m = lonrange_to_distance(corners)*length_pct
    # Convert distance to output distance units
    target_scale_length_out = convert_len(target_scale_length_m,
                                          output_options["length_units"],
                                          "claws.postpro.draw_scale")
    # Find closest distance in the list of acceptable scale distances
    scale_length_out = find_nearest(accepted_distances, target_scale_length_out)
    # Convert this distance back to a longitude range
    scale_length_m = convert_len_to_prog_units(scale_length_out,
                                               output_options["length_units"],
                                               "claws.postpro.draw_scale")
    scale_length_lon = distance_to_lonrange(scale_length_m, corners[2],
                                            corners[3])
    
    # Scale parameters
    scale_length = scale_length_lon*r
    scale_width = aspect_ratio*scale_length
    scale_lonmid = lon_min + 0.975*dlon - scale_length/2.
    
    scale_blhc = np.array([scale_lonmid - scale_length/2.,
                           lat_min + 0.035*dlat])
    lon_breaks = np.linspace(scale_blhc[0], scale_blhc[0] + scale_length, 6)
    
    # Draw scale
    colour = 'white'
    for i, (from_lon, to_lon) in enumerate(list(zip(lon_breaks[:-1],
                                                    lon_breaks[1:]))):
        colour = 'white' if colour == 'black' else 'black'
        rectangle = Rectangle((from_lon, scale_blhc[1]),
                              width=to_lon - from_lon,
                              height=scale_width, facecolor=colour,
                              edgecolor='black', clip_on=False, ls='solid',
                              lw=1, zorder=5)
        ax.add_patch(rectangle)
        
    # Add distance label
    scale_length_out_str = ('%f'%scale_length_out).rstrip('0').rstrip('.')
    ax.text(scale_lonmid, scale_blhc[1] + scale_width + 0.01*dlat,
            '{} {}'.format(scale_length_out_str,
                           output_options["length_units"]),
            verticalalignment='bottom', horizontalalignment='center',
            fontsize=8, zorder=4)
            
def draw_minimap(fig, gs, seeding_locations, label):   
    """Draw a minimap representing the whole of Scotland to help locate the farm
    
    Arguments:
        fig: matplotlib figure object
        
        gs: gridspec location; where to print the minimap on the terrain map
            canvas
            
        seeding_locations: list of Farm objects as implemented in
            claws.Farm; Seeding locations
        
        label: float; scale length expressed as a percentage of the
            terrain map width. default is 12.5%
    """
    # Create another axis to draw the minimap onto. The axis is located at the
    # same position than the right panel, position given by gs
    ax = fig.add_subplot(gs, projection=ccrs.PlateCarree())
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    # Set minimap extent to be that of Scotland
    scotland_extent = (-7.5, -1.5, 53.5, 59.5)
    lon_mid = (scotland_extent[0] + scotland_extent[1])/2.
    ax.set_extent(scotland_extent, ccrs.PlateCarree())
    
    # Add the land and coastlines
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    
    # Mark the farm location
    farm_lon = seeding_locations[0].Site().lon()
    farm_lat = seeding_locations[0].Site().lat()
    farm_loc = (farm_lon, farm_lat)
    circle = Circle(farm_loc, 0.2, facecolor='salmon', edgecolor='darkred',
                    zorder=2)
    ax.add_patch(circle)
    
    # Add text about the location
    ax.text(lon_mid, scotland_extent[2] - 0.75, label,
            verticalalignment='center', horizontalalignment='center',
            fontsize=10, fontstyle ='italic', zorder=4)
            
def plot_series_peak_concentration(working_folder, time, peakconc,
                                   quadtree_peakconc, chemical_obj,
                                   quadtree_obj, pixelsize_m, dz,
                                   time_last_treatment, last_treatment_label):
    """Plot time series of peak concentration
    
    Arguments:
        working_folder: string; working folder
        
        time: numpy array; time series
        
        peakconc: numpy array; peak concentration series using the root bins
        
        quadtree_peakconc: numpy array; peak concentration series using the leaf
            bins
        
        chemical_obj: ChemicalSubstance object as implemented in
            claws.ChemicalSubstance
            
        quadtree_obj: Quadtree object as implemented in
            claws.quadtree
        
        pixelsize_m: float; size of the bins along the longitudinal and
            latitudinal directions, in meters 
            
        dz: float; depth bin size in meters       
        
        time_last_treatment: float; time of latest treatment over all
            seeding locations (units: output_options["time_units"]) 
        
        last_treatment_label: string; message to print on the graph at time
            of last treatment. Default is 'Last treatment'  
    """     
    # Conversion factors from program units to user-defined output units
    cf, lf, tf = get_unit_factors()
    
    fig, ax = plt.subplots(1)
    species_name = chemical_obj.name()
    psize_str = ('%f'%pixelsize_m).rstrip('0').rstrip('.')
    dz_str = ('%f'%dz).rstrip('0').rstrip('.')
    plt.plot(time, peakconc, color='black', linestyle='-', lw=1,
             label='{0} m x {0} m x {1} m'.format(psize_str, dz_str))
    if quadtree_obj.is_active():
        plt.scatter(output_options["time_bins"]*tf,
                    quadtree_peakconc, color='blue', marker='+',
                    label='quadtree (max depth: {})'.format(
                        quadtree_obj.max_quadtree_depth()))         
    if np.isfinite(chemical_obj.eqs_3hour()):
        plt.axhline(y=chemical_obj.eqs_3hour()*cf, color='green',
                    linestyle='dotted', label='3-hour EQS')
    if np.isfinite(chemical_obj.eqs_72hour()):
        plt.axhline(y=chemical_obj.eqs_72hour()*cf, color='green',
                    linestyle='--', label='72-hour EQS')
    if np.isfinite(chemical_obj.mac_72hour()):
        plt.axhline(y=chemical_obj.mac_72hour()*cf, color='red',
                    linestyle='--', label='72-hour MAC')
    ax.set_yscale('log')
    round_logscale_yaxis(ax)
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    print_last_treatment(plt, ax, time_last_treatment,
                         label=last_treatment_label)
    plt.xlabel('Time ({})'.format(output_options["time_units"]))
    plt.ylabel('Peak {} concentration ({})'.format(species_name,
        output_options["concentration_units"]))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=2)
    plt.savefig(working_folder + 'peak_concentration.png')
    plt.close()
    
def plot_series_area_greater_EQS(working_folder, time, area_over_eqs,
                                 quadtree_area_over_eqs, chemical_obj,
                                 quadtree_obj, pixelsize_m, dz,
                                 time_last_treatment, last_treatment_label):
    """Plot time series of area greater than EQS
    
    Arguments:
        working_folder: string; working folder
        
        time: numpy array; time series
        
        area_over_eqs: numpy array; area greater than EQS series using the root
            bins
        
        quadtree_area_over_eqs: numpy array; area greater than EQS series using
            the leaf bins
        
        chemical_obj: ChemicalSubstance object as implemented in
            claws.ChemicalSubstance
            
        quadtree_obj: Quadtree object as implemented in
            claws.quadtree
        
        pixelsize_m: float; size of the bins along the longitudinal and
            latitudinal directions, in meters 
            
        dz: float; depth bin size in meters       
        
        time_last_treatment: float; time of latest treatment over all
            seeding locations (units: output_options["time_units"]) 
        
        last_treatment_label: string; message to print on the graph at time
            of last treatment. Default is 'Last treatment'  
    """    
    # Conversion factors from program units to user-defined output units
    cf, lf, tf = get_unit_factors()
        
    fig, ax = plt.subplots(1)
    species_name = chemical_obj.name()
    eqs_aze = chemical_obj.eqs_aze()*lf**2
    psize_str = ('%f'%pixelsize_m).rstrip('0').rstrip('.')
    dz_str = ('%f'%dz).rstrip('0').rstrip('.')
    plt.plot(time, area_over_eqs, color='black', linestyle='-', lw=1,
             label='{0} m x {0} m x {1} m'.format(psize_str, dz_str))
    if quadtree_obj.is_active():
        plt.scatter(output_options["time_bins"]*tf,
                    quadtree_area_over_eqs,
                    color='blue', marker='+',
                    label='quadtree (max depth: {})'.format(
                        quadtree_obj.max_quadtree_depth()))  
    plt.axhline(y=eqs_aze, color='g', linestyle='--', label='72-hour EQS AZE')
    ax.set_yscale('log')    
    round_logscale_yaxis(ax)  
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    print_last_treatment(plt, ax, time_last_treatment,
                         label=last_treatment_label)
    plt.xlabel('Time ({})'.format(output_options["time_units"]))
    plt.ylabel('Area > EQS ({}$^2$)'.format(output_options["length_units"]))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=2)
    plt.savefig(working_folder + 'area_over_EQS.png')
    plt.close()
    
def plot_series_area_greater_MAC(working_folder, time, area_over_mac,
                                 quadtree_area_over_mac, chemical_obj,
                                 quadtree_obj, pixelsize_m, dz,
                                 time_last_treatment, last_treatment_label):
    """Plot time series of area greater than MAC
    
    Arguments:
        working_folder: string; working folder
        
        time: numpy array; time series
        
        area_over_mac: numpy array; area greater than MAC series using the root
            bins
        
        quadtree_area_over_mac: numpy array; area greater than MAC series using
            the leaf bins
        
        chemical_obj: ChemicalSubstance object as implemented in
            claws.ChemicalSubstance
            
        quadtree_obj: Quadtree object as implemented in
            claws.quadtree
        
        pixelsize_m: float; size of the bins along the longitudinal and
            latitudinal directions, in meters 
            
        dz: float; depth bin size in meters       
        
        time_last_treatment: float; time of latest treatment over all
            seeding locations (units: output_options["time_units"]) 
        
        last_treatment_label: string; message to print on the graph at time
            of last treatment. Default is 'Last treatment'  
    """    
    # Conversion factors from program units to user-defined output units
    cf, lf, tf = get_unit_factors()
        
    fig, ax = plt.subplots(1)
    species_name = chemical_obj.name()
    mac_aze = chemical_obj.mac_aze()*lf**2
    psize_str = ('%f'%pixelsize_m).rstrip('0').rstrip('.')
    dz_str = ('%f'%dz).rstrip('0').rstrip('.')
    plt.plot(time, area_over_mac, color='black', linestyle='-', lw=1,
             label='{0} m x {0} m x {1} m'.format(psize_str, dz_str))
    if quadtree_obj.is_active():
        plt.scatter(output_options["time_bins"]*tf,
                    quadtree_area_over_mac,
                    color='blue', marker='+',
                    label='quadtree (max depth: {})'.format(
                        quadtree_obj.max_quadtree_depth()))  
    plt.axhline(y=mac_aze, color='r', linestyle='--', label='72-hour MAC AZE')
    ax.set_yscale('log')    
    round_logscale_yaxis(ax)  
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    print_last_treatment(plt, ax, time_last_treatment,
                         label=last_treatment_label)
    plt.xlabel('Time ({})'.format(output_options["time_units"]))
    plt.ylabel('Area > MAC ({}$^2$)'.format(output_options["length_units"]))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=2)
    plt.savefig(working_folder + 'area_over_MAC.png')
    plt.close()
    
def plot_vertical_distribution_bar_graph(working_folder, vdist, bar_height,
                                         nmax):
    """Plot vertical distribution bar graph for all output times
    
    Arguments:
        working_folder: string; working folder
        
        vdist: numpy array; particle vertical distribution series
        
        bar_height: float; depth interval in meters
        
        nmax: integer; maximum number of particle per layer at any time
    """
    for tbin in output_options["time_bins"]:
        arr = np.array([100.*v[tbin].item() for v in vdist])
        len_arr = len(arr)
        zbin = np.linspace(start=0, stop=-bar_height*(len_arr-1), num=len_arr)
        
        fig, ax = plt.subplots(1)
        plt.barh(zbin, arr, height=-bar_height, facecolor='navy', alpha=0.35,
                 lw=1, edgecolor='k', align='edge')
        nyticks, mltp = compute_nticks(len_arr)
        plt.yticks(np.linspace(start=0, stop=-bar_height*mltp*(nyticks-1),
                               num=nyticks))
        plt.xlabel('Number of particles (%)')
        plt.ylabel('Depth (m)')
        plt.xlim([0, min(100, round(100.*nmax,0)//5*5.+5.)])
        plt.title(np.datetime64(vdist[0].time[tbin].item(), 'ns').astype(
                  'datetime64[us]').astype(datetime.datetime).strftime(
                  '%d %B %Y, %I:%M:%S %p')) 
        plt.savefig(working_folder + 'vertical_distribution'\
                        '_{:04d}.png'.format(tbin))
        plt.close()

def plot_series_nparticles_in_polygon(working_folder, polygon, time, lons, lats,
                                      statuses):
    """Plot time series of number of particles in a control area (polygon), eg,
    for flushing time calculation, and return the time series of the number of
    particles in the control area, and the flushing time
    
    Arguments:
        working_folder: string; working folder
        
        polygon: GeoJSON polygon; a control area
        
        time: numpy array; time series
        
        lons: numpy array; longitude history for all times and particles
        
        lats: numpy array; latitude history for all times and particles
        
        statuses: numpy array; status history for all times and particles
    """     
    flushing_time = np.nan
    ndt = len(time)
    
    # Conversion factors from program units to user-defined output units
    cf, lf, tf = get_unit_factors()
    
    nparticles_in_polygon = get_nparticles_in_polygon(polygon, lons, lats,
                                                      statuses)
    
    # Normalise
    nparticles_in_polygon *= 100./nparticles_in_polygon[0]
    
    npart_at_ft = 37.0
    npart_at_ft_str = ('%f'%npart_at_ft).rstrip('0').rstrip('.')
    
    # Derive flushing time
    if np.min(nparticles_in_polygon) > npart_at_ft:
        has_flushed = False
    else:
        has_flushed = True
    
    if has_flushed:
        # Find latest time satisfying npart > 37%
        ft_idces = np.where(nparticles_in_polygon < 37.)[0]
        for i in reversed(range(len(ft_idces))):
            if ft_idces[i] != ft_idces[i-1] + 1:
                ft_idx = ft_idces[i]
                break
            
        # Linear interpolation of the flushing time
        ft_idx_m = ft_idx - 1
        ft_idx_p = ft_idx
        val_m = nparticles_in_polygon[ft_idx_m]
        val_p = nparticles_in_polygon[ft_idx_p]
        slope = (val_p - val_m)/(time[ft_idx_p] - time[ft_idx_m])
        flushing_time = time[ft_idx_m] + (npart_at_ft - val_m)/slope
        
    # Plot time series and print the flushing time
    fig, ax = plt.subplots(1)
    plt.plot(time, nparticles_in_polygon, color='black', linestyle='-', lw=1)
    plt.axhline(y=npart_at_ft, color='black', linestyle='dotted')
    if has_flushed:
        plt.axvline(x=flushing_time, color='black', linestyle='dotted')
    xmin, xmax = ax.get_xlim()
    ax.text(xmin + 0.02*(xmax - xmin), npart_at_ft,
            "{}%".format(npart_at_ft_str), va='bottom')
    if has_flushed:
        ymin, ymax = ax.get_ylim()
        time_units_str = output_options["time_units"]
        if flushing_time > 1. and time_units_str == 'day':
            time_units_str += 's'
        ax.text(flushing_time, ymax - 0.02*(ymax - ymin), 
                "Flushing time ~ {:.2f} {}".format(flushing_time,
                time_units_str), rotation=90, ha='right', va='top')
                
    ax.set_yscale('log')
    plt.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter('')
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(ticker.ScalarFormatter())
    plt.xlabel('Time ({})'.format(output_options["time_units"]))
    plt.ylabel('Number of particles (%)')
    plt.savefig(working_folder + 'flushing_time.png')
    plt.close()
    
    flushing_time = convert_time_to_prog_units(
        flushing_time,
        output_options["time_units"],
        "claws.postpro.plot_series_nparticles_in_polygon")
    
    return nparticles_in_polygon, flushing_time

def round_logscale_axis(arr, leave_top=False, leave_bottom=False,
                        set_min=-np.inf, set_max=np.inf):
    """Round the y-axis min and max values to the nearest power of ten
    
    Arguments:
        arr: list; python list or numpy array containing the min and max axis
            values
        
        leave_top: boolean; whether to not round the axis maximum value. default
            is False
            
        leave_bottom: boolean; whether to not round the axis minimum value.
            default is False
        
        set_min: float; minimum value to impose, needs to be greater than the 
            mnimum axis value to take effect. default is -inf    
        
        set_max: float; maximum value to impose, needs to be lesser than the 
            mnimum axis value to take effect. default is +inf       
    """   
    ax_min = max(arr[0], set_min)
    ax_max = min(arr[1], set_max)
    
    int_max = np.ceil(np.log10(ax_max))
    if ax_min == 0.0:
        int_min = int_max - 4
    else:
        int_min = np.floor(np.log10(ax_min))
        
    if not leave_top and not leave_bottom:
        ax_min = np.power(10, int_min)
        ax_max = np.power(10, int_max)
    elif leave_top:
        ax_min = np.power(10, int_min)
    elif leave_bottom:
        ax_max = np.power(10, int_max)
    return [ax_min, ax_max]
    
def round_logscale_yaxis(ax, leave_top=False, leave_bottom=False,
                         set_min=-np.inf, set_max=np.inf):
    """Round the y-axis min and max values to the nearest power of ten
    
    Arguments:
        ax: matplotlib axis; axis object as returned by plt.subplots
        
        leave_top: boolean; whether to not round the axis maximum value. default
            is False
            
        leave_bottom: boolean; whether to not round the axis minimum value.
            default is False
            
        set_min: float; minimum value to impose, needs to be greater than the 
            mnimum axis value to take effect. default is -inf    
        
        set_max: float; maximum value to impose, needs to be lesser than the 
            mnimum axis value to take effect. default is +inf   
    """   
    rounded_bounds = round_logscale_axis(ax.get_ylim(), leave_top, leave_bottom,
                                         set_min, set_max)
    ax.set_ylim(rounded_bounds)
    return ax
    
def compute_nticks(nbars):
    """Compute the number of ticks from the number of bars on a histogram and
    attempts to reduce their number when there are too many
    
    Arguments:
        nbars: integer; number of intervals
    """
    nticks = nbars + 1
    multiplicator = 1
    
    if nticks > 10:
        if nticks%2 == 0:
            nticks = nticks//2
            multiplicator = 2
        elif nticks%3 == 0:
            nticks = nticks//3
            multiplicator = 3
        elif nticks%5 == 0:
            nticks = nticks//5
            multiplicator = 5
    return [nticks, multiplicator]

def print_last_treatment(plt, ax, t_last_treatment, label=None):
    """Round the y-axis min and max values to the nearest power of ten
    
    Arguments:
        plt: matplotlib plot object
        
        ax: matplotlib axis object
        
        t_last_treatment: float; time of the last treatment. Units are that
            defined by the user in output_options, time_units
            
        label: string; message to print on plots. default is None, which prints
            'Last treatment'  
    """ 
    if label is None:
        label = 'Last treatment'
    ax_min, ax_max = ax.get_ylim()
    plt.axvline(x=t_last_treatment, color='black', linestyle='--', lw=0.5)
    ax.text(1.03*t_last_treatment, 1.05*ax_min, label, rotation=90, va='bottom')
    return [plt, ax]
    
def animate_concentration_terrain(working_folder):
    """Animate concentration on terrain map to produce of mp4 video from png
    images
    
    Arguments:
        wf: string; working folder
    """
    if not animation["animate"]: return
    # A video cannot be created from a single frame, vlc would crash
    if len(output_options["time_bins"]) == 1: return
    
    # Shorthands
    imrate = animation["image_rate"]
    frate = animation["frame_rate"]
    ov = 'concentration'
    os.system("ffmpeg -y -framerate {} -pattern_type glob -i " \
              "'{}concentration_*.png' -c:v libx264 -r {} " \
              "-pix_fmt yuv420p {}.mp4 > /dev/null 2>&1".format(imrate,
              working_folder, frate, working_folder + ov))
                                                                
def animate_vertical_distribution(working_folder):
    """Animate particle vertical distribution to produce of mp4 video from png
    images
    
    Arguments:
        working_folder: string; working folder
    """
    if not animation["animate"]: return
    # A video cannot be created from a single frame, vlc would crash
    if len(output_options["time_bins"]) == 1: return
    
    # Shorthands
    imrate = animation["image_rate"]
    frate = animation["frame_rate"]
    ov = 'vertical_distribution'
    os.system("ffmpeg -y -framerate {} -pattern_type glob -i " \
              "'{}vertical_distribution_*.png' -c:v libx264 -r {} " \
              "-pix_fmt yuv420p {}.mp4 > /dev/null 2>&1".format(
              imrate, working_folder, frate, working_folder + ov))
                                                                
def print_time_series_to_file(filename, input_list):
    """Print time series provided in the form of a list of tuples to a .dat
    file. The first element of the input list is the name of the series and the
    second element of the tuple is the list
    
    Arguments:
        filename: string; output file name (absolute path)
        
        input_list: list of tuples; format as given in the description
    """

    if not output_options["write_time_series_to_file"]: return
    
    np.savetxt(fname=filename,
               X=np.array([i[1] for i in input_list]).transpose(),
               header='\t'.join(str(i[0]) for i in input_list))
