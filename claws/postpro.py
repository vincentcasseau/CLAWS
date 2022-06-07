#!/usr/bin/env python
"""Post-processing functions.
"""

# Import modules
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from matplotlib.ticker import FormatStrFormatter

from opendrift.readers import reader_global_landmask
import cartopy.crs as ccrs
from cartopy.crs import CRS, Globe
import cartopy.io.img_tiles as cimgt

from claws.claws import output_options,\
    animation, get_unit_factors

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "hystrath@gmail.com"
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
    
def create_terrain(corners, seeding_locations, tile_style='',
                   display_locations=False, abbreviate_locations=True,
                   display_probe_markers=False, probes=None,
                   custom_crs=ccrs.Geodetic(),
                   globe=None):
    """Create a terrain map using a Stamen tile of size domain_extent. Seeding
    sources are labelled on the map.
    
    Arguments:
        corners: list of 4 floats: [lon_min, lon_max, lat_min, lat_max];   
            impose the extent of the histogram
            
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
    # Add small margins (2%) around the domain's corners
    lon_min, lon_max, lat_min, lat_max = corners
    lon_offset = (lon_max - lon_min)/50.
    lat_offset = (lat_max - lat_min)/50.
    domain_w_margins = [lon_min - lon_offset, lon_max + lon_offset,
                        lat_min - lat_offset, lat_max + lat_offset]
                        
    # Create figure
    fig = plt.figure(figsize=(8., 8.))
            
    # Create a Stamen instance
    available_tile_styles = ['terrain', 'terrain-background', 'watercolor',
                             'toner-lite', 'toner-background']
    if tile_style in available_tile_styles:
        tiler = cimgt.Stamen(tile_style)
        # Create a GeoAxes in the tile's projection
        ax = fig.add_subplot(1, 1, 1, projection=tiler.crs)
        globe = tiler.crs.globe
    else:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
        globe = ccrs.Mercator().globe
    ax.set_extent(domain_w_margins, crs=custom_crs)

    # Draw gridlines
    gl = ax.gridlines(ccrs.PlateCarree(globe=globe), draw_labels=True,
                      linewidth=0.5, color='grey', alpha=0.75, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    
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
    plt.plot([lon_min, lon_min], [lat_min, lat_max], color='grey', lw=0.5,
             alpha=0.5, transform=custom_crs)
    plt.plot([lon_max, lon_max], [lat_min, lat_max], color='grey', lw=0.5,
             alpha=0.5, transform=custom_crs)
    plt.plot([lon_min, lon_max], [lat_min, lat_min], color='grey', lw=0.5,
             alpha=0.5, transform=custom_crs)
    plt.plot([lon_min, lon_max], [lat_max, lat_max], color='grey', lw=0.5,
             alpha=0.5, transform=custom_crs)
             
#    reader_global_landmask.plot_land(ax, domain_w_margins[0],
#                               domain_w_margins[2],
#                               domain_w_margins[1],
#                               domain_w_margins[3], False,
#                               lscale='full', globe=ccrs.Mercator().globe)

    # Use the cartopy interface to create a matplotlib transform object
    # for the Geodetic coordinate system. 
    geodetic_transform = custom_crs._as_mpl_transform(ax)
    
    if display_probe_markers:
        if probes is not None:
            for P in probes:    
                # Print probe marker: '+'
                plt.plot(P.Station().lon(), P.Station().lat(), marker='+',
                         color='black', markersize=4, alpha=0.7,
                         transform=custom_crs)
    
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
        plt.plot(lon, lat, marker='o', color='black', markersize=4, alpha=0.7,
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
                
            text_transform = offset_copy(geodetic_transform, units='dots', x=nx,
                                         y=ny)
            plt.text(lon, lat, loc_text,
                     verticalalignment=va, horizontalalignment=ha,
                     transform=text_transform, fontsize=9,
                     bbox=dict(facecolor='sandybrown', alpha=0.5,
                               boxstyle='round'))

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
    
def animate_concentration_terrain(wf):
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
              "-pix_fmt yuv420p {}.mp4 > /dev/null 2>&1".format(imrate, wf,
                                                                frate, wf + ov))
                                                                
def animate_vertical_distribution(wf):
    """Animate particle vertical distribution to produce of mp4 video from png
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
    ov = 'vertical_distribution'
    os.system("ffmpeg -y -framerate {} -pattern_type glob -i " \
              "'{}vertical_distribution_*.png' -c:v libx264 -r {} " \
              "-pix_fmt yuv420p {}.mp4 > /dev/null 2>&1".format(imrate, wf,
                                                                frate, wf + ov))
                                                                
def print_time_series_to_file(wf, input_list):
    """Print time series provided in the form of a list of tuples to a .dat
    file. The first element is the name of the file to print, minus its path and
    extension, and the second element of the tuple is the list
    
    Arguments:
        wf: string; working folder
        
        input_list: list of tuples; format as given in the description
    """

    if not output_options["write_time_series_to_file"]: return
    
    for i in input_list:
        filename, list2print = i
        if len(np.shape(list2print)) == 2:
            list2print = list2print.transpose()
        np.savetxt(fname=wf + filename + '.dat', X=list2print)