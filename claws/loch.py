#!/usr/bin/env python
"""Root class defining a Loch
"""

# Import modules
from os.path import exists as file_exists
import re
import json
import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

from claws.custom_exceptions import InputError
from claws.helpers import *
from claws.claws import output_options
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
# Functions 
# ---------------------------------------------------------------------------- #

def tidal_range(working_folder, sea_water_elevation_file, input_len_units='m'):
    """Compute the min/mean/max tidal ranges over a time period

    Tidal range: vertical difference in height between consecutive high and low
    waters over a tidal cycle
    
    Arguments:
        working_folder: string; working folder
        
        sea_water_elevation_file: string; full path to the sea water elevation
            file. It is a single-column file with no header line
            
        input_len_units: string; input elevation units. default is meters
    """    
    # Check if the input sea water elevation file exists. If not, do nothing
    if not(sea_water_elevation_file and file_exists(sea_water_elevation_file)):
        print("Sea water elevation file not provided")
        return [None, None, None]
    
    # Read the input file
    elevation = np.genfromtxt(sea_water_elevation_file)
    # Convert elevation to program units
#    elevation = convert_len_to_prog_units(elevation, input_len_units,
#                                          'claws.loch.tidal_range')
    npoints = len(elevation)

    # Find positions of all local maxima
    flow_tide_elevation = argrelextrema(elevation, np.greater)[0]

    # Find positions of all local minima
    ebb_tide_elevation = argrelextrema(elevation, np.less)[0]

    # Count number of flow and ebb tides, and the number of tides as the minimum
    # of the two
    n_flow_tides = len(flow_tide_elevation)
    n_ebb_tides = len(ebb_tide_elevation)

    n_tides = min(n_flow_tides, n_ebb_tides)
    flow_tide_elevation = flow_tide_elevation[:n_tides]
    ebb_tide_elevation = ebb_tide_elevation[:n_tides]

    # Over the time period that is that given in the input file, if the ebb tide
    # comes first, the tidal range will be defined as the difference between the
    # ebb tide elevation and the consecutive high tide. If a flow tide comes
    # first, it will be the difference between the flow tide elevation and the
    # elevation of the subsequent ebb tide.
    tidal_range = \
        elevation[flow_tide_elevation] - elevation[ebb_tide_elevation]
        
    min_tidal_range = np.min(tidal_range)
    mean_tidal_range = np.mean(tidal_range)
    max_tidal_range = np.max(tidal_range)
    
    # Plot time series of the sea water elevation
    fig, ax = plt.subplots(1)
    plt.plot(np.arange(0, npoints), elevation, color='black', linestyle='-',
             lw=1, label='Elevation')
    plt.plot(flow_tide_elevation, elevation[flow_tide_elevation], color='red',
             marker='+', label='Flow tide')
    plt.plot(ebb_tide_elevation, elevation[ebb_tide_elevation], color='blue',
             marker='x', label='Ebb tide')
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=3)         
    plt.xlabel('Timestamp')
    plt.ylabel('Sea water elevation (m)')
    plt.savefig(working_folder + 'sea_water_elevation.png')
    plt.close()

    # Plot time series of the tidal range
    fig, ax = plt.subplots(1)
    plt.plot(np.arange(0, n_tides), tidal_range, color='black', linestyle='-',
             lw=1)
    plt.axhline(y=min_tidal_range, color='blue', linestyle='dashed',
                label='min. = {:.3f}'.format(min_tidal_range))
    plt.axhline(y=mean_tidal_range, color='green', linestyle='dashed',
                label='mean = {:.3f}'.format(mean_tidal_range))
    plt.axhline(y=max_tidal_range, color='red', linestyle='dashed',
                label='max. = {:.3f}'.format(max_tidal_range))
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand",
               borderaxespad=0, ncol=3)               
    plt.xlabel('Number of flow/ebb tides')
    plt.ylabel('Tidal range (m)')
    plt.savefig(working_folder + 'tidal_range.png')
    plt.close()
    return [min_tidal_range, mean_tidal_range, max_tidal_range]


# ---------------------------------------------------------------------------- #
# Classes 
# ---------------------------------------------------------------------------- #

class Loch(object):
    def __init__(self, area=np.inf, tidal_range=np.inf, volume=np.inf,
                 mean_depth=np.inf, geojson_file="", existing_biomass=0.0,
                 input_len_units='m', input_area_units='km^2',
                 input_vol_units='M m^3', input_mass_units='tonne',
                 reference="", name=""):
        # __name: string; Loch name. default is class name
        self.__name = name
        # __area: float; Loch's low water (LW) area (m^2). default is +inf
        self.__LWarea = area
        # __tidal_range: float; Loch's tidal range (m). default is +inf
        self.__tidal_range = tidal_range
        # __LWvolume: float; Loch's low water (LW) volume (m^3).
        # default is +inf
        self.__LWvolume = volume
        # __LWmean_depth: float; Loch's mean depth (m). default is +inf
        self.__LWmean_depth = mean_depth
        # __reference: string; Reference. default is empty string
        self.__reference = reference
        
        # __geojson_file: string; File storing a GeoJSON polygon to define the 
        # geographical extent of the loch. default is empty string
        self.__geojson_file = geojson_file
        self.__geojson_polygon = None
        # __existing_biomass: float; Existing nitrogen biomass in kg.
        # default is 0
        self.__existing_biomass = existing_biomass
        
        self._sanitize(input_len_units, input_area_units, input_vol_units,
                       input_mass_units)
        # Compute the flushing time (h) - it can be reset by the user
        self._compute_flushing_time()
        
    def __str__(self):
        return """{}\n\tLW Area (km^2) = {}\n\tLW Volume (M m^3) = {}\
            \n\tLW Mean depth (m) = {}\n\tTidal range (m) = {}\
            \n\tFlushing time (day) = {}\
            \n\tExisting nitrogen biomass (tonne) = {}\
            \n\tReference = {}""".format( self.name(), self.area('km^2'),
            self.volume('M m^3'), self.__LWmean_depth, self.__tidal_range,
            self.flushing_time('day'), self.existing_biomass('tonne'),
            self.__reference)
    
    def name(self):
        return self.__name
        
    def area(self, units='m^2'):
        return convert_area(self.__LWarea, units, self.name())
        
    def tidal_range(self, units='m'):
        return convert_len(self.__tidal_range, units, self.name())
        
    def volume(self, units='m^3'):
        return convert_vol(self.__LWvolume, units, self.name())
        
    def mean_depth(self, units='m'):
        return convert_len(self.__LWmean_depth, units, self.name())
        
    def geojson_polygon(self):
        return self.__geojson_polygon
        
    def geojson_polygon_geometry(self):
        if self.__geojson_polygon is not None:
            return self.__geojson_polygon['geometry']
        return None
        
    def geojson_polygon_geometry_coords(self):
        if self.__geojson_polygon is not None:
            return self.__geojson_polygon['geometry']['coordinates'][0]
        return None
        
    def averaged_coordinates(self):
        if self.__geojson_polygon is not None:
            coords = self.geojson_polygon_geometry_coords()
            lons = [pt[0] for pt in coords]
            lats = [pt[1] for pt in coords]
            return Site(name=self.name(), longitude=np.mean(lons),
                        latitude=np.mean(lats),
                        reference="{}'s GeoJSON polygon".format(self.name()))
        return None
        
    def geojson_polygon_properties(self):
        if self.__geojson_polygon is not None:
            return self.__geojson_polygon['properties']
        return None
    
    def existing_biomass(self, units='tonne'):
        return convert_mass(self.__existing_biomass, units, self.name())
        
    def set_geojson_polygon_properties(self, input_dict):
        if self.__geojson_polygon is not None:
            for key, value in input_dict.items():
                self.__geojson_polygon['properties'][key] = value
        
    def _sanitize(self, input_len_units, input_area_units, input_vol_units,
                  input_mass_units):
        if not self.__name:
            self.__name = ' '.join(re.findall('([A-Z][a-z]+)',
                                   type(self).__name__))
        
        if np.isfinite(self.__LWarea):
            assert(type(self.__LWarea) is float)
            if self.__LWarea <= 0.0:
                raise InputError(self.__LWarea,
                    "Loch's area should be a positive number")
                    
        assert(type(self.__tidal_range) is float)
        if self.__tidal_range <= 0.0:
            raise InputError(self.__tidal_range,
                "Loch's tidal range should be a positive number")
                    
        if np.isinf(self.__LWvolume) and np.isinf(self.__LWmean_depth):
            pass
        elif np.isinf(self.__LWvolume):
            assert(type(self.__LWmean_depth) is float)
            if self.__LWmean_depth < 0.0:
                self.__LWmean_depth *= -1.0
            elif self.__LWmean_depth == 0.0:
                raise InputError(self.__LWmean_depth,
                    "Loch's mean depth should be non-zero")
        elif np.isinf(self.__LWmean_depth):
            assert(type(self.__LWvolume) is float)
            if self.__LWvolume <= 0.0:
                raise InputError(self.__LWvolume,
                    "Loch's volume should be a positive number")
        else:
            assert(type(self.__LWvolume) is float)
            if self.__LWvolume <= 0.0:
                raise InputError(self.__LWvolume,
                    "Loch's volume should be a positive number")
            assert(type(self.__LWmean_depth) is float)
            if self.__LWmean_depth <= 0.0:
                raise InputError(self.__LWmean_depth,
                    "Loch's mean depth should be a positive number")
        
        # Area stored in m^2, volume in m^3, mean depth and tidal range in m
        self.__LWmean_depth = convert_len_to_prog_units(self.__LWmean_depth,
            input_len_units, self.name())
        self.__tidal_range = convert_len_to_prog_units(self.__tidal_range,
            input_len_units, self.name())
            
        self.__LWarea = convert_area_to_prog_units(self.__LWarea,
            input_area_units, self.name())
        self.__LWvolume = convert_vol_to_prog_units(self.__LWvolume,
            input_vol_units, self.name())
            
        if np.isfinite(self.__LWarea):
            # When the volume is not provided, it is derived from the mean depth
            # and the area - units are consistent: m^2 * m -> m^3
            if np.isinf(self.__LWvolume):
                self.__LWvolume = self.__LWarea*self.__LWmean_depth
            # When the mean depth is not provided, it is derived from the volume
            # and the area - units are consistent: m^3 / m^2 -> m
            if np.isinf(self.__LWmean_depth):
                self.__LWmean_depth = self.__LWvolume/self.__LWarea
                
        # Check if the GeoJSON file exists
        if self.__geojson_file:
            if not file_exists(self.__geojson_file):
                raise InputError(self.__geojson_file,
                                 "{}'s GeoJSON file does not exist".format(
                                    self.name()))
            else:
              # Load GeoJSON polygon
              with open(self.__geojson_file) as f:
                  js = json.load(f)
              self.__geojson_polygon = js['features'][0]
              
        assert(type(self.__existing_biomass) in [int, float])
        self.__existing_biomass = float(self.__existing_biomass)
        if self.__existing_biomass < 0.0:
            raise InputError(self.__existing_biomass,
                "Loch's existing biomass should be a positive number")
        self.__existing_biomass = convert_mass_to_prog_units(
            self.__existing_biomass, input_mass_units, self.name())
        
    def set_GeoJSON_properties(self, units='m'):
        if self.__geojson_polygon is not None:
            self.__geojson_polygon['properties']
    
    def aze(self, units='m^2'):
        """Allowable zone of effect (AZE) of the Loch
        
        In SEPA's 2008 Annex G: "Models for assessing the use of medicines in
        bath treatments", pp. 2-3, the AZE is defined as:
            AZE = min(0.5 km^2, 2% of loch area)
            
        This function returns: 2% of loch area
        
        Arguments:
            units: string; Output units. default: m^2 (program units)
        
        Reference:
            https://www.sepa.org.uk/media/113498/fish-farm-manual-annex-g.pdf
        """
        return self.area(units)*0.02
    
    def _compute_flushing_time(self):
        """Compute the flushing time of the Loch, T_F, in hours (program units)
        The flushing time of the Loch is considered to be infinite (very large
        water area when the volume is not set)
        
        References:
          - https://consultation.sepa.org.uk/permits/loch-long-salmon-beinn-
                reithe-car-application/user_uploads/210113-bnrt-par-nutrient-
                modelling_redacted.pdf, page 14
                
          - Gillibrand, P.A., Gubbins, G.J., Greathead, C., Davies, I.M. (2002),
            "Scottish executive locational guidelines for fish farming:
            predicted levels of nutrient enhancement and benthic impact",
            Scottish Fisheries Research Report Number 63 / 2002, page 5  
        """
        nhours_per_tiday_cycle=0.52*24.
        self.__flushing_time = nhours_per_tiday_cycle/0.7*self.__LWvolume/(
            self.__LWarea*self.__tidal_range)
        
    def set_flushing_time(self, value, units='h'):
        """Set the flushing time of the Loch, T_F
        
        Arguments:
            units: string; Input units. default: h (program units) 
        """
        self.__flushing_time = convert_time_to_prog_units(value, units,
                                                          self.name())
        
    def flushing_time(self, units='h'):
        """Return the flushing time of the Loch, T_F.
        The flushing may either have been evaluated using _compute_flushing_time
        or set using set_flushing_time
        
        Arguments:
            units: string; Output units. default: h (program units) 
        """
        return convert_time(self.__flushing_time, units, self.name())
        
    def flushing_rate(self, units='m^3/h'):
        """Return the flushing rate of the Loch, Q
        
        Arguments:
            units: string; Output units. default: m^3/h (program units)
        
        Reference:
          - Gillibrand, P.A., Gubbins, G.J., Greathead, C., Davies, I.M. (2002),
            "Scottish executive locational guidelines for fish farming:
            predicted levels of nutrient enhancement and benthic impact",
            Scottish Fisheries Research Report Number 63 / 2002, page 7  
        """
        q = self.__LWvolume/self.__flushing_time
        return convert_flowrate(q, units, self.name())
