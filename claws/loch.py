#!/usr/bin/env python
"""Root class defining a Loch
"""

# Import modules
from os.path import exists as file_exists
import re
import json
import numpy as np

from claws.custom_exceptions import InputError
from claws.helpers import *

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
        return self.__geojson_polygon['geometry']
        
    def geojson_polygon_properties(self):
        return self.__geojson_polygon['properties']
    
    def existing_biomass(self, units='tonne'):
        return convert_mass(self.__existing_biomass, units, self.name())  
        
    def set_geojson_polygon_properties(self, input_dict):
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
