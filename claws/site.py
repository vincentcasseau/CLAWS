#!/usr/bin/env python
"""Root class defining a site.
"""

# Import modules
import re
import numpy as np
from pyproj import Proj

from claws.helpers import indent
from claws.custom_exceptions import InputError

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

class Site(object):
    def __init__(self, longitude=None, latitude=None, x=None, y=None,
                 proj_params=None, reference="", name=""):
        # self.__name: string; Site name. default is class name
        self.__name = name
        # __lon: float; Site longitude. default is None
        self.__lon = longitude
        # __lat: float; Site latitude. default is None
        self.__lat = latitude
        # __lon: float; Site x-position. default is None
        self.__x = x
        # __lat: float; Site y-position. default is None
        self.__y = y
        # __reference: string; Reference. default: empty string
        self.__reference = reference
        
        self._sanitize(proj_params)
        
    def __str__(self, indent_lvl=1):
        if self.__input_coordinate_system == 'Geographic':
            return """{1}\n{0}Longitude = {2}\n{0}Latitude = {3}\
                \n{0}Reference = {4}""".format(indent(indent_lvl), self.name(),
                self.__lon, self.__lat, self.__reference)
        else:
            return """{1}\n{0}x-position = {2}\n{0}y-position = {3}\
                \n{0}Reference = {4}""".format(indent(indent_lvl), self.name(),
                self.__x, self.__y, self.__reference)
    
    def coordinate_system(self):
        return self.__input_coordinate_system
        
    def lon(self):
        return self.__lon
        
    def lat(self):
        return self.__lat
        
    def x(self):
        return self.__x
        
    def y(self):
        return self.__y
      
    def name(self):
        return self.__name
        
    def position(self):
        if self.__input_coordinate_system == 'Geographic':
            return (self.__lon, self.__lat)
        return (self.__x, self.__y)
        
    def _sanitize(self, proj_params):
        if not self.__name:
            self.__name = ' '.join(re.findall('([A-Z][a-z]+)',
                                   type(self).__name__))
        
        if self.__lon is None and self.__lat is not None:
            raise InputError(self.__lon, 
                "Longitude of Site {} must be given when the latitude is"
                    .format(self.name()))
        elif self.__lat is None and self.__lon is not None:
            raise InputError(self.__lat, 
                "Latitude of Site {} must be given when the longitude is"
                    .format(self.name()))
                    
        if self.__x is None and self.__y is not None:
            raise InputError(self.__x, 
                "x-position of Site {} must be given when the y-position is"
                    .format(self.name()))
        elif self.__y is None and self.__x is not None:
            raise InputError(self.__y, 
                "y-position of Site {} must be given when the x-position is"
                    .format(self.name()))
        
        if ((self.__x is None and self.__y is None) and 
                (self.__lon is None and self.__lat is None)):
             raise InputError(self.__lon, 
                  "Site {}: Please provide the site location either as (x,y) "\
                      "or as (lon,lat)".format(self.name()))          
        
        if self.__lon is not None:
            self.__input_coordinate_system = 'Geographic'
            
            assert(type(self.__lon) in [int, float])
            if abs(self.__lon) > 90.0:
                raise InputError(self.__lon, 
                    "Longitude of Site {} must be comprised between -90 deg "\
                        "and +90 deg".format(self.name()))
            
            assert(type(self.__lat) in [int, float])
            if abs(self.__lat) > 90.0:
                raise InputError(self.__lat, 
                    "Latitude of Site {} must be comprised between -90 deg "\
                        "and +90 deg".format(self.name()))
        
        if self.__x is not None:
            self.__input_coordinate_system = 'Cartesian'
            
            assert(type(self.__x) in [int, float])
            assert(type(self.__y) in [int, float])
            
            if proj_params is None:
                raise InputError("proj_params", 
                  "Site {}: Please provide projection parameters to the site "\
                      "location when defined as (x,y)".format(self.name())) 
            
            # Projection and transformation
            xy2lonlat = Proj(proj_params, preserve_units=False)
            self.__lon, self.__lat = xy2lonlat(self.__x, self.__y, inverse=True)
