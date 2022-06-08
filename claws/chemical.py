#!/usr/bin/env python
"""Abstract class to define a chemical substance in terms of
    - Environmental Quality Standard (EQS)
    - Maximum Allowable Concentration (MAC)
    - Allowable Zone of Effect (AZE)
"""

# Import modules
import re
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

class ChemicalSubstance(object):
    def __init__(self, name="", eqs_3hr=np.inf, eqs_6hr=np.inf, eqs_12hr=np.inf,
                 eqs_24hr=np.inf, eqs_48hr=np.inf, eqs_72hr=np.inf,
                 mac_72hr=np.inf, eqs_aze=0.5, mac_aze=0.5, half_life=1.e10,
                 Loch=None, input_conc_units='ng/L', input_area_units='km^2',
                 input_time_units='day', reference=""):
        # self.__name: string; Name of the chemical substance. default is empty
        # string, ie., name of the derived class
        self.__name = name
        # __eqs_Xhr: float; X-hour Environmental Quality Standard (EQS) in ng/L.
        # default: infinity
        self.__eqs_3hr = eqs_3hr
        self.__eqs_6hr = eqs_6hr
        self.__eqs_12hr = eqs_12hr
        self.__eqs_24hr = eqs_24hr
        self.__eqs_48hr = eqs_48hr
        self.__eqs_72hr = eqs_72hr
        # __eqs_aze: float; Allowable Zone of Effect (AZE) in m^2 for the
        # 72-hour EQS. default: 0.5 km^2 = 5.e5 m^2
        self.__eqs_aze = eqs_aze
        # __mac_72hr: float; Maximum Allowable Concentration (MAC) standard in
        # ng/L. default: infinity
        self.__mac_72hr = mac_72hr
        # __mac_aze: float; Allowable Zone of Effect (AZE) in m^2 for the
        # 72-hour MAC. default: 0.5 km^2 = 5.e5 m^2
        self.__mac_aze = mac_aze
        #__half_life: float; Chemical half-life in hours. default is 1.e10 (inf)
        self.__half_life = half_life
        # __reference: string; Link to input data. default: empty string
        self.__reference = reference
        
        self._sanitize(input_conc_units, input_area_units, input_time_units)
        self._update_aze(Loch)
        
    def __str__(self):
        return """{}\n\t3-hour EQS (ng/L) = {}\n\t72-hour EQS (ng/L) = {}\
            \n\t72-hour EQS AZE (km^2) = {}\n\t72-hour MAC (ng/L) = {}\
            \n\t72-hour MAC AZE (km^2) = {}\n\tHalf-life (h) = {}
            \n\tReference = {}""".format(self.name(), self.__eqs_3hr,
            self.__eqs_72hr, self.eqs_aze('km^2'), self.__mac_72hr,
            self.mac_aze('km^2'), self.__half_life, self.__reference)
            
    def __eq__(self, other):
        if type(other).__name__ == type(self).__name__:
            return (self.__eqs_3hr == other.__eqs_3hr and
                    self.__eqs_6hr == other.__eqs_6hr and
                    self.__eqs_12hr == other.__eqs_12hr and
                    self.__eqs_24hr == other.__eqs_24hr and
                    self.__eqs_48hr == other.__eqs_48hr and
                    self.__eqs_72hr == other.__eqs_72hr and
                    self.__eqs_aze == other.__eqs_aze and
                    self.__mac_72hr == other.__mac_72hr and
                    self.__mac_aze == other.__mac_aze and
                    self.__half_life == other.__half_life and
                    self.__reference == other.__reference)
        return False        
            
    def eqs_3hour(self):
        return self.__eqs_3hr
        
    def eqs_6hour(self):
        return self.__eqs_6hr 
        
    def eqs_12hour(self):
        return self.__eqs_12hr  
        
    def eqs_24hour(self):
        return self.__eqs_24hr
        
    def eqs_48hour(self):
        return self.__eqs_48hr 
        
    def eqs_72hour(self):
        return self.__eqs_72hr
        
    def eqs_aze(self, units='m^2'):
        return convert_area(self.__eqs_aze, units, self.name())
      
    def mac_72hour(self):
        return self.__mac_72hr
        
    def mac_aze(self, units='m^2'):
        return convert_area(self.__mac_aze, units, self.name())
        
    def half_life(self, units='s'):
        return convert_time(self.__half_life, units, self.name())
      
    def name(self):
        return self.__name
    
    def _sanitize(self, input_conc_units, input_area_units, input_time_units):
        if not self.__name:
            self.__name = ' '.join(re.findall('([A-Z][a-z]+)',
                                   type(self).__name__))
        
        if np.isfinite(self.__eqs_3hr):
            assert(type(self.__eqs_3hr) is float)
            if self.__eqs_3hr <= 0.0:
                raise InputError(self.__eqs_3hr,
                    '3-hour EQS should be a positive number')
            self.__eqs_3hr = convert_conc_to_prog_units(self.__eqs_3hr,
                input_conc_units, self.name())
                
        if np.isfinite(self.__eqs_6hr):
            assert(type(self.__eqs_6hr) is float)
            if self.__eqs_6hr <= 0.0:
                raise InputError(self.__eqs_6hr,
                    '6-hour EQS should be a positive number')
            self.__eqs_6hr = convert_conc_to_prog_units(self.__eqs_6hr,
                input_conc_units, self.name())
                
        if np.isfinite(self.__eqs_12hr):
            assert(type(self.__eqs_12hr) is float)
            if self.__eqs_12hr <= 0.0:
                raise InputError(self.__eqs_12hr,
                    '12-hour EQS should be a positive number')
            self.__eqs_12hr = convert_conc_to_prog_units(
                self.__eqs_12hr, input_conc_units, self.name())
                
        if np.isfinite(self.__eqs_24hr):
            assert(type(self.__eqs_24hr) is float)
            if self.__eqs_24hr <= 0.0:
                raise InputError(self.__eqs_24hr,
                    '24-hour EQS should be a positive number')
            self.__eqs_24hr = convert_conc_to_prog_units(
                self.__eqs_24hr, input_conc_units, self.name())
                
        if np.isfinite(self.__eqs_48hr):
            assert(type(self.__eqs_48hr) is float)
            if self.__eqs_48hr <= 0.0:
                raise InputError(self.__eqs_48hr,
                    '48-hour EQS should be a positive number')
            self.__eqs_48hr = convert_conc_to_prog_units(
                self.__eqs_48hr, input_conc_units, self.name())
                
        if np.isfinite(self.__eqs_72hr):
            assert(type(self.__eqs_72hr) is float)
            if self.__eqs_72hr <= 0.0:
                raise InputError(self.__eqs_72hr,
                    '72-hour EQS should be a positive number')
            self.__eqs_72hr = convert_conc_to_prog_units(
                self.__eqs_72hr, input_conc_units, self.name())
                
        if np.isfinite(self.__mac_72hr):
            assert(type(self.__mac_72hr) is float)
            if self.__mac_72hr <= 0.0:
                raise InputError(self.__mac_72hr,
                    '72-hour MAC should be a positive number')
            self.__mac_72hr = convert_conc_to_prog_units(
                self.__mac_72hr, input_conc_units, self.name())
                
        assert(type(self.__eqs_aze) is float)
        if self.__eqs_aze <= 0.0:
            raise InputError(self.__eqs_aze,
                '72-hour EQS Allowable Zone of Effect (AZE) '\
                'should be a positive number')        
                
        assert(type(self.__mac_aze) is float)
        if self.__mac_aze <= 0.0:
            raise InputError(self.__mac_aze,
                '72-hour MAC Allowable Zone of Effect (AZE) '\
                'should be a positive number')        
                
        # Areas are stored in m^2
        self.__eqs_aze = convert_area_to_prog_units(self.__eqs_aze,
            input_area_units, self.name())
        self.__mac_aze = convert_area_to_prog_units(self.__mac_aze,
            input_area_units, self.name())
        # These AZEs have yet to be corrected using a Loch's AZE    
        self.__aze_updated = False
        
        # Half-life stored in hours
        assert(type(self.__half_life) is float)
        if self.__half_life <= 0.0:
            raise InputError(self.__half_life,
                             'Chemical half-life should be a positive number')
                
        self.__half_life = convert_time_to_prog_units(self.__half_life,
            input_time_units, self.name())
        
    def _update_aze(self, loch_obj=None):
        """Allowable zone of effect (AZE) in m^2 (program units)
        
        In SEPA's 2008 Annex G: "Models for assessing the use of medicines in
        bath treatments", pp. 2-3, the AZE is defined as:
            AZE = min(0.5 km^2, 2% of loch area)
            
        Reference:
            https://www.sepa.org.uk/media/113498/fish-farm-manual-annex-g.pdf
            
        Arguments:
            loch_obj: Loch object. default is None    
        """
        # To prevent multiple updates
        if not self.__aze_updated:
            # If the Loch object is not passed as an argument, self.__eqs_aze
            # isn't edited
            if loch_obj is not None:
                self.__eqs_aze = min(self.__eqs_aze, loch_obj.aze())
                self.__mac_aze = min(self.__eqs_aze, loch_obj.aze())
                self.__aze_updated = True
