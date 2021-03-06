#!/usr/bin/env python
"""Classes to define chemicals"""

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

class Chemical(object):
    def __init__(self, name="", half_life=1.e10, input_time_units='day',
                 reference=""):
        # __name: string; Name of the chemical. default is class name
        self.__name = name
        #__half_life: float; Chemical half-life in hours. default is 1.e10 (inf)
        self.__half_life = half_life
        # __reference: string; Link to input data. default: empty string
        self.__reference = reference
        
        self._sanitize_chemical(input_time_units)
        
    def __str__(self):
        return """{}\n\tHalf-life (h) = {}\n\tReference = {}""".format(
            self.name(), self.__half_life, self.__reference)
            
    def __eq__(self, other):
        if type(other).__name__ == type(self).__name__:
            return (self.__half_life == other.__half_life and
                    self.__reference == other.__reference)
        return False        
            
    def name(self):
        return self.__name
        
    def half_life(self, units='s'):
        return convert_time(self.__half_life, units, self.name())
    
    def _sanitize_chemical(self, input_time_units):
        if not self.__name:
            self.__name = ' '.join(re.findall('([A-Z][a-z]+)',
                                   type(self).__name__))
        
        # Half-life stored in hours
        assert(type(self.__half_life) is float)
        if self.__half_life <= 0.0:
            raise InputError(self.__half_life,
                             'Chemical half-life should be a positive number')
                
        self.__half_life = convert_time_to_prog_units(self.__half_life,
            input_time_units, self.name())
            
class ChemicalSubstance(Chemical):
    def __init__(self, name="", eqs_3hr=np.inf, eqs_6hr=np.inf, eqs_12hr=np.inf,
                 eqs_24hr=np.inf, eqs_48hr=np.inf, eqs_72hr=np.inf,
                 mac_72hr=np.inf, eqs_aze=0.5, mac_aze=0.5, half_life=1.e10,
                 Loch=None, input_conc_units='ng/L', input_area_units='km^2',
                 input_time_units='day', reference=""):
        super().__init__(name=name, half_life=half_life,
                         input_time_units=input_time_units, reference=reference)
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
        
        self._sanitize(input_conc_units, input_area_units)
        self._update_aze(Loch)
        
    def __str__(self):
        chemical_substance_str = """\n\t3-hour EQS (ng/L) = {}\
            \n\t72-hour EQS (ng/L) = {}\n\t72-hour EQS AZE (km^2) = {}\
            \n\t72-hour MAC (ng/L) = {}\n\t72-hour MAC AZE (km^2) = {}\
            """.format(self.__eqs_3hr, self.__eqs_72hr, self.eqs_aze('km^2'),
            self.__mac_72hr, self.mac_aze('km^2'))
        return super().__str__() + chemical_substance_str
            
    def __eq__(self, other):
        if type(other).__name__ == type(self).__name__:
            chemical_substance_eq_op = \
                (self.__eqs_3hr == other.__eqs_3hr and
                 self.__eqs_6hr == other.__eqs_6hr and
                 self.__eqs_12hr == other.__eqs_12hr and
                 self.__eqs_24hr == other.__eqs_24hr and
                 self.__eqs_48hr == other.__eqs_48hr and
                 self.__eqs_72hr == other.__eqs_72hr and
                 self.__eqs_aze == other.__eqs_aze and
                 self.__mac_72hr == other.__mac_72hr and
                 self.__mac_aze == other.__mac_aze)
            return super().__eq__(other) and chemical_substance_eq_op
        return False        
            
    def eqs_3hour(self, units='ng/L'):
        return convert_conc(self.__eqs_3hr, units, self.name())
        
    def eqs_6hour(self, units='ng/L'):
        return convert_conc(self.__eqs_6hr, units, self.name()) 
        
    def eqs_12hour(self, units='ng/L'):
        return convert_conc(self.__eqs_12hr, units, self.name()) 
        
    def eqs_24hour(self, units='ng/L'):
        return convert_conc(self.__eqs_24hr, units, self.name())
        
    def eqs_48hour(self, units='ng/L'):
        return convert_conc(self.__eqs_48hr, units, self.name())
        
    def eqs_72hour(self, units='ng/L'):
        return convert_conc(self.__eqs_72hr, units, self.name())
        
    def eqs_aze(self, units='m^2'):
        return convert_area(self.__eqs_aze, units, self.name())
      
    def mac_72hour(self, units='ng/L'):
        return convert_conc(self.__mac_72hr, units, self.name())
        
    def mac_aze(self, units='m^2'):
        return convert_area(self.__mac_aze, units, self.name())
    
    def _sanitize(self, input_conc_units, input_area_units):
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
                
class Particle(Chemical):
    def __init__(self, name="", diameter=0.0, grain_density=0.0,
                 settling_velocity_mean=0.0, settling_velocity_deviation=0.0,
                 half_life=1.e10, input_conc_units='kg/m^3',
                 input_len_units='m', input_time_units='day', 
                 input_vel_units='m/s', reference=""):
        super().__init__(name=name, half_life=half_life,
                         input_time_units=input_time_units, reference=reference)
        # __diameter: float; Particle diameter in meters. default is 0
        self.__diameter = diameter
        # __grain_density: float; Grain_density stored in ng/L. default is 0
        self.__grain_density = grain_density
        # __settling_velocity_mu: float; Still-water mean settling velocity,
        # in m/s. default is 0
        self.__settling_velocity_mu = settling_velocity_mean
        # ____settling_velocity_std: float; Deviation from the still-water
        # settling velocity, in m/s. default is 0
        self.__settling_velocity_std = settling_velocity_deviation
        
        self._sanitize(input_conc_units, input_len_units, input_vel_units)
        
    def __str__(self):
        particle_str = """\n\tDiameter (m) = {}\
            \n\tGrain density (kg/m^3) = {}\
            \n\tStill-water mean settling velocity (m/s) = {}\
            \n\tDeviation from the still-water settling velocity (m/s) = {}\
            """.format( self.__diameter, self.grain_density(),
            self.__settling_velocity_mu, self.__settling_velocity_std)
        return super().__str__() + particle_str
            
    def __eq__(self, other):
        if type(other).__name__ == type(self).__name__:
            particle_str = \
                (self.__diameter == other.__diameter and
                 self.__grain_density == other.__grain_density and
                 self.__settling_velocity_mu == other.__settling_velocity_mu and
                 self.__settling_velocity_std == other.__settling_velocity_std)
            return super().__eq__(other) and particle_str
        return False        
            
    def diameter(self, units='m'):
        return convert_len(self.__diameter, units, self.name())
        
    def grain_density(self, units='kg/m^3'):
        return convert_conc(self.__grain_density, units, self.name())
        
    def settling_velocity_mean(self, units='m/s'):
        return convert_vel(self.__settling_velocity_mu, units, self.name())
        
    def settling_velocity_deviation(self, units='m/s'):
        return convert_vel(self.__settling_velocity_std, units,
                           self.name())
    
    def _sanitize(self, input_conc_units, input_len_units, input_vel_units):
        assert(type(self.__diameter) is float)
        if self.__diameter <= 0.0:
            raise InputError(self.__diameter,
                'Diameter should be a positive number')
        self.__diameter = convert_len_to_prog_units(self.__diameter,
            input_len_units, self.name())
        
        assert(type(self.__grain_density) is float)
        if self.__grain_density <= 0.0:
            raise InputError(self.__grain_density,
                'Grain density should be a positive number')
        self.__grain_density = convert_conc_to_prog_units(self.__grain_density,
            input_conc_units, self.name())
            
        assert(type(self.__settling_velocity_mu) is float)
        self.__settling_velocity_mu = convert_vel_to_prog_units(
            self.__settling_velocity_mu, input_vel_units, self.name())
        
        assert(type(self.__settling_velocity_std) is float)
        if self.__settling_velocity_std < 0.0:
            raise InputError(self.__settling_velocity_std,
                'Still-water deviation from the settling velocity should be '\
                'positive or null')
        self.__settling_velocity_std = convert_vel_to_prog_units(
            self.__settling_velocity_std, input_vel_units, self.name())
