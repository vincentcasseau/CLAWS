#!/usr/bin/env python
"""Derived classes defining chemical substances.

Append a chemical to the already existing list of chemicals at the bottom of
this file or use the CustomSubstance class in the user scripts.
"""

# Import modules
import numpy as np
from claws.chemical import Chemical, ChemicalSubstance, Particle

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

class Dye(Chemical):
    pass

class Azamethiphos(ChemicalSubstance):
    def __init__(self, half_life=8.9, Loch=None, input_time_units='day',
                 reference='https://www.sepa.org.uk/regulations/water/'\
                           'aquaculture/environmental-standards/'):
        super().__init__(eqs_3hr=250.0, eqs_72hr=40.0, mac_72hr=100.0,
                         half_life=half_life, Loch=Loch,
                         input_time_units=input_time_units, reference=reference)
        
class Cypermethrin(ChemicalSubstance):
    def __init__(self, half_life=1.e10, Loch=None, input_time_units='day',
                 reference='https://www.sepa.org.uk/media/113498/'\
                         'fish-farm-manual-annex-g.pdf'):
        super().__init__(eqs_6hr=16.0, half_life=half_life, Loch=Loch,
                         input_time_units=input_time_units,
                         reference=reference)
        
class Deltamethrin(ChemicalSubstance):
    def __init__(self, half_life=1.e10, Loch=None, input_time_units='day',
                 reference='https://www.sepa.org.uk/regulations/water/'\
                           'aquaculture/environmental-standards/'):
        super().__init__(eqs_3hr=9.0, eqs_6hr=6.0, eqs_12hr=4.0, eqs_24hr=2.0,
                         eqs_48hr=1.0, half_life=half_life, Loch=Loch,
                         input_time_units=input_time_units, reference=reference)
                         
class FeedPellet(Particle):
    def __init__(self, diameter=0.0, grain_density=1127.0,
                 settling_velocity_mean=-0.032,
                 settling_velocity_deviation=0.0032, input_len_units='m',
                 input_conc_units='kg/m^3', input_vel_units='m/s',
                 reference='Cromey, C.J., Nickell, T.D., and Black, '\
                         'K.D., 2002, DEPOMOD – Modelling the Deposition and '\
                         'Biological Effects of Waste Solids From Marine Cage '\
                         'Farms. Aquaculture, 214: 211-239.'):
        super().__init__(diameter=diameter, grain_density=grain_density,
                         settling_velocity_mean=settling_velocity_mean,
                         settling_velocity_deviation=settling_velocity_deviation,
                         input_len_units=input_len_units,
                         input_conc_units=input_conc_units,
                         input_vel_units=input_vel_units, reference=reference)
