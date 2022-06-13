#!/usr/bin/env python
"""Derived classes defining Lochs.

Append a Loch to the already existing list of Lochs at the bottom of this file
or use the UndefinedLoch class in the absence of available data.

Lochs' data can be found in
    Edwards, A. and Sharples, F. (1986), Scottish Sea Lochs: a Catalogue,
    Scottish Marine Biological Association/Nature Conservancy Council

UndefinedLoch assumes that the Loch's area is infinitely large. As such, the
allowable zone of effect (AZE) as defined by SEPA

    AZE = min(0.5 km^2, 2% of loch area)

will always equal 0.5 km^2.      
        
Reference:
    https://www.sepa.org.uk/media/113498/fish-farm-manual-annex-g.pdf
"""

# Import modules
from claws.loch import Loch

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Variable 
# ---------------------------------------------------------------------------- #

_refEdwardsSharples = 'Edwards, A. and Sharples, F. (1986), '\
                      '"Scottish Sea Lochs: a Catalogue", Scottish Marine '\
                      'Biological Association / Nature Conservancy Council, '\

# ---------------------------------------------------------------------------- #
# Classes 
# ---------------------------------------------------------------------------- #

class UndefinedLoch(Loch):
    def __init__(self, name=""):
        super().__init__(name=name)
        
class LochAlsh(Loch):
    def __init__(self):
        super().__init__(area=27.5, tidal_range=4.6, volume=1134.3,
                         mean_depth=41.2,
                         reference=_refEdwardsSharples + "p. 47") 
            
class LochCreran(Loch):
    def __init__(self):
        super().__init__(area=13.3, tidal_range=3.3, volume=177.6,
                         mean_depth=13.4,
                         reference=_refEdwardsSharples + "p. 84") 
                         
class LochDuich(Loch):
    def __init__(self):
        super().__init__(area=11.6, tidal_range=4.6, volume=512.2,
                         mean_depth=44.2,
                         reference=_refEdwardsSharples + "p. 93")
            
class LochHourn(Loch):
    def __init__(self):
        super().__init__(area=33.7, tidal_range=4.2, volume=2005.5,
                         mean_depth=59.5,
                         reference=_refEdwardsSharples + "p. 133")
            
class LochLinnhe(Loch):
    def __init__(self):
        super().__init__(area=31.7, tidal_range=3.7, volume=1344.7,
                         mean_depth=42.5,
                         reference=_refEdwardsSharples + "p. 151")
                         
class LochLong(Loch):
    def __init__(self):
        super().__init__(area=44.0, tidal_range=3.1, volume=1758.0,
                         mean_depth=40.0,
                         reference="Gillibrand, P.A., Gubbins, G.J., "\
                         "Greathead, C., Davies, I.M. (2002), Scottish "\
                         "executive locational guidelines for fish farming: "\
                         "predicted levels of nutrient enhancement and benthic"\
                         " impact, Scottish Fisheries Research Report Number "\
                         "63 / 2002")

class LochMelfort(Loch):
    def __init__(self):
        super().__init__(area=9.3, tidal_range=2.3, volume=260.5,
                         mean_depth=27.9,
                         reference=_refEdwardsSharples + "p. 160")
