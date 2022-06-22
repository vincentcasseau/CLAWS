#!/usr/bin/env python
"""Test script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"

from claws.lochs import *

def print_loch_info(loch_obj):
    print(loch_obj)
    print("\n{}'s allowable zone of effect (AZE) (km^2): {:.6f}".format(
          loch_obj.name(), loch_obj.aze('km^2')))
    print("{}'s flushing rate (M m^3/yr): {:.6f}".format(loch_obj.name(),
          loch_obj.flushing_rate('M m^3/yr')))
    
if __name__ == "__main__":
    print_loch_info(LochCreran())
    print("")
    print_loch_info(UndefinedLoch('Firth of Clyde'))
