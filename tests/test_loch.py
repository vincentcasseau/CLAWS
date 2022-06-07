#!/usr/bin/env python
"""Test script"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "hystrath@gmail.com"
__status__ = "Production"

from claws.lochs import *

if __name__ == "__main__":
    loch = LochCreran()
    print(loch)
    print("{}'s allowable zone of effect (AZE) (km^2): {:.6f}".format(
          loch.name(), loch.aze('km^2')))
    print("{}'s flushing time (day): {:.6f}".format(loch.name(),
          loch.flushing_time('day')))
    print("\n{}'s flushing rate (M m^3/yr): {:.6f}".format(loch.name(),
          loch.flushing_rate('M m^3/yr')))
