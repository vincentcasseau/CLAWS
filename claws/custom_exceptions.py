#!/usr/bin/env python
"""Defines custom exceptions.
"""

__author__ = "Vincent Casseau and Tom Scanlon"
__copyright__ = "Copyright 2022, MTS-CFD Ltd."
__credits__ = ["Vincent Casseau", "Tom Scanlon"]
__license__ = "GPLv2"
__version__ = "1.2"
__maintainer__ = "Vincent Casseau"
__email__ = "claws.scot@gmail.com"
__status__ = "Production"


# ---------------------------------------------------------------------------- #
# Custom exceptions 
# ---------------------------------------------------------------------------- #

class InputError(Exception):
    """Exception raised to deal with user-input errors.

    Attributes:
        expression: string; input expression in which the error occurred
        message: string; explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
