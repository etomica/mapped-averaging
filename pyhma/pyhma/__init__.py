###############################################################################
# pyhma: A Python library for HMA method 
# 
# Copyright (c) 2019
# 
# Authors: Sabry Moustafa, Andrew Schultz, and David Kofke 
# 
# pyhma is a free software
###############################################################################


"""The pyhma Python package is an implementation of the harmonically-mapped averaging (HMA) method
to measure anharmonic properties of crystalline systems from MD codes (VASP, currently) outputs.

"""

__author__ = "Sabry Moustafa, Andrew Schultz, and David Kofke"
__license__ = "Mozilla Public License"
__email__ = "sabrygad@buffalo.edu, ajs42@buffalo.edu, kofke@buffalo.edu"

from pyhma.vasp_reader   import read
from pyhma.processor     import Processor

