"""PyAtoms package

(c) James Kermode <jrk33@cam.ac.uk> 2008
    http://www.srcf.ucam.org/~jrk33/PyAtoms"""

import atoms
from atoms import *

import units
from units import *

import castep

from paramreader import ParamReader
from ordereddict import OrderedDict

# util module is distinctly out of date at the moment!
#import util 


__all__ = ['Atoms','OrderedDict','ParamReader','frame_reader','norm',
           'norm2','rms_diff','count','diamond','castep']

# export all the constants in the units module
for u in units.__dict__:
   if not u.startswith('__'):
      __all__.append(u)

