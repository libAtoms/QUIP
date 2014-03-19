"""PyAtoms package

(c) James Kermode <jrk33@cam.ac.uk> 2008
    http://www.srcf.ucam.org/~jrk33/PyAtoms"""
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

