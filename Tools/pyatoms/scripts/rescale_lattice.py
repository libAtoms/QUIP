#!/usr/bin/env python

from pyatoms import *
from numpy import *

import sys, os

factor = float(eval(sys.argv[1]))

for filename in sys.argv[2:]:
   root, ext = os.path.splitext(filename)
   a = Atoms(filename)

   a.lattice *= factor
   a.pos *= factor

   a.write_xyz(root+'_rescaled.xyz')
