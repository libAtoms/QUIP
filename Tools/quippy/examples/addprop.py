"""Simple example of looping over an input file, adding properties
   and writing to an output file.
"""

from quippy import *
import sys

if len(sys.argv[1:]) != 2:
   print 'Usage: addprop.py <INPUT> <OUTPUT>'
   sys.exit(1)

frames = AtomsList(sys.argv[1])

for at in frames[1:]:
  at.add_property('p', 0.0, n_cols=3)
  at.p[:] = frames[0].pos[:]

(frames[1:]).write(sys.argv[2])
