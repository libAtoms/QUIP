"""Simple example of looping over an input file, adding properties
   and writing to an output file.
"""

from quippy import *
import sys

if len(sys.argv[1:]) != 2:
   print 'Usage: addprop.py <INPUT> <OUTPUT>'
   sys.exit(1)

fr = FrameReader(sys.argv[1])
outf = CInOutput(sys.argv[2], OUTPUT)

at0 = fr.next() # first frame
for at in fr:   # rest of frames
  at.add_property('p', 0.0, n_cols=3)
  at.p[:] = at0.pos[:]
  outf.write(at)

outf.close()
