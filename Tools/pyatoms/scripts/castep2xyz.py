from pyatoms import *
import sys

xyz = open(sys.argv[1],'w')
files = sys.argv[2:]

i = 0
for f in files:
   try:
      a = castep.read_castep_output(f)
      a.write_xyz(xyz)
      print i
      i += 1
   except ValueError:
      continue
