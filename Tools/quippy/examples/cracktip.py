from quippy import *

if len(sys.argv[1:]) != 3:
   print """Usage: cracktip RADIUS INFILE OUTFILE

RADIUS  - radius around cracktip (defined by <CrackPos,0,0)
INFILE  - input file, in XYZ or NetCDF format
OUTFILE - output file, in XYZ or NetCDF format
"""
   sys.exit(1)

def filter_cracktip(frames, radius):
   for at in frames:
      crackpos = farray([[at.params['CrackPos'],0.0,0.0]])
      yield at.select(mask=((at.pos - crackpos).norm() < radius))

frames    = AtomsList(sys.argv[2])
tipframes = AtomsList(filter_cracktip(frames, float(sys.argv[1])))
tipframes.write(sys.argv[3], progress=True)
