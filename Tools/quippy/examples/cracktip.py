from quippy import *

if len(sys.argv[1:]) != 3:
   print """Usage: cracktip RADIUS INFILE OUTFILE

RADIUS  - radius around cracktip (defined by <CrackPos,0,0)
INFILE  - input file, in XYZ or NetCDF format
OUTFILE - output file, in XYZ or NetCDF format
"""
   sys.exit(1)

radius = float(sys.argv[1])
infile = CInOutput(sys.argv[2], INPUT)
outfile = CInOutput(sys.argv[3], OUTPUT)

infile.query()

for fr in range(infile.n_frame):
   at = infile.read(frame=fr, zero=True)

   crackpos = farray([[at.params['CrackPos'],0.0,0.0]]).T
   subat = at.select(mask=((at.pos - crackpos).norm() < radius))

   outfile.write(subat)

   sys.stdout.write('%s: frame %d/%d -- selected %d/%d atoms.               \r' % (sys.argv[2], fr, infile.n_frame, subat.n, at.n))
   sys.stdout.flush()

print

infile.close()
outfile.close()
