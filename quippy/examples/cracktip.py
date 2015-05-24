# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from quippy import *

if len(sys.argv[1:]) != 3:
   print """Usage: cracktip RADIUS INFILE OUTFILE

RADIUS  - radius around cracktip (defined by <CrackPos,0,0>)
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
