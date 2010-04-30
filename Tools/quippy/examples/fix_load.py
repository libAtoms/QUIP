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

from pylab import *
from quippy import *
import sys

if len(sys.argv[1:]) != 2:
   print 'Usage: fix_load <crack_relax_loading_movie> <outfile>'
   sys.exit(1)

al = AtomsList(sys.argv[1])

df2 = al.df2[:]

# Find first minimum in df2
min1 = 2
while df2[min1] <= df2[1]:
   min1 += 1
min1 -= 1
relaxed = al[min1].copy()
relaxed.load[:] = al[min1].pos - al[-1].pos

relaxed.write(sys.argv[2])

