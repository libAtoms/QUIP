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

