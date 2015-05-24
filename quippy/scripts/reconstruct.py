#!/usr/bin/env python
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
from numpy import *
import operator, sys, os

def find_tip(crack, ytol, bondlength):

#   w = where(logical_and(crack.nn == 3, crack.edge_mask == 0))[0]

   # Find left-most fully coordinated atoms closer than ytol to y=0
   w = where(logical_and(abs(crack.pos[2,:]) < ytol,crack.nn == 4))[0]
   L = [ (crack.pos[1,w[i]], crack.pos[2,w[i]], w[i]) for i in frange(len(w)) ]
   L.sort()

   def bonded(pair):
      (xi, yi, i), (xj, yj, j) = pair
      sep = crack.pos[i] - crack.pos[j]
      return sqrt(dot(sep,sep)) < bondlength

   def bottom(a,b):
      (xi, yi, i), (xj, yj, j) = a, b
      if yi > yj:
         return a
      else:
         return b

   def top(a,b):
      (xi, yi, i), (xj, yj, j) = a, b
      if yi > yj:
         return b
      else:
         return a

   pairs = [ (bottom(a,b),top(a,b)) for a, b in zip(L[::2],L[1::2]) ] 
   pairs = filter(bonded, pairs)

   tip_lower, tip_upper = pairs[0][0][2], pairs[0][1][2]

   print 'Crack tip is %d-%d bond' % (tip_lower, tip_upper)
   return tip_lower, tip_upper


def reconstruct(crack, tip_lower, tip_upper, bondlength):

   crack.set_cutoff(bondlength)
   crack.calc_connect()
   #neartip = neighbours(crack, tip_lower, 3.0*bondlength, remove_self=False)

#   def below_left(neighb):
#      i, dist, diff = neighb
#      return diff[0] <= 0.0 and diff[1] >= 0.0#
#
#   neartip = filter(below_left, neartip)

   # calc connectivity in this region
   #neighb = dict([ (i,neighbours(crack, i, bondlength))
   #                for i in map(lambda x: x[0],neartip) ])

   # locate atoms a, b

   idx = {}
   
   # start at tip_lower, try to step down and left

   for neighb in crack.neighbours[tip_lower]:
      if neighb.diff[1] < 0.0 and neighb.diff[2] > 0.0:  break
   else:      raise ValueError("Can't find atom b")
   idx['b'] = neighb.j

   # start at b, try to step left
   for neighb in crack.neighbours[idx['b']]:
      print neighb.j,neighb.distance,neighb.diff
      if neighb.diff[1] < -bondlength*.4: break
   else:
      raise ValueError("Can't find atom a")
   idx['a'] = neighb.j

   disp = {}
   disp['a'] = array([.7*bondlength/2.6,.8*bondlength/2.6,0])
   disp['b'] = array([0,-.8*bondlength/2.6,0])

   for k in idx:
      print k, idx[k], disp[k]
      crack.pos[:,idx[k]] += disp[k]

if __name__ == '__main__':
   
   if (len(sys.argv) < 3):
      print 'Usage: %s inputfile.xyz outputfile.xyz [tip_lower tip_upper]' % sys.argv[0]
      sys.exit(1)

   tip_lower = tip_upper = None
   if len(sys.argv) > 3:
      tip_lower = int(sys.argv[3])
      tip_upper = int(sys.argv[4])
      print 'Crack tip manually set to %d--%d bond' % (tip_lower, tip_upper)

   crack = Atoms(sys.argv[1])

   if all(crack.z == 14):
      bondlength = 3.0 # silicon
   elif all(crack.z == 6):
      bondlength = 1.8
   else:
      print "Don't know bondlength for this system"
      sys.exit(1)

   if tip_lower is None or tip_upper is None:
      tip_lower, tip_upper = find_tip(crack,ytol=2.0,bondlength=bondlength)
   reconstruct(crack, tip_lower, tip_upper, bondlength=bondlength)
   crack.write(sys.argv[2])
   
