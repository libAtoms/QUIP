#!/usr/bin/env python
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from pyatoms import *
from numpy import *
import operator, sys, os

def neighbours(at, i, c, remove_self=True):
   p = at.pos[i,:]
   diff = at.pos - p
   dist = norm(diff)
   neighb = where(dist < c)[0]
   if remove_self:
      neighb = neighb[neighb != i]

   return zip(neighb, dist[neighb], diff[neighb])


def find_tip(crack, ytol, bondlength):

#   w = where(logical_and(crack.nn == 3, crack.edge_mask == 0))[0]

   # Find left-most fully coordinated atoms closer than ytol to y=0
   w = where(logical_and(abs(crack.pos[:,1]) < ytol,crack.nn == 4))[0]
   L = [ (crack.pos[w[i],0], crack.pos[w[i],1], w[i]) for i in range(len(w)) ]
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

   neartip = neighbours(crack, tip_lower, 3.0*bondlength, remove_self=False)

#   def below_left(neighb):
#      i, dist, diff = neighb
#      return diff[0] <= 0.0 and diff[1] >= 0.0#
#
#   neartip = filter(below_left, neartip)

   # calc connectivity in this region
   neighb = dict([ (i,neighbours(crack, i, bondlength))
                   for i in map(lambda x: x[0],neartip) ])

   # locate atoms a, b

   idx = {}
   
   # start at tip_lower, try to step down and left

   for i,dist,diff in neighb[tip_lower]:
      if diff[0] < 0.0 and diff[1] > 0.0:  break
   else:      raise ValueError("Can't find atom b")
   idx['b'] = i

   # start at b, try to step left
   for i,dist,diff in neighb[idx['b']]:
      print i,dist,diff
      if diff[0] < -bondlength*.4: break
   else:
      raise ValueError("Can't find atom a")
   idx['a'] = i

   disp = {}
   disp['a'] = array([.7*bondlength/2.6,.8*bondlength/2.6,0])
   disp['b'] = array([0,-.8*bondlength/2.6,0])

   for k in idx:
      print k, idx[k], disp[k]
      crack.pos[idx[k],:] += disp[k]

   # Make new QM zone
   crack.embed[:] = 0
   qmzone = neighbours(crack, tip_lower, 3.5*bondlength, remove_self=False)
   for i,dist,diff in qmzone:
      crack.embed[i] = 1
   print 'Got %d QM atoms' % count(crack.embed==1)


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

   unrec = crack.copy() # save a copy

   if all(crack.species == 'Si'):
      bondlength = 3.0 # silicon
   elif all(crack.species == 'C'):
      bondlength = 1.8
   else:
      print "Don't know bondlength for this system"
      sys.exit(1)

   if tip_lower is None or tip_upper is None:
      tip_lower, tip_upper = find_tip(crack,ytol=2.0,bondlength=bondlength)
   reconstruct(crack, tip_lower, tip_upper, bondlength=bondlength)
   crack.write_xyz(sys.argv[2])

   # copy recon embed list back into original file
   unrec.embed[:] = crack.embed

   print 'count(crack.embed) = ', count(crack.embed)
   print 'count(unrec.embed) = ', count(unrec.embed)

   unrec.write_xyz(sys.argv[1])
   
