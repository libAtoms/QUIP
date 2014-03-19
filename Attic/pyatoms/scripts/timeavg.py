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

# read in a bunch of frames and do a linear average of all their real properties

from pyatoms import *

def temperature(a):
   return a.v2.mean()*a.mass[0]/(3.0*BOLTZMANN_K)

def add_speed(frames):
   for a in frames:
      a.add_property('v',norm(a.velo))
      a.add_property('v2',norm2(a.velo))
      yield a

def timeavg(frames):
   first = frames.next()
   avg = first.copy()

   n = 1
   for a in frames:
      print n, temperature(a), temperature(avg)
      avg.real[:] = (n*avg.real + a.real)/float(n+1.0)
      n = n + 1

   return avg
