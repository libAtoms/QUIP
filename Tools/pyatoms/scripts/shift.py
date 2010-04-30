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

# centre on a given atom

from numpy import *

def center_on_atom(at, ref):
   ref_pos = at.pos[ref,:]
   for i in range(at.n):
      if at.pos[i,0] > ref_pos[0]:
         at.pos[i,:] = at.pos[i,:] - ref_pos
      else:
         at.pos[i,:] = at.pos[i,:] + array((at.lattice[0,0]-ref_pos[0],-ref_pos[1],-ref_pos[2]))


def shift(at, shift):
   for i in range(at.n):
      at.pos[i,:] += shift
