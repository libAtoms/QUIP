#!/usr/bin/env python

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
