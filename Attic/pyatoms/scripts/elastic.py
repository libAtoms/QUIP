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

def elastic_energy(at):
   # assumes periodic along z, open surfaces on x and y
   w = at.pos[:,0].max() - at.pos[:,0].min()
   h = at.pos[:,1].max() - at.pos[:,1].min()
   d = at.lattice[2,2]
   vol_per_atom = w*h*d/at.n # A^3/atom
   e = zeros(at.n)
   eps = zeros(6)
   sig = zeros(6)
   for i in range(at.n):
      eps[:] = (at.S_xx_sub1[i],at.S_yy_sub1[i],at.S_zz_sub1[i],
                at.S_xy[i], at.S_xz[i], at.S_yz[i])

      sig[:] = (at.Sig_xx[i], at.Sig_yy[i], at.Sig_zz[i],
                at.Sig_xy[i], at.Sig_xz[i], at.Sig_yz[i])

      e[i] = 0.5*dot(sig,eps)

   return e*vol_per_atom/GPA # convert to eV/atom


def add_elastic_energy(at):
   at.add_property('ee',elastic_energy(at))


def elastic_region(at,R,tip_pos):
   return logical_not(sqrt((at.pos[:,0]-tip_pos[0])**2 + (at.pos[:,1]-tip_pos[1])**2) < R)
 
def sweep_elastic_energy(at1,at2,tip_pos,Rmax,dR=1):
   data = []
   for R in arange(0.0,Rmax,dR):
      data.append( (R, at1.ee[elastic_region(at1,R,tip_pos)].sum() -
                       at2.ee[elastic_region(at2,R,tip_pos)].sum() ) )
      
   return array(data)

def sweep_elastic_energy_tot(at1,tip_pos,Rmax,dR=1):
   data = []
   for R in arange(0.0,Rmax,dR):
      data.append( (R, at1.ee[elastic_region(at1,R,tip_pos)].sum()))
      
   return array(data)

def elastic_constant_matrix(C11,C12,C44):
   c = zeros((6,6),'d')
   c[0,0] = c[1,1] = c[2,2] = C11
   c[0,1] = c[1,0] = c[0,2] = c[2,0] = c[1,2] = c[2,1] = C12
   c[3,3] = c[4,4] = c[5,5] = C44
   return c
