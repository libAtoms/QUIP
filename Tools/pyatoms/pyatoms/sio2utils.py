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

spind  = {'Si': 2, 'O': 1}
revspind = dict(zip(spind.values(), spind.keys()))

def write_pos_cel(f1, f2=None, posfilename='pos.in', celfilename='cel.in', dt=1.0):

   if f2 is None:
      f2 = f1.copy()
      f2.pos += f1.velo*dt

   if f1.n != f2.n:
      print 'frame1.n (%d) != frame2.n (%d)' % (f1.n, f2.n)
      sys.exit(1)

   it = 0

   posf = open(posfilename, 'w')
   celf = open(celfilename, 'w')

   objind = 1

   for (step_name, at) in (('MSTEP', f1), ('0STEP', f2)):
      posf.write(' %s %d\n' % (step_name, it))
      for i in range(at.n):
         p = at.pos[i,:]/BOHR
         posf.write('%20.10e%20.10e%20.10e%4d%4d X\n' % (p[0], p[1], p[2], spind[at.species[i]], objind))

      celf.write(' %s %d\n' % (step_name, it))
      for a in at.lattice:
         L = a/BOHR
         celf.write('%20.10e%20.10e%20.10e\n' % (L[0], L[1], L[2]))

   posf.close()
   celf.close()


def read_pos_cel(posfilename='pos.in', celfilename='cel.in'):
   pos = open(posfilename).readlines()
   cel = open(celfilename).readlines()

   n = (len(pos)-2)/2

   lattice1 = array([ [float(x) for x in L.split()] for L in cel[1:4] ])*BOHR
   lattice2 = array([ [float(x) for x in L.split()] for L in cel[5:8] ])*BOHR

   a1 = Atoms(n=n, lattice=lattice1)
   a2 = Atoms(n=n, lattice=lattice2)

   pos1 = array([ [float(x) for x in L.split()[0:3] ] for L in pos[1:n+1] ])*BOHR
   pos2 = array([ [float(x) for x in L.split()[0:3] ] for L in pos[n+2:2*n+2] ])*BOHR

   species1 = array([ revspind[int(L.split()[3])] for L in pos[1:n+1] ])
   species2 = array([ revspind[int(L.split()[3])] for L in pos[n+2:2*n+2] ])

   if not all(species1 == species2):
      print 'Mismatch in species between MSTEP and 0STEP'
      sys.exit(1)

   a1.pos[:] = pos1
   a1.species[:] = species1

   a2.pos[::] = pos2
   a2.species[:] = species2

   return a1, a2




