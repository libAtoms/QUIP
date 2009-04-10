#!/usr/bin/env python

from pyatoms import *
import sys

a = Atoms(cell=sys.argv[1])

a.lattice = a.lattice * (3.56683/5.43)
a.species[:] = 'C'

for i in range(a.n):
   a.pos[i,:] = a.pos[i,:] * (3.56683/5.43)

a.write_xyz(sys.argv[2]+'.xyz')

cell = castep.CastepCell()
a.to_cell(cell)

cell.write(sys.argv[2]+'.cell')
