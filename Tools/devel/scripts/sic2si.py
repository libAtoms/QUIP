#!/usr/bin/env python

from atoms import *
import sys

a = Atoms(cell=sys.argv[1])

a.lattice = a.lattice * (5.43/4.385)
a.element[:] = 'Si'

for i in range(a.n):
   a.pos[i,:] = a.pos[i,:] * (5.43/4.385)

a.write(sys.argv[2]+'.xyz')

cell = CastepCell()
a.write_cell(cell)

cell.write(sys.argv[2]+'.cell')
