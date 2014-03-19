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
