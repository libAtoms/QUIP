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

from numpy import *
from pyatoms import *

import sys

selection_edge_tol = 10.0

if len(sys.argv) != 3:
   print 'Usage: %s <old file.xyz> <new file.xyz>' % sys.argv[0]
   sys.exit(1)

at = Atoms(sys.argv[1])

print 'Found properties: %r' % at.properties.keys()

# Rename some properties
try:
   at.properties.rename('embed_mask','embed')
   at.properties.rename('fit_mask','fit')
   at.properties.rename('constrain_mask', 'move_mask')
   at.repoint()
except ValueError, message:
   print message
   sys.exit(1)

# Negate move_mask since constrain_mask was other way round
at.move_mask[:] = 1 - at.move_mask[:]

print 'Fixing %d atoms' % count(at.move_mask == 0)

at.add_property('old_nn',at.nn)
at.add_property('md_old_changed_nn',at.changed_nn)
at.add_property('edge_mask',0)

# Set edge_mask to 1 for atoms closer than edge_mask_tol to an edge

minx = at.pos[:,0].min() + selection_edge_tol
maxx = at.pos[:,0].max() - selection_edge_tol
miny = at.pos[:,1].min() + selection_edge_tol
maxy = at.pos[:,1].max() - selection_edge_tol

at.edge_mask[at.pos[:,0] < minx] = 1
at.edge_mask[at.pos[:,0] > maxx] = 1
at.edge_mask[at.pos[:,1] < miny] = 1
at.edge_mask[at.pos[:,1] > maxy] = 1

print 'Setting edge_mask=1 for %d atoms' % count(at.edge_mask == 1)

required_properties = 'pos:embed:fit:nn:changed_nn:old_nn:md_old_changed_nn:edge_mask:load'.split(':')

for p in required_properties:
   if not p in at.properties.keys():
      print '%s missing' % p

at.write(sys.argv[2])
