# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from quippy import *
from numpy import *

import unittest, quippy
from quippytest import *

class TestDynamicalSystem(QuippyTestCase):
   def setUp(self):
       d = diamond(5.44, 14)
       self.at = supercell(d, 2, 2, 2)
       self.ds = DynamicalSystem(self.at)

   def test_atoms(self):
       self.assertEqual(self.ds.atoms, self.at)

   def test_atoms_pointer(self):
       self.assert_(all(self.ds.atoms._fpointer == self.at._fpointer))

   def test_avgpos(self):
       self.assertArrayAlmostEqual(self.ds.atoms.avgpos, self.ds.atoms.pos)

   def test_avg_ke(self):
       self.assertArrayAlmostEqual(self.ds.atoms.avg_ke, 0.5*self.at.mass*self.at.velo.norm2())


if __name__ == '__main__':
   unittest.main()
