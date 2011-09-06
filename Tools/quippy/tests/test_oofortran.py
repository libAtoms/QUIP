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
import unittest
from quippytest import *

class TestOOFortran(QuippyTestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)

   def tearDown(self):
      del self.dia

   def testdia(self):
      self.assertEqual(self.dia.n, 8)
      self.assert_((self.dia.z == 14).all())
      self.assert_((self.dia.species.stripstrings() == 'Si').all())

   def testcopy(self):
      cp = self.dia.copy()
      self.assert_(not all(cp._fpointer == self.dia._fpointer))
      self.assertEqual(self.dia.n, cp.n)
      self.assert_((self.dia.lattice == cp.lattice).all())
      self.assertEqual(self.dia.properties.keys(), cp.properties.keys())
      self.assertEqual(self.dia.params.keys(), cp.params.keys())

   def testgetset(self):
      self.dia.cutoff = 1.3
      self.assertEqual(self.dia.cutoff, 1.3)

   def testroutine(self):
      self.assert_(self.dia.has_property('pos'))

   def testroutinekwargs(self):
      self.dia.calc_connect(own_neighbour=True)

   def testroutine2(self):
      sup = supercell(self.dia, 2, 2, 2)
      ds = DynamicalSystem(sup)
      self.assertRaises(TypeError, DynamicalSystem, ds)
      self.assertRaises(TypeError, DynamicalSystem, atoms_in=ds)

   def testinterface(self):
      self.assertRaises(TypeError, self.dia.add_atoms, 1)
      self.dia.add_atoms([1.0,1.0,1.0], 1)
      self.assertEqual(self.dia.n, 9)
      self.assertEqual(self.dia.z[9], 1)
      self.assertEqual(list(self.dia.pos[9]), [1.0, 1.0, 1.0])

   def testinterfacekwargs(self):
      self.dia.add_atoms(pos=farray([1.0,1.0,1.0]), z=1)
      self.assertEqual(self.dia.n, 9)
      self.assertEqual(self.dia.z[9], 1)
      self.assertEqual(list(self.dia.pos[9]), [1.0, 1.0, 1.0])
      
   def testfpointernone(self):
      self.dia._fpointer = None
      self.assertRaises(ValueError, getattr, self.dia, 'n')
      self.assertRaises(ValueError, setattr, self.dia, 'n', 0)

   def testabort(self):
      import quippy._atoms
      self.assertRaises(RuntimeError, quippy._atoms.Atoms.select, self.dia, self.dia)

   def testoptional(self):
      self.dia.calc_connect() # without optional argument
      self.dia.calc_connect(self.dia.connect) # optional argument by position
      self.dia.calc_connect(own_neighbour=1)


if __name__ == '__main__':
   unittest.main()
