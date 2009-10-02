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
      self.assert_(cp._fpointer != self.dia._fpointer)
      self.assertEqual(self.dia.n, cp.n)
      self.assert_((self.dia.lattice == cp.lattice).all())
      self.assertEqual(self.dia.properties.keys(), cp.properties.keys())
      self.assertEqual(self.dia.params.keys(), cp.params.keys())
      self.assertEqual(self.dia.data, cp.data)

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
      self.assertRaises(ValueError, self.dia._get_n)
      self.assertRaises(ValueError, self.dia._set_n, 0)

   def testabort(self):
      self.assertRaises(RuntimeError, self.dia.read_xyz, '')

   def testoptional(self):
      self.dia.calc_connect() # without optional argument
      self.dia.calc_connect(1) # optional argument by position
      self.dia.calc_connect(own_neighbour=1)


if __name__ == '__main__':
   unittest.main()
