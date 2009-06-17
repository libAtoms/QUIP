from quippy import *
import unittest


class TestOOFortran(unittest.TestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)

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
      cp = self.dia.copy()
      cp.cutoff = 1.3
      self.assertEqual(cp.cutoff, 1.3)

   def testroutine(self):
      self.assert_(self.dia.has_property('pos'))

   def testroutinekwargs(self):
      cp = self.dia.copy()
      cp.calc_connect(own_neighbour=True)

   def testroutine2(self):
      sup = supercell(self.dia, 2, 2, 2)
      ds = DynamicalSystem(sup)
      self.assertRaises(TypeError, DynamicalSystem, ds)
      self.assertRaises(TypeError, DynamicalSystem, atoms_in=ds)

   def testinterface(self):
      cp = self.dia.copy()
      self.assertRaises(TypeError, cp.add_atoms, 1)
      cp.add_atoms([1.0,1.0,1.0], 1)
      self.assertEqual(cp.n, 9)
      self.assertEqual(cp.z[9], 1)
      self.assertEqual(list(cp.pos[9]), [1.0, 1.0, 1.0])

   def testinterfacekwargs(self):
      cp = self.dia.copy()
      cp.add_atoms(pos=farray([1.0,1.0,1.0]), z=1)
      self.assertEqual(cp.n, 9)
      self.assertEqual(cp.z[9], 1)
      self.assertEqual(list(cp.pos[9]), [1.0, 1.0, 1.0])
      
   def testfpointernone(self):
      cp = self.dia.copy()
      cp._fpointer = None
      self.assertRaises(ValueError, cp._get_n)
      self.assertRaises(ValueError, cp._set_n, 0)


def getTestSuite():
   return unittest.TestLoader().loadTestsFromTestCase(TestOOFortran)

if __name__ == '__main__':
   suite = getTestSuite()
   unittest.TextTestRunner(verbosity=2).run(suite)

