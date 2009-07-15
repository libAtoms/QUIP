from quippy import *
import unittest, itertools, sys, quippy
from quippytest import *

class TestAtomsList(QuippyTestCase):

   def setUp(self):
      self.listsrc = [diamond(5.44+0.01*x, 14) for x in range(5)]
      self.listal = AtomsList(self.listsrc)
      
      self.gensrc = (diamond(5.44+0.01*x,14) for x in range(5))
      self.genal = AtomsList(self.gensrc)
      
   def testgetitem(self):
      a0 = self.listal[0]
      a4 = self.listal[4]
      self.assertEqual(a0, diamond(5.44,14))
      self.assertEqual(a4, diamond(5.44+0.04,14))

   def testlazygetitem(self):
      # before items are accessed, AtomsList._list = None
      self.assert_(self.listal._list[0] is None)
      a0 = self.listal[0]
      self.assertEqual(self.listal._list[0], a0)
      self.assert_(self.listal._list[0] is a0)

   def testgetslice(self):
      sl = self.listal[1:]
      self.assert_(isinstance(sl, AtomsList))
      self.assertEqual(list(sl), list(self.listal)[1:])

   def testgetitemnegindex(self):
      last = self.listal[-1]
      self.assert_(self.listal._list[-1] is last)

   def testdelitem(self):
      self.listal.loadall()
      del self.listal[0]
      self.assert_(self.listal._list[0] is None)
      self.assertEqual(self.listal[0], diamond(5.44, 14))
      
   def testgetitemoutofrange(self):
      self.assertRaises(IndexError, self.listal.__getitem__, 6)

   def testtuple(self):
      tupal = AtomsList(tuple(self.listsrc))
      self.assertEqual(list(tupal), list(self.listal))

   def testloadall(self):
      self.listal.loadall()
      nonlazy = AtomsList([diamond(5.44+0.01*x, 14) for x in range(5)], lazy=False)
      self.assertEqual(self.listal._list, nonlazy._list)

   def testtolist(self):
      L1 = list(self.listal)
      L2 = list(self.listal)
      self.assertEqual(L1, L2)
      self.assert_(L1[0] is L2[0])
      self.assert_(L1 is not L2)

   def testreversed(self):
      self.assertEqual(list(reversed(self.listal)), list(reversed(list(self.listal))))

   def testgen(self):
      self.assertEqual(len(self.genal), 0)
      a0 = self.genal[0]
      self.assertEqual(len(self.genal), 1)
      a1 = self.genal[1]
      self.assertEqual(len(self.genal), 2)
      i = iter(self.genal)
      self.assert_(a0 is i.next())
      self.assert_(a1 is i.next())
      self.assertEqual(len(self.genal), 2)
      a2 = i.next()
      self.assertEqual(len(self.genal), 3)

   def testgentolist(self):
      self.assertEqual(list(self.genal), list(self.listal))

   def testgennonlazy(self):
      nonlazy = AtomsList((diamond(5.44+0.01*x,14) for x in range(5)), lazy=False)
      self.assertEqual(nonlazy._list, list(self.genal))

   def testgenloadall(self):
      self.assertEqual(len(self.genal), 0)
      self.genal.loadall()
      self.assertEqual(len(self.genal), 5)      
      self.assertEqual(list(self.genal), self.genal._list)

   def testgengetitem(self):
      a0 = self.genal[0]
      a4 = self.genal[4]
      self.assertEqual(a0, diamond(5.44,14))
      self.assertEqual(a4, diamond(5.44+0.04,14))

   def testgengetslice(self):
      sl = self.genal[0:2]
      self.assert_(isinstance(sl, AtomsList))
      self.assertEqual(list(sl), list(self.genal)[0:2])

   def testgengetsliceopen(self):
      sl = self.genal[0:]
      self.assert_(isinstance(sl, AtomsList))
      self.assertEqual(list(sl), list(self.genal)[0:])

   def testgengetitemoutofrange(self):
      self.assertRaises(IndexError, self.genal.__getitem__, 5)

   def testrandomaccess(self):
      self.assert_(not self.genal._randomaccess)
      self.assert_(self.listal._randomaccess)

   def testgenreverse(self):
      revgen = reversed(self.genal)
      self.assertRaises(ValueError, list, revgen)
      L1 = list(self.genal)
      L2 = list(reversed(self.genal))
      self.assertEqual(list(reversed(L1)), L2)

   def testwrite(self):

      def testwriter():
         # Dummy atoms write generator which returns number of atoms
         at = yield None
         while True:
            at = yield at.n

      g = testwriter()
      g.next()
      self.assertEqual(self.listal.write(g), [8, 8, 8, 8, 8])

      g = testwriter()
      g.next()
      self.assertEqual(self.genal.write(g),  [8, 8, 8, 8, 8])


   def testfromat(self):
      nl = AtomsList(self.listal[0])
      self.assertEqual(list(nl), [self.listal[0]])


      
class TestAtomsListCInOutput(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44,14), 2,2,2)
      self.al = AtomsList([ supercell(diamond(5.44+0.01*x,14),2,2,2) for x in range(5) ])

      self.xyz_ref = ['64\n',
                      'Lattice="10.880000 0.000000 0.000000 0.000000 10.880000 0.000000 0.000000 0.000000 10.880000" Properties=species:S:1:pos:R:3:Z:I:1\n',
                      'Si              0.00000000      0.00000000      0.00000000      14\n',
                      'Si              1.36000000      1.36000000      1.36000000      14\n',
                      'Si              2.72000000      2.72000000      0.00000000      14\n',
                      'Si              4.08000000      4.08000000      1.36000000      14\n',
                      'Si              2.72000000      0.00000000      2.72000000      14\n',
                      'Si              4.08000000      1.36000000      4.08000000      14\n',
                      'Si              0.00000000      2.72000000      2.72000000      14\n',
                      'Si              1.36000000      4.08000000      4.08000000      14\n',
                      'Si              0.00000000      0.00000000      5.44000000      14\n',
                      'Si              1.36000000      1.36000000      6.80000000      14\n',
                      'Si              2.72000000      2.72000000      5.44000000      14\n',
                      'Si              4.08000000      4.08000000      6.80000000      14\n',
                      'Si              2.72000000      0.00000000      8.16000000      14\n',
                      'Si              4.08000000      1.36000000      9.52000000      14\n',
                      'Si              0.00000000      2.72000000      8.16000000      14\n',
                      'Si              1.36000000      4.08000000      9.52000000      14\n',
                      'Si              0.00000000      5.44000000      0.00000000      14\n',
                      'Si              1.36000000      6.80000000      1.36000000      14\n',
                      'Si              2.72000000      8.16000000      0.00000000      14\n',
                      'Si              4.08000000      9.52000000      1.36000000      14\n',
                      'Si              2.72000000      5.44000000      2.72000000      14\n',
                      'Si              4.08000000      6.80000000      4.08000000      14\n',
                      'Si              0.00000000      8.16000000      2.72000000      14\n',
                      'Si              1.36000000      9.52000000      4.08000000      14\n',
                      'Si              0.00000000      5.44000000      5.44000000      14\n',
                      'Si              1.36000000      6.80000000      6.80000000      14\n',
                      'Si              2.72000000      8.16000000      5.44000000      14\n',
                      'Si              4.08000000      9.52000000      6.80000000      14\n',
                      'Si              2.72000000      5.44000000      8.16000000      14\n',
                      'Si              4.08000000      6.80000000      9.52000000      14\n',
                      'Si              0.00000000      8.16000000      8.16000000      14\n',
                      'Si              1.36000000      9.52000000      9.52000000      14\n',
                      'Si              5.44000000      0.00000000      0.00000000      14\n',
                      'Si              6.80000000      1.36000000      1.36000000      14\n',
                      'Si              8.16000000      2.72000000      0.00000000      14\n',
                      'Si              9.52000000      4.08000000      1.36000000      14\n',
                      'Si              8.16000000      0.00000000      2.72000000      14\n',
                      'Si              9.52000000      1.36000000      4.08000000      14\n',
                      'Si              5.44000000      2.72000000      2.72000000      14\n',
                      'Si              6.80000000      4.08000000      4.08000000      14\n',
                      'Si              5.44000000      0.00000000      5.44000000      14\n',
                      'Si              6.80000000      1.36000000      6.80000000      14\n',
                      'Si              8.16000000      2.72000000      5.44000000      14\n',
                      'Si              9.52000000      4.08000000      6.80000000      14\n',
                      'Si              8.16000000      0.00000000      8.16000000      14\n',
                      'Si              9.52000000      1.36000000      9.52000000      14\n',
                      'Si              5.44000000      2.72000000      8.16000000      14\n',
                      'Si              6.80000000      4.08000000      9.52000000      14\n',
                      'Si              5.44000000      5.44000000      0.00000000      14\n',
                      'Si              6.80000000      6.80000000      1.36000000      14\n',
                      'Si              8.16000000      8.16000000      0.00000000      14\n',
                      'Si              9.52000000      9.52000000      1.36000000      14\n',
                      'Si              8.16000000      5.44000000      2.72000000      14\n',
                      'Si              9.52000000      6.80000000      4.08000000      14\n',
                      'Si              5.44000000      8.16000000      2.72000000      14\n',
                      'Si              6.80000000      9.52000000      4.08000000      14\n',
                      'Si              5.44000000      5.44000000      5.44000000      14\n',
                      'Si              6.80000000      6.80000000      6.80000000      14\n',
                      'Si              8.16000000      8.16000000      5.44000000      14\n',
                      'Si              9.52000000      9.52000000      6.80000000      14\n',
                      'Si              8.16000000      5.44000000      8.16000000      14\n',
                      'Si              9.52000000      6.80000000      9.52000000      14\n',
                      'Si              5.44000000      8.16000000      8.16000000      14\n',
                      'Si              6.80000000      9.52000000      9.52000000      14\n']

   def tearDown(self):
      if os.path.exists('test.xyz'): os.remove('test.xyz')
      if os.path.exists('test.nc'): os.remove('test.nc')
      if os.path.exists('test.xyz.idx'): os.remove('test.xyz.idx')
      if os.path.exists('test2.xyz'): os.remove('test2.xyz')
      if os.path.exists('test2.xyz.idx'): os.remove('test2.xyz.idx')

   def testsinglexyz(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz')
      self.assertEqual(self.at, at)
      self.assertEqual(self.xyz_ref, open('test.xyz', 'r').readlines())

   def testsinglenc(self):
      self.at.write('test.nc')
      at = Atoms('test.nc')
      self.assertEqual(self.at, at)
      
   def testmultixyz(self):
      self.al.write('test.xyz')
      al = AtomsList('test.xyz')
      self.assertEqual(len(al), 5)
      self.assertEqual(len(self.al), len(al))
      self.assertEqual(list(self.al), list(al))

   def testmultinc(self):
      self.al.write('test.nc')
      al = AtomsList('test.nc')
      self.assertEqual(list(self.al), list(al))      

   def testxyzlowlevel(self):
      cio = CInOutput("test.xyz", OUTPUT, append=False)
      for a in self.al:
         cio.write(a)
      cio.close()

      cio = CInOutput("test.xyz")
      a = []
      for i in range(5):
         a.append(cio.read())
      self.assertEqual(a, list(self.al))      

   def testnclowlevel(self):
      cio = CInOutput("test.nc", OUTPUT)
      for a in self.al:
         cio.write(a)
      cio.close()

      cio = CInOutput("test.nc")
      a = []
      for i in range(5):
         a.append(cio.read())
      self.assertEqual(a, list(self.al))      

   def teststdinstdout(self):
      self.al.write('test.xyz')
      os.system("python -c 'from quippy import *; al = AtomsList(\"stdin\"); al.write(\"stdout\",)' < test.xyz > test2.xyz")
      al = AtomsList('test2.xyz')
      self.assertEqual(list(al), list(self.al))


class TestPythonNetCDF(QuippyTestCase):
   def setUp(self):
      self.at = supercell(diamond(5.44,14), 2,2,2)
      self.al = AtomsList([ supercell(diamond(5.44+0.01*x,14),2,2,2) for x in range(5) ])
      self.al.write('test3.nc', netcdf4=False)

   def tearDown(self):
      if os.path.exists('test3.nc'): os.remove('test3.nc')
      if os.path.exists('dataset.nc'): os.remove('dataset.nc')

   def testpupynere_read(self):
      from pupynere import netcdf_file
      nc = netcdf_file('test3.nc', 'r')
      al = AtomsList(nc, format=quippy.netcdf_file)
      self.assertEqual(list(self.al), list(al))
      nc.close()
      
   def testnetcdf4_read(self):
      from netCDF4 import Dataset
      nc = Dataset('test3.nc','r')
      al = AtomsList(nc)
      self.assertEqual(list(self.al), list(al))
      nc.close()

   def testnetcdf4_write(self):
      from netCDF4 import Dataset
      nc = Dataset('dataset.nc','w')
      al2 = AtomsList(self.al)
      al2.write(nc)
      nc.close()
      al = AtomsList('dataset.nc')
      self.assertEqual(list(self.al), list(al))

class TestNetCDFAtomsList(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 2, 2, 2)
      self.al = AtomsList(self.at for x in range(5))
      self.al.write('test.nc')

   def tearDown(self):
      os.remove('test.nc')

   def testread(self):
      at = netcdf_read('test.nc').next()
      self.assertEqual(at, self.at)

   def testgetvars(self):
      nc = NetCDFAtomsList('test.nc')
      self.assertEqual(nc.variables.keys(), ['cell_angular', 'cell_lengths', 'coordinates', 'spatial', 'cell_spatial', 'Z', 'cell_angles', 'species'])


def getTestSuite():
   tl = unittest.TestLoader()
   return unittest.TestSuite([tl.loadTestsFromTestCase(TestAtomsList),
                              tl.loadTestsFromTestCase(TestAtomsListCInOutput),
                              tl.loadTestsFromTestCase(TestPythonNetCDF),
                              tl.loadTestsFromTestCase(TestNetCDFAtomsList)])

if __name__ == '__main__':
   suite = getTestSuite()
   unittest.TextTestRunner(verbosity=2).run(suite)
