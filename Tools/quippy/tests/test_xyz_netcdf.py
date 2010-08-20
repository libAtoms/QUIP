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
from quippy.xyz import *
from quippy.netcdf import *
import unittest, itertools, sys, quippy
from quippytest import *


class TestCInOutput(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44,14), 2,2,2)
      self.at.add_property('log', False)
      self.at.params['real'] = 1.0
      self.at.params['int'] = 2
      self.at.params['int_a'] = [1,2,3]
      self.at.params['real_a'] = [1.0,2.0,3.0]
      self.at.params['int_a2'] = farray([1,2,3,4,5,6,7,8,9]).reshape(3,3)
      self.at.params['real_a2'] = farray([1.0,2,3,4,5,6,7,8,9]).reshape(3,3)
      self.at.params['log_param'] = True
      self.at.params['log_a'] = [True, True, False]
      self.at.params['string'] = 'string'
      self.at.params['string2'] = 'string with spaces'
      self.al = AtomsList([ supercell(diamond(5.44+0.01*x,14),2,2,2) for x in range(5) ])
      for at in self.al:
         at.params.update(self.at.params)

      self.xyz_ref =  ['64\n', 'real=1.00000000 int=2 int_a="1       2       3" real_a="1.00000000      2.00000000      3.00000000" int_a2="1        4        7        2        5        8        3        6        9" real_a2="1.00000000       4.00000000       7.00000000       2.00000000       5.00000000       8.00000000       3.00000000       6.00000000       9.00000000" log_param=T log_a="T T F" string=string string2="string with spaces" Lattice="10.88000000       0.00000000       0.00000000       0.00000000      10.88000000       0.00000000       0.00000000       0.00000000      10.88000000" Properties=species:S:1:pos:R:3:Z:I:1:log:L:1\n',
                       'Si              0.00000000      0.00000000      0.00000000      14    F\n',
                       'Si              1.36000000      1.36000000      1.36000000      14    F\n',
                       'Si              2.72000000      2.72000000      0.00000000      14    F\n',
                       'Si              4.08000000      4.08000000      1.36000000      14    F\n',
                       'Si              2.72000000      0.00000000      2.72000000      14    F\n',
                       'Si              4.08000000      1.36000000      4.08000000      14    F\n',
                       'Si              0.00000000      2.72000000      2.72000000      14    F\n',
                       'Si              1.36000000      4.08000000      4.08000000      14    F\n',
                       'Si              0.00000000      0.00000000      5.44000000      14    F\n',
                       'Si              1.36000000      1.36000000      6.80000000      14    F\n',
                       'Si              2.72000000      2.72000000      5.44000000      14    F\n',
                       'Si              4.08000000      4.08000000      6.80000000      14    F\n',
                       'Si              2.72000000      0.00000000      8.16000000      14    F\n',
                       'Si              4.08000000      1.36000000      9.52000000      14    F\n',
                       'Si              0.00000000      2.72000000      8.16000000      14    F\n',
                       'Si              1.36000000      4.08000000      9.52000000      14    F\n',
                       'Si              0.00000000      5.44000000      0.00000000      14    F\n',
                       'Si              1.36000000      6.80000000      1.36000000      14    F\n',
                       'Si              2.72000000      8.16000000      0.00000000      14    F\n',
                       'Si              4.08000000      9.52000000      1.36000000      14    F\n',
                       'Si              2.72000000      5.44000000      2.72000000      14    F\n',
                       'Si              4.08000000      6.80000000      4.08000000      14    F\n',
                       'Si              0.00000000      8.16000000      2.72000000      14    F\n',
                       'Si              1.36000000      9.52000000      4.08000000      14    F\n',
                       'Si              0.00000000      5.44000000      5.44000000      14    F\n',
                       'Si              1.36000000      6.80000000      6.80000000      14    F\n',
                       'Si              2.72000000      8.16000000      5.44000000      14    F\n',
                       'Si              4.08000000      9.52000000      6.80000000      14    F\n',
                       'Si              2.72000000      5.44000000      8.16000000      14    F\n',
                       'Si              4.08000000      6.80000000      9.52000000      14    F\n',
                       'Si              0.00000000      8.16000000      8.16000000      14    F\n',
                       'Si              1.36000000      9.52000000      9.52000000      14    F\n',
                       'Si              5.44000000      0.00000000      0.00000000      14    F\n',
                       'Si              6.80000000      1.36000000      1.36000000      14    F\n',
                       'Si              8.16000000      2.72000000      0.00000000      14    F\n',
                       'Si              9.52000000      4.08000000      1.36000000      14    F\n',
                       'Si              8.16000000      0.00000000      2.72000000      14    F\n',
                       'Si              9.52000000      1.36000000      4.08000000      14    F\n',
                       'Si              5.44000000      2.72000000      2.72000000      14    F\n',
                       'Si              6.80000000      4.08000000      4.08000000      14    F\n',
                       'Si              5.44000000      0.00000000      5.44000000      14    F\n',
                       'Si              6.80000000      1.36000000      6.80000000      14    F\n',
                       'Si              8.16000000      2.72000000      5.44000000      14    F\n',
                       'Si              9.52000000      4.08000000      6.80000000      14    F\n',
                       'Si              8.16000000      0.00000000      8.16000000      14    F\n',
                       'Si              9.52000000      1.36000000      9.52000000      14    F\n',
                       'Si              5.44000000      2.72000000      8.16000000      14    F\n',
                       'Si              6.80000000      4.08000000      9.52000000      14    F\n',
                       'Si              5.44000000      5.44000000      0.00000000      14    F\n',
                       'Si              6.80000000      6.80000000      1.36000000      14    F\n',
                       'Si              8.16000000      8.16000000      0.00000000      14    F\n',
                       'Si              9.52000000      9.52000000      1.36000000      14    F\n',
                       'Si              8.16000000      5.44000000      2.72000000      14    F\n',
                       'Si              9.52000000      6.80000000      4.08000000      14    F\n',
                       'Si              5.44000000      8.16000000      2.72000000      14    F\n',
                       'Si              6.80000000      9.52000000      4.08000000      14    F\n',
                       'Si              5.44000000      5.44000000      5.44000000      14    F\n',
                       'Si              6.80000000      6.80000000      6.80000000      14    F\n',
                       'Si              8.16000000      8.16000000      5.44000000      14    F\n',
                       'Si              9.52000000      9.52000000      6.80000000      14    F\n',
                       'Si              8.16000000      5.44000000      8.16000000      14    F\n',
                       'Si              9.52000000      6.80000000      9.52000000      14    F\n',
                       'Si              5.44000000      8.16000000      8.16000000      14    F\n',
                       'Si              6.80000000      9.52000000      9.52000000      14    F\n']

   def tearDown(self):
      if os.path.exists('test.xyz'): os.remove('test.xyz')
      if os.path.exists('test.nc'): os.remove('test.nc')
      if os.path.exists('test.xyz.idx'): os.remove('test.xyz.idx')
      if os.path.exists('test2.xyz'): os.remove('test2.xyz')
      if os.path.exists('test2.xyz.idx'): os.remove('test2.xyz.idx')
      if os.path.exists('quartz.xyz'): os.remove('quartz.xyz')
      if os.path.exists('quartz.xyz.idx'): os.remove('quartz.xyz.idx')
      if os.path.exists('quartz.nc'): os.remove('quartz.nc')

   def testsinglexyz(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz')
      self.assertAtomsEqual(self.at, at)
      self.assertEqual(self.xyz_ref, open('test.xyz', 'r').readlines())

   def testsinglexyzprefix(self):
      self.at.write('test.xyz', prefix='PREFIX')
      lines = open('test.xyz').readlines()
      self.assert_(all([line[:len('PREFIX')] == 'PREFIX' for line in lines]))
      lines_without_prefix = [line[len('PREFIX '):] for line in lines]
      at = Atoms(''.join(lines_without_prefix), format='string')
      self.assertAtomsEqual(self.at, at)

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

   def testmultixyzprefix(self):
      self.al.write('test.xyz', prefix='PREFIX')
      lines = open('test.xyz').readlines()
      self.assert_(all([line[:len('PREFIX')] == 'PREFIX' for line in lines]))

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

   def testxyzlowlevel2(self):
      cio = CInOutput("test.xyz", OUTPUT, append=False)
      for a in self.al:
         a.write(cio)
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

   def testwritecio(self):
      cio = CInOutput("test2.xyz", OUTPUT)
      self.al.write(cio)
      cio.close()
      al = AtomsList("test2.xyz")
      self.assertEqual(list(al), list(self.al))

   def testreadcio(self):
      self.al.write("test.xyz")
      cio = CInOutput("test.xyz", INPUT)
      al = AtomsList(cio)
      self.assertEqual(list(al), list(self.al))
      cio.close()

   def testframe_random_access(self):
      self.al.write("test.xyz")
      cio = CInOutput("test.xyz", INPUT)
      error = farray(0)
      at = cio.read(frame=4)
      self.assertArrayAlmostEqual(at.lattice, 2*(5.44+0.01*4)*fidentity(3))
      cio.close()

   def testframe_out_of_range(self):
      self.al.write("test.xyz")
      cio = CInOutput("test.xyz", INPUT)
      self.assertRaises(RuntimeError, cio.read, frame=5)
      cio.close()

   def testwrite_single_xyz_properties(self):
      self.at.write('test.xyz', properties=['species','pos'])
      at = Atoms('test.xyz')
      self.assertEqual(sorted(at.properties.keys()), sorted(['species', 'pos', 'Z']))
      self.assertArrayAlmostEqual(self.at.pos, at.pos)
      self.assertEqual(list(self.at.z), list(at.z))

   def testwrite_multi_xyz_properties(self):
      self.al.write('test.xyz', properties=['species','pos'])
      al = AtomsList('test.xyz')
      for at,at_ref in zip(al, self.al):
         self.assertEqual(sorted(at.properties.keys()), sorted(['species', 'pos', 'Z']))
         self.assertArrayAlmostEqual(at.pos, at_ref.pos)
         self.assertEqual(list(at.z), list(at_ref.z))

   def test_non_orthorhombic_xyz(self):
      from quippy.sio2 import quartz_params
      aq1 = alpha_quartz(**quartz_params['ASAP_JRK'])
      aq1.write('quartz.xyz')
      aq2 = Atoms('quartz.xyz')
      orig_params = get_lattice_params(aq1.lattice)
      xyz_params  = get_lattice_params(aq2.lattice)
      self.assertArrayAlmostEqual(orig_params, xyz_params)

   def test_non_orthorhombic_nc(self):
      from quippy.sio2 import quartz_params
      aq1 = alpha_quartz(**quartz_params['ASAP_JRK'])
      aq1.map_into_cell()
      aq1.write('quartz.nc', netcdf4=False)
      aq2 = Atoms('quartz.nc')
      orig_params = get_lattice_params(aq1.lattice)
      nc_params   = get_lattice_params(aq2.lattice)
      self.assertArrayAlmostEqual(orig_params, nc_params)

   def test_read_string(self):
      s = ''.join(self.xyz_ref)
      cio = CInOutput()
      at = cio.read(str=s)
      self.assertAtomsEqual(at, self.at)

   def test_read_ext_string(self):
      es = Extendable_str(''.join(self.xyz_ref))
      cio = CInOutput()
      at = cio.read(estr=es)
      self.assertAtomsEqual(at, self.at)

   def test_read_loop(self):
      import resource
      max_open_files_soft, max_open_files_hard = resource.getrlimit(resource.RLIMIT_NOFILE)
      self.at.write('test.xyz')
      for i in range(2*max_open_files_soft):
         a = Atoms('test.xyz', frame=0)

   def test_write_loop(self):
      import resource
      max_open_files_soft, max_open_files_hard = resource.getrlimit(resource.RLIMIT_NOFILE)
      for i in range(2*max_open_files_soft):
         self.at.write('test.xyz')

   def test_read_bad_range_1(self):
      self.at.write('test.xyz')
      self.assertRaises(RuntimeError, Atoms, 'test.xyz', range=[1,self.at.n+1])

   def test_read_bad_range_2(self):
      self.at.write('test.xyz')
      self.assertRaises(RuntimeError, Atoms, 'test.xyz', range=[12, 10])

   def test_read_xyz_range_all(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', range=[1,64])
      self.assertAtomsEqual(at, self.at)

   def test_read_xyz_range_subset(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', range=[1,32])
      sub = self.at.select(list=frange(1,32), orig_index=False)
      self.assertAtomsEqual(at, sub)

   def test_read_nc_range_all(self):
      self.at.write('test.nc')
      at = Atoms('test.nc', range=[1,64])
      self.assertAtomsEqual(at, self.at)

   def test_read_nc_range_subset(self):
      self.at.write('test.nc')
      at = Atoms('test.nc', range=[1,32])
      sub = self.at.select(list=frange(1,32), orig_index=False)
      self.assertAtomsEqual(at, sub)

   def test_write_ext_string(self):
      es = Extendable_str()
      self.at.write('', estr=es, format='xyz')
      self.assertEqual(str(es), ''.join(self.xyz_ref))

   def test_cinoutput_query_frame(self):
      al = AtomsList([Atoms(n=n, lattice=fidentity(3)*n) for n in (0,1,2,3,4)])
      al.write('test.xyz')
      self.assertEqual([CInOutput('test.xyz', frame=n).n_atom for n in range(5) ], [0, 1, 2, 3, 4])

   def test_properties_parse_bug(self):
      cio = CInOutput()
      at = cio.read(str="""12
Lattice="3.679445 -0.181265 -0.219095 -0.140937 3.759568 0.175386 -0.223013 -0.011764 9.715448" Properties=species:S:1:pos:R:3:frac_pos:R:3:f:R:3 castep_run_time=40430.63 energy=-9949.25932697 virial="6.38880278817 4.17889707917 7.15229427215 4.17889707917 -0.740653778339 -2.86812348957 7.15229427215 -2.86812348957 -4.60133902559"
O       0.07892115      1.81915384      1.95631150      0.05184500      0.48697900      0.19373900     -0.08661000      1.65235000      0.57685000
O       0.03867388      1.98421766     -1.99187494      0.01775400      0.52796400     -0.21415200     -0.53478000     -0.09838000      1.11345000
O      -0.00808325     -0.07057365      4.35164490      0.02438200     -0.01619200      0.44875200      1.66233000      0.71894000      0.73628000
O      -0.07597036      0.06054846      0.43978436     -0.01735400      0.01540800      0.04459700      0.55106000     -0.33373000     -0.94889000
O       1.88194480      1.90132724      5.13029557      0.56414400      0.53459200      0.53112700     -0.50613000     -0.96190000      2.26224000
O       1.86840060     -0.05023179      2.83158150      0.52666200      0.01298000      0.30309400     -0.37249000      0.16802000     -0.89396000
O       1.83382247     -0.00862026     -2.78189090      0.48244900      0.02010500     -0.27582000      0.48621000      0.09282000     -1.55320000
O       1.90053677      1.80865713     -0.42277565      0.53347300      0.50667500     -0.04063200     -2.57147000     -0.07581000     -1.09444000
Ti      -0.09878726      1.99062810      0.05586401     -0.00682000      0.52914200     -0.00395600      2.77217000     -0.23594000      0.91766000
Ti       0.08835527      0.05867827      2.34945896      0.03940400      0.01826600      0.24238600     -0.43629000     -0.90315000      1.27819000
Ti       1.90982467      1.87755198     -2.42094136      0.52416000      0.52390600     -0.24682200     -0.51711000     -0.31319000     -1.29767000
Ti       1.88176829      0.00352974      4.80843526      0.54323600      0.02871600      0.50665900     -0.44688000      0.28996000     -1.09650000""")
      self.assertArrayAlmostEqual(at.pos,   [[0.07892115 ,    1.81915384 ,    1.95631150],
                                             [ 0.03867388,     1.98421766,    -1.99187494],
                                             [-0.00808325,    -0.07057365,     4.35164490],
                                             [-0.07597036,     0.06054846,     0.43978436],
                                             [ 1.88194480,     1.90132724,     5.13029557],
                                             [ 1.86840060,    -0.05023179,     2.83158150],
                                             [ 1.83382247,    -0.00862026,    -2.78189090],
                                             [ 1.90053677,     1.80865713,    -0.42277565],
                                             [-0.09878726,     1.99062810,     0.0558640],
                                             [ 0.08835527,     0.05867827,     2.3494589],
                                             [ 1.90982467,     1.87755198,    -2.4209413],
                                             [ 1.88176829,     0.00352974,     4.8084352]])




                                  
      
            

try:
   import netCDF4
   got_netcdf4 = True
except ImportError:
   got_netcdf4 = False

class TestPythonNetCDF(QuippyTestCase):
   def setUp(self):
      self.at = supercell(diamond(5.44,14), 2,2,2)
      self.al = AtomsList([ supercell(diamond(5.44+0.01*x,14),2,2,2) for x in range(5) ])
      self.al.write('test3.nc', netcdf4=False)

   def tearDown(self):
      if os.path.exists('test3.nc'): os.remove('test3.nc')
      if os.path.exists('dataset.nc'): os.remove('dataset.nc')

   def testpupynere_read(self):
      from quippy.pupynere import netcdf_file
      nc = netcdf_file('test3.nc', 'r')
      al = AtomsList(nc, format=quippy.netcdf.netcdf_file)
      self.assertEqual(list(self.al), list(al))
      nc.close()

   if got_netcdf4:
      
      def testnetcdf4_read(self):
         from netCDF4 import Dataset
         nc = Dataset('test3.nc','r')
         al = AtomsList(nc)
         for a, b in zip(self.al, al):
            self.assertAtomsEqual(a, b)
         nc.close()

      def testnetcdf4_write(self):
         from netCDF4 import Dataset
         nc = Dataset('dataset.nc','w')
         al2 = AtomsList(self.al)
         al2.write(nc)
         nc.close()
         al = AtomsList('dataset.nc')
         self.assertEqual(list(self.al), list(al))

class TestPuPyXYZ(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 2, 2, 2)
      self.at.add_property('log', False)
      self.at.params['real'] = 1.0
      self.at.params['int'] = 2
      self.at.params['int_a'] = [1,2,3]
      self.at.params['real_a'] = [1.0,2.0,3.0]
      self.at.params['int_a2'] = farray([1,2,3,4,5,6,7,8,9]).reshape(3,3)
      self.at.params['real_a2'] = farray([1.0,2,3,4,5,6,7,8,9]).reshape(3,3)
      self.at.params['log_param'] = True
      self.at.params['log_a'] = [True, True, False]
      self.al = AtomsList(self.at for x in range(5))
      self.al = AtomsList([ supercell(diamond(5.44+0.01*x,14),2,2,2) for x in range(5) ])
      for a in self.al:
         a.params.update(self.at.params)
         a.add_property('log', False)


   def tearDown(self):
      if os.path.exists('test.xyz'): os.remove('test.xyz')

   def testsinglexyz(self):
      self.at.write(PuPyXYZWriter('test.xyz'))
      at = Atoms(PuPyXYZReader('test.xyz'))
      self.assertAtomsEqual(at, self.at)

   def testmultixyz(self):
      self.al.write(PuPyXYZWriter('test.xyz'))
      al = AtomsList(PuPyXYZReader('test.xyz'), lazy=False)
      self.assertEqual(len(al), 5)
      self.assertEqual(len(self.al), len(al))
      self.assertEqual(list(self.al), list(al))

   def teststring(self):
      s = self.at.write('string')
      a = Atoms(s, format='string')
      self.assertEqual(a, self.at)

   def testmultistring(self):
      s = self.al.write('string')
      al = AtomsList(s, format='string')
      self.assertEqual(list(al), list(self.al))

   def testwrite_single_xyz_properties(self):
      self.at.write('test.xyz', properties=['species','pos'], format='pupyxyz')
      at = Atoms('test.xyz', format='pupyxyz')
      self.assertEqual(sorted(at.properties.keys()), sorted(['species', 'pos', 'Z']))
      self.assertArrayAlmostEqual(self.at.pos, at.pos)
      self.assertEqual(list(self.at.z), list(at.z))

   def testwrite_multi_xyz_properties(self):
      self.al.write('test.xyz', properties=['species','pos'], format='pupyxyz')
      al = AtomsList('test.xyz', format='pupyxyz')
      for at,at_ref in zip(al, self.al):
         self.assertEqual(sorted(at.properties.keys()), sorted(['species', 'pos', 'Z']))
         self.assertArrayAlmostEqual(at.pos, at_ref.pos)
         self.assertEqual(list(at.z), list(at_ref.z))

if __name__ == '__main__':
   unittest.main()
