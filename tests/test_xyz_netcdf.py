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
from quippy.io import *
from quippy.netcdf import *
import unittest, itertools, sys, quippy
import numpy as np
from quippytest import *
import os
import glob

class TestCInOutput(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44,14), 2,2,2)
      self.at.set_cutoff(3.0)
      self.at.add_property('log', False)
      self.at.log[1] = True
      self.at.params['real'] = 1.0
      self.at.params['int'] = 2
      self.at.params['neg_int'] = -3
      self.at.params['bad_neg'] = '3-4'
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

      self.xyz_ref =  ['64\n', 'real=1.00000000 int=2 neg_int=-3 bad_neg=3-4 int_a="1       2       3" real_a="1.00000000      2.00000000      3.00000000" int_a2="1        4        7        2        5        8        3        6        9" real_a2="1.00000000       4.00000000       7.00000000       2.00000000       5.00000000       8.00000000       3.00000000       6.00000000       9.00000000" log_param=T log_a="T T F" string=string string2="string with spaces" cutoff=3.00000000 nneightol=1.20000000 pbc="T T T" Lattice="10.88000000       0.00000000       0.00000000       0.00000000      10.88000000       0.00000000       0.00000000       0.00000000      10.88000000" Properties=species:S:1:pos:R:3:Z:I:1:log:L:1\n',
                       'Si              0.00000000      0.00000000      0.00000000      14    T\n',
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
      if os.path.exists('empty.xyz'): os.remove('empty.xyz')
      for t in glob.glob('test*.xyz*'):
         os.remove(t)

   def testsinglexyz(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz')
      self.assertEqual(self.at, at)
      self.assertEqual(self.xyz_ref, open('test.xyz', 'r').readlines())

   def testsinglexyzprefix(self):
      self.at.write('test.xyz', prefix='PREFIX')
      lines = open('test.xyz').readlines()
      self.assert_(all([line[:len('PREFIX')] == 'PREFIX' for line in lines]))
      lines_without_prefix = [line[len('PREFIX '):] for line in lines]
      at = Atoms(''.join(lines_without_prefix), format='string')
      self.assertEqual(self.at, at)

   def testxyz_negative_integer(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz')
      self.assertEqual(type(at.neg_int), type(1))

   def testxyz_bad_negative_integer(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz')
      self.assertEqual(type(at.bad_neg), type(''))

   if 'netcdf' in available_modules:
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

   def testmultixyz_prefix_write(self):
      self.al.write('test.xyz', prefix='PREFIX')
      lines = open('test.xyz').readlines()
      self.assert_(all([line[:len('PREFIX')] == 'PREFIX' for line in lines]))

   def testmultixyz_prefix_read(self):
      self.al.write('test.xyz', prefix='PREFIX')
      al = AtomsList('test.xyz')
      self.assert_(al[0].xyz_prefix == 'PREFIX')
      for at in al:
         del at.params['xyz_prefix']
      self.assertEqual(list(self.al), list(al))

   if 'netcdf' in available_modules:
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

   if 'netcdf' in quippy.available_modules:

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
      self.assertRaises(EOFError, cio.read, frame=5)
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
      from quippy.structures import quartz_params
      aq1 = alpha_quartz(**quartz_params['ASAP_JRK'])
      aq1.write('quartz.xyz')
      aq2 = Atoms('quartz.xyz')
      self.assertEqual(aq1, aq2)

   if 'netcdf' in available_modules:
      def test_non_orthorhombic_nc(self):
         from quippy.structures import quartz_params
         aq1 = alpha_quartz(**quartz_params['ASAP_JRK'])
         aq1.write('quartz.nc', netcdf4=False)
         aq2 = Atoms('quartz.nc')
         self.assertEqual(aq1, aq2)

   def test_read_string(self):
      s = ''.join(self.xyz_ref)
      cio = CInOutput()
      at = cio.read(str=s)
      self.assertEqual(at, self.at)

   def test_read_ext_string(self):
      es = Extendable_str(''.join(self.xyz_ref))
      cio = CInOutput()
      at = cio.read(estr=es)
      self.assertEqual(at, self.at)

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
      self.assertEqual(at, self.at)

   def test_read_xyz_range_subset(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', range=[2,32])
      sub = self.at.select(list=frange(2,32), orig_index=False)
      self.assertEqual(at, sub)

   def test_read_xyz_range_empty(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', range='empty')
      sub = self.at.select(list=[], orig_index=False)
      self.assertEqual(at, sub)

   if 'netcdf' in available_modules:
      def test_read_nc_range_all(self):
         self.at.write('test.nc')
         at = Atoms('test.nc', range=[1,64])
         self.assertEqual(at, self.at)

      def test_read_nc_range_subset(self):
         self.at.write('test.nc')
         at = Atoms('test.nc', range=[2,32])
         sub = self.at.select(list=frange(2,32), orig_index=False)
         self.assertEqual(at, sub)

      def test_read_nc_range_empty(self):
         self.at.write('test.nc')
         at = Atoms('test.nc', range='empty')
         sub = self.at.select(list=[], orig_index=False)
         self.assertEqual(at, sub)

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
      self.assertArrayAlmostEqual(at.pos.T,   [[0.07892115 ,    1.81915384 ,    1.95631150],
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


   def test_long_properties_string_bug(self):
      # caused segmentation fault due to buffer overrun in xyz.c when splitting properties string into fields
      # fixed by ensuring we don't run off end of fields array, and increasing MAX_FIELD_COUNT
      cio = CInOutput()
      at = cio.read(str="""57
OrigWidth=526.11583369 OrigHeight=170.89163941 YoungsModulus=129.04197204 PoissonRatio_yx=0.14345535 PoissonRatio_yz=0.43228896 G=3.07200000 CrackPosx=-115.66261165 CrackPosy=-1.23372306 OrigCrackPos=-120.55791685 Time=1072.00000000 Temp=12.08187913 LastStateChangeTime=0.00000000 LastMDIntervalTime=0.00000000 LastCalcConnectTime=1071.00000000 State=MD_LOADING LastPrintTime=1065.00000000 LastCheckpointTime=1050.00000000 Lattice="25.200000 0.000000 0.000000 0.000000 25.200000 0.000000 0.000000 0.000000 19.152000" Properties=species:S:1:pos:R:3:Z:I:1:travel:I:3:map_shift:I:3:hybrid:I:1:hybrid_mark:I:1:changed_nn:I:1:move_mask:I:1:nn:I:1:old_nn:I:1:md_old_changed_nn:I:1:edge_mask:I:1:crack_surface:L:1:crack_front:L:1:load:R:3:local_e:R:1:mass:R:1:damp_mask:I:1:thermostat_region:I:1:avg_ke:R:1:velo:R:3:acc:R:3:avgpos:R:3:oldpos:R:3:force:R:3:qm_force:R:3:mm_force:R:3:weight_region1:R:1:modified_hybrid_mark:I:1:orig_index:I:1:index:I:1:shift:I:3:termindex:I:1:rescale:R:1:cluster_ident:S:1
Si       8.71395610    -11.96756847     -1.23127991      14      -5       0       0       0       0       0       1       1       0       1       3       3       0       1    F    F      0.00000000     -0.05631387      0.00000000     -2.13989447   2910.85773629       1       1      0.00779406      0.00146777      0.00039089     -0.00096292      0.00007544     -0.00002422     -0.00000000   -116.35243582    -11.86333217     -1.17083136   -116.34701002    -11.74808703     -0.94746363      0.21960161     -0.07048883     -0.00000286      0.00000000      0.00000000      0.00000000      0.21960161     -0.07048883     -0.00000286      1.00000000       1    6917    6917       0       0       0       0      1.00000000 h_active
Si       4.30478574     10.87644762      0.17727867      14      -6      -1       0       0       0       0       0       5       0       1       4       4       0       1    F    F      0.00000000     -0.06746728      0.00000000     -4.31534112   2910.85773629       1       1      0.01013096      0.00069680      0.00021642      0.00200411     -0.00000167     -0.00001051     -0.00008246   -746.89000076    -14.21863402      0.06912598   -747.01781783    -14.25631961     -0.08767505     -0.00485419     -0.03060039     -0.24001896      0.00000000      0.00000000      0.00000000     -0.00485419     -0.03060039     -0.24001896      0.00000000       5    6558    6558       0       0       0       0      1.00000000 h_outer_l
Si       8.24970338      6.64448083      0.20086606      14      -6      -1       0       0       0       0       0       5       0       1       4       4       0       1    F    F      0.00000000     -0.08623099      0.00000000     -4.31444149   2910.85773629       1       1      0.02900658      0.00051252     -0.00353710      0.00069644      0.00000302      0.00000109     -0.00009738   -742.95634925    -18.20859743      0.15292815   -742.99480221    -17.80811792      0.06466226      0.00877686      0.00318345     -0.28345732      0.00000000      0.00000000      0.00000000      0.00877686      0.00318345     -0.28345732      0.00000000       5    6734    6734       0       0       0       0      1.00000000 h_outer_l
Si       6.22714153      9.42009846      2.86410960      14      -5      -1       0       0       0       0       0       5       0       1       4       4       0       1    F    F      0.00000000     -0.07683360      0.00000000     -4.31059612   2910.85773629       1       1      0.02632547      0.00112384      0.00305052      0.00053796      0.00000095      0.00001112     -0.00017437   -118.88297454    -15.80754992      2.75417035   -119.00689027    -16.10494844      2.67744696      0.00276956      0.03236137     -0.50757693      0.00000000      0.00000000      0.00000000      0.00276956      0.03236137     -0.50757693      0.00000000       5    6735    6735       0       0       0       0      1.00000000 h_outer_l
Si       8.76887118      9.24036128     -0.83260069      14      -5      -1       0       0       0       0       0       2       0       1       3       3       0       1    F    F      0.00000000     -0.07625862      0.00000000     -2.15281270   2910.85773629       1       1      0.00305664      0.00045064      0.00018728      0.00119890      0.00000450      0.00000624      0.00003874   -116.31665248    -15.82031222     -0.87896734   -116.31869102    -15.75430588     -1.03291962      0.01310705      0.01815300      0.11277990      0.00000000      0.00000000      0.00000000      0.01310705      0.01815300      0.11277990      0.00000000       2    6913    6913       0       0       0       0      1.00000000 h_buffer
Si       3.96048525    -10.42335741     -0.06137142      14      -6       0       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.04758470      0.00000000     -4.21495814   2910.85773629       1       1      0.00442838      0.00067761      0.00145321     -0.00052192      0.00004989      0.00001432     -0.00012611   -747.16026074    -10.33259345     -0.07169561   -747.12659791    -10.28616791     -0.09415054      0.14522785      0.04168196     -0.36708851      0.00000000      0.00000000      0.00000000      0.14522785      0.04168196     -0.36708851      0.00000000       5    6562    6562       0       0       0       0      1.00000000 h_outer_l
Si       8.22103073     11.86639581     -0.08813685      14      -6      -1       0       0       0       0       1       2       0       1       4       4       0       1    F    F      0.00000000     -0.06628625      0.00000000     -4.30160038   2910.85773629       1       1      0.03389340      0.00055341      0.00427216     -0.00149836      0.00000875      0.00001817     -0.00012102   -742.98676013    -13.43112253     -0.02146731   -743.08033073    -13.82577056      0.17519905      0.02548120      0.05289868     -0.35226923      0.00000000      0.00000000      0.00000000      0.02548120      0.05289868     -0.35226923      0.00000000       2    6738    6738       0       0       0       0      1.00000000 h_buffer
Si       6.08949575    -12.18701655      2.54510076      14      -5       0       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.05691993      0.00000000     -4.27704940   2910.85773629       1       1      0.00572211      0.00074019      0.00123365     -0.00241284     -0.00001656      0.00002967     -0.00019042   -118.99769996    -12.11018795      2.56356325   -119.10437854    -12.14059705      2.70180516     -0.04819948      0.08637297     -0.55429896      0.00000000      0.00000000      0.00000000     -0.04819948      0.08637297     -0.55429896      0.00000000       5    6739    6739       0       0       0       0      1.00000000 h_outer_l
Si       6.29748788     11.34259538      1.32877748      14      -5      -1       0       0       0       0       1       2       0       1       4       4       0       1    F    F      0.00000000     -0.06687677      0.00000000     -4.30003543   2910.85773629       1       1      0.01603555      0.00142988      0.00243740      0.00027116      0.00000572      0.00002026      0.00014613   -118.82828754    -13.86332113      1.32883694   -119.01712234    -14.11458198      1.36228272      0.01665546      0.05897411      0.42537476      0.00000000      0.00000000      0.00000000      0.01665546      0.05897411      0.42537476      0.00000000       2    6740    6740       0       0       0       0      1.00000000 h_buffer
Si       7.84906673    -10.44177452     -0.31736195      14      -6       0       0       0       0       0       1       2       0       1       4       4       0       1    F    F      0.00000000     -0.04634150      0.00000000     -3.30618439   2910.85773629       1       1      0.03997914     -0.00197616     -0.00742133      0.00413833     -0.00011553     -0.00038404      0.00017345   -743.28897556    -10.17038604     -0.25382676   -743.17088574     -9.73599759      0.11691124     -0.33628823     -1.11789726      0.50489062      0.00000000      0.00000000      0.00000000     -0.33628823     -1.11789726      0.50489062      0.00000000       2    6742    6742       0       0       0       0      1.00000000 h_buffer
Si       5.98713586     -9.92119355      1.08560365      14      -5       0       0       0       0       0       1       2       0       1       4       4       0       1    F    F      0.00000000     -0.04696310      0.00000000     -4.66963964   2910.85773629       1       1      0.02383987     -0.00067738      0.00303798     -0.00142843     -0.00030403     -0.00025233      0.00021366   -119.07316224     -9.99043490      1.15955172   -119.08823909    -10.22269395      1.29900355     -0.88499322     -0.73449236      0.62193190      0.00000000      0.00000000      0.00000000     -0.88499322     -0.73449236      0.62193190      0.00000000       2    6744    6744       0       0       0       0      1.00000000 h_buffer
Si      11.22363885      9.45966354     -0.83917718      14      -5      -1       0       0       0       0       1       5       0       1       3       3       0       1    F    F      0.00000000     -0.07510866      0.00000000     -2.14978554   2910.85773629       1       1      0.00296022      0.00049440      0.00049279      0.00137614     -0.00000107      0.00004616      0.00006302   -113.88983290    -15.60761646     -0.87725268   -113.92986303    -15.57193400     -0.97714986     -0.00312593      0.13435155      0.18344644      0.00000000      0.00000000      0.00000000     -0.00312593      0.13435155      0.18344644      0.00000000       5    7093    7093       0       0       0       0      1.00000000 h_outer_l
Si       9.95045330     11.19797357      1.75273442      14      -5      -1       0       0       0       0       1       2       0       1       4       4       0       1    F    F      0.00000000     -0.06569573      0.00000000     -4.30020605   2910.85773629       1       1      0.00750617      0.00000564     -0.00113389      0.00155563     -0.00000950      0.00001654      0.00019337   -115.14265935    -13.80948975      1.71327816   -115.19828960    -13.64629804      1.68443551     -0.02766016      0.04815927      0.56286846      0.00000000      0.00000000      0.00000000     -0.02766016      0.04815927      0.56286846      0.00000000       2    6920    6920       0       0       0       0      1.00000000 h_buffer
Si       9.85841239      9.33874157      3.16264090      14      -5      -1       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.07568364      0.00000000     -4.31037030   2910.85773629       1       1      0.01091468     -0.00037563     -0.00035816      0.00089415     -0.00000303      0.00001118     -0.00018709   -115.20222823    -15.72218466      3.02629839   -115.18652178    -15.71284006      2.87408467     -0.00883202      0.03255002     -0.54459956      0.00000000      0.00000000      0.00000000     -0.00883202      0.03255002     -0.54459956      0.00000000       5    6915    6915       0       0       0       0      1.00000000 h_outer_l
Si       9.81327870    -11.84422120      3.05122402      14      -5       0       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.05570781      0.00000000     -4.24984898   2910.85773629       1       1      0.00683173      0.00041979     -0.00130011      0.00177176      0.00000736     -0.00005516     -0.00011338   -115.27721752    -11.71372308      2.95517169   -115.30445907    -11.60915538      2.90894955      0.02141517     -0.16055717     -0.33003037      0.00000000      0.00000000      0.00000000      0.02141517     -0.16055717     -0.33003037      0.00000000       5    6919    6919       0       0       0       0      1.00000000 h_outer_l
Si      11.66890056     11.41384494      0.53307967      14      -6      -1       0       0       0       0       1       2       0       1       4       4       0       1    F    F      0.00000000     -0.06510521      0.00000000     -4.28624802   2910.85773629       1       1      0.01969188     -0.00037557     -0.00077550      0.00337366     -0.00002527     -0.00003538     -0.00011090   -739.53002452    -13.62438837      0.37241695   -739.51338756    -13.52095976      0.18848401     -0.07356122     -0.10297195     -0.32280376      0.00000000      0.00000000      0.00000000     -0.07356122     -0.10297195     -0.32280376      0.00000000       2    6918    6918       0       0       0       0      1.00000000 h_buffer
Si      11.38645753    -11.82752608     -0.74100243      14      -5       0       0       0       0       0       1       2       0       1       3       3       0       1    F    F      0.00000000     -0.05510176      0.00000000     -2.12964481   2910.85773629       1       1      0.02069508      0.00042249     -0.00248485      0.00147799     -0.00005302     -0.00004027     -0.00006819   -113.77536250    -11.65840519     -0.84966902   -113.98200015    -11.47149100     -0.92516973     -0.15432238     -0.11720947     -0.19848409      0.00000000      0.00000000      0.00000000     -0.15432238     -0.11720947     -0.19848409      0.00000000       2    7097    7097       0       0       0       0      1.00000000 h_buffer
Si       5.70022470     -8.52232737      2.77242842      14      -5       0       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.03700627      0.00000000     -4.62469133   2910.85773629       1       1      0.08724770     -0.00142425     -0.01143400      0.00162045     -0.00034840     -0.00028492      0.00005331   -119.29740956     -8.22466810      2.66107840   -119.30175934     -8.38795899      2.74475577     -1.01413346     -0.82935654      0.15517584      0.00000000      0.00000000      0.00000000     -1.01413346     -0.82935654      0.15517584      0.00000000       5    6743    6743       0       0       0       0      1.00000000 h_outer_l
Si       7.85418145     -7.84729771      1.34453271      14      -5       0       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.03636913      0.00000000     -1.74747897   2910.85773629       1       1      1.22206433     -0.01957119     -0.00310448      0.02707448      0.00060825      0.00000821     -0.00051521   -116.64193052     -7.71616759      0.23121038   -116.50148519     -7.65868548     -1.05520941      1.77051542      0.02390028     -1.49969834      1.77051542      0.02390028     -1.49969834      0.50897018      0.61143311     -5.48834847      0.00000000       5    6921    6921       0       0       0       0      1.00000000 h_outer_l
Si      11.41712809     -9.46096197     -0.34410070      14      -6       0       0       0       0       0       1       2       0       1       3       3       0       1    F    F      0.00000000     -0.04509830      0.00000000     -4.28332711   2910.85773629       1       1      0.06593045      0.00231900     -0.00308341     -0.00456046      0.00008667     -0.00018885      0.00009085   -739.77412910     -9.32805683     -0.08357541   -739.66050457     -9.34677559      0.14993922      0.25227775     -0.54970697      0.26443978      0.27247178     -0.32464584      0.29180001      0.25227775     -0.54970697      0.26443978      0.00000000       2    6922    6922       0       0       0       0      1.00000000 h_buffer
Si       9.98544849     -9.52722500      2.14731514      14      -5       0       0       0       0       0       1       5       0       1       5       5       0       1    F    F      0.00000000     -0.04571990      0.00000000     -4.15344423   2910.85773629       1       1      0.26521356      0.01307194     -0.00886248      0.01870463      0.00039493     -0.00027661      0.00037549   -115.30331152     -9.37823805      1.72966257   -115.32908631     -9.54570985      1.68955998      1.14958177     -0.80517544      1.09299622      0.26170265      0.05871369      0.20478374      1.14958177     -0.80517544      1.09299622      0.00000000       5    6924    6924       0       0       0       0      1.00000000 h_outer_l
Si      11.48771523     -7.57364361     -1.61081637      14      -5       0       0       0       0       0       1       5       0       1       1       1       0       1    F    F      0.00000000     -0.03509485      0.00000000     -1.60800483   2910.85773629       1       1      0.22686516      0.01378965     -0.00167090     -0.00062736     -0.00054110     -0.00015493      0.00011478   -113.97542357     -7.45774315     -1.46795582   -114.11043672     -7.30599968     -1.15447153     -1.57506038     -0.45098774      0.33411710     -1.57506038     -0.45098774      0.33411710      0.15999617      0.16114721     -0.08548830      0.00000000       5    7101    7101       0       0       0       0      1.00000000 h_outer_l
Si     -11.03033501     11.75219435      1.29130564      14      -4      -1       0       0       0       0       1       5       0       1       4       4       0       1    F    F      0.00000000     -0.06503454      0.00000000     -4.31586931   2910.85773629       1       1      0.03759277      0.00335830     -0.00051664     -0.00095021      0.00001096      0.00003063      0.00022685   -111.12841423    -13.35616121      1.39539547   -111.49993821    -13.41710190      1.61425010      0.03191026      0.08916064      0.66031378      0.00000000      0.00000000      0.00000000      0.03191026      0.08916064      0.66031378      0.00000000       5    7100    7100       0       0       0       0      1.00000000 h_outer_l
Si       4.30141000    -12.50936964     -1.10699120      14      -5       0       0       0       0       0       0       0       0       1       2       2       0       1    F    F      0.00000000     -0.05752599      0.00000000     -2.15617401   2910.85773629       1       1      0.01278483      0.00150931     -0.00089190      0.00155151      0.00000126      0.00004012      0.00002524   -120.81685599    -12.35376198     -1.19034796   -120.99904636    -12.26412407     -1.39878588      0.00366528      0.11679367      0.07346097      0.00000000      0.00000000      0.00000000      0.00366528      0.11679367      0.07346097      0.00000000       2    6737    6737       0       0       0       0      1.00000000 n_cut_bond
Si       6.13777577     -5.92233824      1.95214599      14      -5       0       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.02704943      0.00000000     -4.89075486   2910.85773629       1       1      0.08187159      0.00556532      0.00174036     -0.00271442     -0.00013072      0.00004224      0.00008725   -119.12253238     -6.04156107      1.82922149   -119.46941711     -6.57344474      1.32673569     -0.38051410      0.12294128      0.25397784     -0.38051410      0.12294128      0.25397784     -1.39817712      1.62585343     -0.47993452      0.00000000       2    6748    6748       0       0       0       0      1.00000000 clash
Si       7.42133286     -7.73647704      4.14880166      14      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    T    T      0.00000000     -0.03636913      0.00000000     -3.74571525   2910.85773629       1       1      0.01286734     -0.00188239      0.00363703      0.00141049      0.00002775      0.00015554      0.00007423   -117.54659085     -7.75024158      4.10533238   -117.44667891     -7.90024304      4.16909004      0.08078982      0.45274337      0.21608694      0.08078982      0.45274337      0.21608694      0.39881032     -0.13025234      3.17559729      0.00000000       2   31581   31581       0       0       0       0      1.00000000 clash
Si       7.99530988     -5.11845503      0.53795829      14      -6       0       0       0       0       0       1       0       0       1       3       3       0       1    F    F      0.00000000     -0.02639675      0.00000000     -3.71365327   2910.85773629       1       1      0.11187492      0.00190241      0.00794512      0.00834210     -0.00010123      0.00007942      0.00042385   -743.37800119     -5.37407621      0.38696704   -743.71284391     -5.80815305      0.13575982     -0.29467038      0.23119306      1.23377531     -0.29467038      0.23119306      1.23377531     -0.34948146      0.16150864      1.22087116      0.00000000       2    6746    6746       0       0       0       0      1.00000000 clash
Si       9.53459971     -7.25580075      3.01507012      14      -5       0       0       0       0       0       1       0       0       1       5       5       0       1    F    F      0.00000000     -0.03573199      0.00000000     -4.61910967   2910.85773629       1       1      0.04010370      0.00014906     -0.00069275      0.00338637     -0.00017563     -0.00035196     -0.00017357   -115.54721317     -7.32615814      2.85418130   -115.51964973     -7.48221089      2.90549397     -0.51124018     -1.02450749     -0.50523828     -0.51124018     -1.02450749     -0.50523828      0.75387553     -0.93058256     -0.62857407      0.00000000       2    6923    6923       0       0       0       0      1.00000000 clash
H        4.21816809      9.74753283     -0.63963980       1      -5      -1       0       0       0       0       0       0       0       1       2       2       0       1    F    F      0.00000000     -0.07740858      0.00000000     -2.15723222   2910.85773629       1       1      0.01267813      0.00020428      0.00199591      0.00215515      0.00000208      0.00000895      0.00004048   -120.86600975    -16.04036967     -1.18291903   -120.84066948    -16.25311200     -1.38924439      0.00604728      0.02605144      0.11783933      0.00000000      0.00000000      0.00000000      0.00604728      0.02605144      0.11783933      0.00000000       0    6733    6733       0       0       0       2      0.64414414 term
H        3.05661919     10.73654065      0.96954628       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.06805780      0.00000000     -4.31379667   2910.85773629       1       1      0.00758973      0.00118633      0.00029801      0.00167791      0.00000249     -0.00001844      0.00012965   -122.71735753    -14.43222576      1.33420216   -122.86729374    -14.43934968      1.22760771      0.00724924     -0.05366642      0.37740533      0.00000000      0.00000000      0.00000000      0.00724924     -0.05366642      0.37740533      0.00000000       0    6560    6560       0       0       0       2      0.64414414 term
H        6.68443388      6.91519875      1.16081358       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.08679043      0.00000000     -4.31253572   2910.85773629       1       1      0.01756779     -0.00166669      0.00072773      0.00286961     -0.00000072     -0.00000079      0.00016075   -119.12685112    -18.01914126      1.56585328   -118.88270456    -18.05336279      1.32013417     -0.00208781     -0.00230899      0.46793137      0.00000000      0.00000000      0.00000000     -0.00208781     -0.00230899      0.46793137      0.00000000       0    6736    6736       0       0       0       3      0.64414414 term
H        9.42117765      5.78856427     -0.27646085       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.09620337      0.00000000     -2.15686210   2910.85773629       1       1      0.14933162      0.00738870      0.00053052      0.00426835      0.00000119     -0.00000037      0.00003388   -115.42429157    -19.75303214     -0.75709853   -116.35039074    -19.79223431     -1.19350032      0.00347482     -0.00107448      0.09863139      0.00000000      0.00000000      0.00000000      0.00347482     -0.00107448      0.09863139      0.00000000       0    6909    6909       0       0       0       3      0.64414414 term
H        9.32232285      7.02744195      1.42718730       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.08567155      0.00000000     -4.31484916   2910.85773629       1       1      0.01671089     -0.00022165     -0.00085511      0.00340673      0.00000256      0.00000160      0.00017180   -115.14739341    -17.77663034      1.94128628   -115.12545534    -17.72249855      1.59959992      0.00745367      0.00466664      0.50009371      0.00000000      0.00000000      0.00000000      0.00745367      0.00466664      0.50009371      0.00000000       0    6916    6916       0       0       0       3      0.64414414 term
H        9.43227633      5.68148488     -0.29443867       1      -5      -1       0       0       0       0       0       0       0       1       3       3       0       1    F    F      0.00000000     -0.09511557      0.00000000     -2.15491707   2910.85773629       1       1      0.14656804     -0.00740069     -0.00180744      0.00355504      0.00000218      0.00002415      0.00003508   -114.57747503    -19.79736790     -0.74610216   -113.66532058    -19.61901192     -1.11134152      0.00633677      0.07029401      0.10211938      0.00000000      0.00000000      0.00000000      0.00633677      0.07029401      0.10211938      0.00000000       0    7089    7089       0       0       0       3      0.64414414 term
H        5.96469340      7.90291855      2.10854441       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.08679043      0.00000000     -4.31253572   2910.85773629       1       1      0.01756779     -0.00166669      0.00072773      0.00286961     -0.00000072     -0.00000079      0.00016075   -119.12685112    -18.01914126      1.56585328   -118.88270456    -18.05336279      1.32013417     -0.00208781     -0.00230899      0.46793137      0.00000000      0.00000000      0.00000000     -0.00208781     -0.00230899      0.46793137      0.00000000       0    6736    6736       0       0       0       4      0.64414414 term
H        4.91411539      8.93437353      3.63168261       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       0    F    F      0.00000000     -0.07740858      0.00000000     -4.31560604   2910.85773629       1       1      0.00368257      0.00061326     -0.00105576      0.00221356     -0.00000178     -0.00001457      0.00008608   -120.87863570    -16.34974938      4.00388594   -120.95481014    -16.21785219      4.00461151     -0.00519205     -0.04239856      0.25055600      0.00000000      0.00000000      0.00000000     -0.00519205     -0.04239856      0.25055600      0.00000000       0   31393   31393       0       0       0       4      0.64414414 term
H        7.42876151      9.21229359      3.80529495       1      -5      -1       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.07625862      0.00000000     -4.31757643   2910.85773629       1       1      0.00599390      0.00101847     -0.00076784      0.00314533     -0.00000304     -0.00001813      0.00010837   -117.02922450    -15.94794557      4.23164295   -117.15491649    -15.88848484      4.14710696     -0.00884007     -0.05276345      0.31545655      0.00000000      0.00000000      0.00000000     -0.00884007     -0.05276345      0.31545655      0.00000000       0   31573   31573       0       0       0       4      0.64414414 term
H        2.75002990    -10.68524444      0.83904276       1      -5       0       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.04820630      0.00000000     -4.30791785   2910.85773629       1       1      0.00648023      0.00088743     -0.00140633      0.00136876      0.00001198      0.00003041      0.00009493   -122.95658266    -10.66552419      1.28236479   -123.01031487    -10.51700876      1.17209188      0.03486098      0.08850727      0.27632895      0.00000000      0.00000000      0.00000000      0.03486098      0.08850727      0.27632895      0.00000000       0    6564    6564       0       0       0       6      0.64414414 term
H        4.07390408     -9.18439222     -0.87526970       1      -5       0       0       0       0       0       1       0       0       1       2       2       0       1    F    F      0.00000000     -0.03764340      0.00000000     -1.82898998   2910.85773629       1       1      0.00880836      0.00216046      0.00011270      0.00056599      0.00003324     -0.00007319     -0.00002621   -120.96268987     -8.42716822     -1.36594846   -121.09312291     -8.31741878     -1.44777236      0.09675248     -0.21303812     -0.07628997      0.00000000      0.00000000      0.00000000      0.09675248     -0.21303812     -0.07628997      0.00000000       0    6741    6741       0       0       0       6      0.64414414 term
H        4.90118276    -12.38071205      3.42260272       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.05752599      0.00000000     -4.31500239   2910.85773629       1       1      0.00714379      0.00163463     -0.00031498      0.00112650      0.00001563      0.00001222      0.00011408   -120.87881476    -12.35721595      3.90261228   -121.08844973    -12.29690211      3.98549987      0.04550978      0.03556356      0.33205761      0.00000000      0.00000000      0.00000000      0.04550978      0.03556356      0.33205761      0.00000000       0   31397   31397       0       0       0       8      0.64414414 term
H        7.18469979    -12.02352660      3.69117092       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.05631387      0.00000000     -4.31244222   2910.85773629       1       1      0.00317310     -0.00025195      0.00006627      0.00220801     -0.00001190      0.00002159      0.00013694   -117.26272217    -11.84987875      4.24542500   -117.26036194    -11.88914865      4.17816773     -0.03463140      0.06285148      0.39861459      0.00000000      0.00000000      0.00000000     -0.03463140      0.06285148      0.39861459      0.00000000       0   31577   31577       0       0       0       8      0.64414414 term
H       11.54245515      7.76681718     -0.17975758       1      -6      -1       0       0       0       0       0       0       0       1       5       5       0       1    F    F      0.00000000     -0.08511212      0.00000000     -4.31422758   2910.85773629       1       1      0.03235268     -0.00033850     -0.00379030      0.00003999      0.00000135     -0.00001131     -0.00011072   -739.46969709    -18.00482853      0.16627425   -739.43172609    -17.57916740      0.14593390      0.00393066     -0.03292344     -0.32229356      0.00000000      0.00000000      0.00000000      0.00393066     -0.03292344     -0.32229356      0.00000000       0    6914    6914       0       0       0      12      0.64414414 term
H        9.89479138      7.98621043      2.48115224       1      -5      -1       0       0       0       0       0       0       0       1       4       4       0       1    F    F      0.00000000     -0.08567155      0.00000000     -4.31484916   2910.85773629       1       1      0.01671089     -0.00022165     -0.00085511      0.00340673      0.00000256      0.00000160      0.00017180   -115.14739341    -17.77663034      1.94128628   -115.12545534    -17.72249855      1.59959992      0.00745367      0.00466664      0.50009371      0.00000000      0.00000000      0.00000000      0.00745367      0.00466664      0.50009371      0.00000000       0    6916    6916       0       0       0      14      0.64414414 term
H        8.72097053      9.18334227      3.91152906       1      -5      -1       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.07625862      0.00000000     -4.31757643   2910.85773629       1       1      0.00599390      0.00101847     -0.00076784      0.00314533     -0.00000304     -0.00001813      0.00010837   -117.02922450    -15.94794557      4.23164295   -117.15491649    -15.88848484      4.14710696     -0.00884007     -0.05276345      0.31545655      0.00000000      0.00000000      0.00000000     -0.00884007     -0.05276345      0.31545655      0.00000000       0   31573   31573       0       0       0      14      0.64414414 term
H       10.98757946      9.28317328      4.09082914       1      -5      -1       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.07510866      0.00000000     -4.31999565   2910.85773629       1       1      0.01655784     -0.00130542     -0.00111060      0.00440899     -0.00001054     -0.00000413      0.00013554   -113.41330756    -15.74314929      4.45006850   -113.29170666    -15.57763185      4.24508799     -0.03068543     -0.01203417      0.39455141      0.00000000      0.00000000      0.00000000     -0.03068543     -0.01203417      0.39455141      0.00000000       0   31753   31753       0       0       0      14      0.64414414 term
H        8.50982978    -11.90154086      3.87127785       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.05631387      0.00000000     -4.31244222   2910.85773629       1       1      0.00317310     -0.00025195      0.00006627      0.00220801     -0.00001190      0.00002159      0.00013694   -117.26272217    -11.84987875      4.24542500   -117.26036194    -11.88914865      4.17816773     -0.03463140      0.06285148      0.39861459      0.00000000      0.00000000      0.00000000     -0.03463140      0.06285148      0.39861459      0.00000000       0   31577   31577       0       0       0      15      0.64414414 term
H       10.95736777    -11.74164293      3.97495123       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.05510176      0.00000000     -4.31755133   2910.85773629       1       1      0.01019699     -0.00086051      0.00005667      0.00326526      0.00003987      0.00003291      0.00012061   -113.46284069    -11.55376750      4.36109899   -113.40730988    -11.46777510      4.27883965      0.11606251      0.09580697      0.35108022      0.11050532      0.05435168      0.23013977      0.11606251      0.09580697      0.35108022      0.00000000       0   31757   31757       0       0       0      15      0.64414414 term
H        4.51322901     -8.47258601      3.52349225       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.03764340      0.00000000     -3.85005193   2910.85773629       1       1      0.00790893     -0.00014837     -0.00312480      0.00115689      0.00004979     -0.00001542      0.00001694   -121.17166388     -8.35004364      3.89734516   -121.29797249     -8.42534374      3.98221952      0.14494210     -0.04489052      0.04931330      0.00000000      0.00000000      0.00000000      0.14494210     -0.04489052      0.04931330      0.00000000       0   31401   31401       0       0       0      18      0.64414414 term
H      -10.09474754     12.07649194      0.46826709       1      -5      -1       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.06503454      0.00000000     -4.29174547   2910.85773629       1       1      0.03279935     -0.00031782      0.00382059     -0.00145587      0.00001705      0.00005455     -0.00016630   -735.59858429    -13.01281504      0.08061985   -735.59315149    -13.39547878      0.34375724      0.04962752      0.15878447     -0.48407994      0.00000000      0.00000000      0.00000000      0.04962752      0.15878447     -0.48407994      0.00000000       0    7098    7098       0       0       0      23      0.64414414 term
H      -11.48968957    -12.39450342      2.58788014       1      -4       0       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.05502922      0.00000000     -4.30005828   2910.85773629       1       1      0.04128402     -0.00024373     -0.00279770      0.00420177     -0.00004893      0.00008236     -0.00020742   -111.62990037    -11.56561926      3.08739028   -111.53444061    -11.34508945      2.84248115     -0.14243203      0.23975081     -0.60375602     -0.06652223      0.23083075     -0.42755097     -0.14243203      0.23975081     -0.60375602      0.00000000       0    7099    7099       0       0       0      23      0.64414414 term
H        4.79902731     -5.97365752      0.79717935       1      -6       0       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.02770211      0.00000000     -3.94410216   2910.85773629       1       1      0.04951243      0.00388260      0.00214916      0.00210994      0.00015937     -0.00023196      0.00004122   -747.20835279     -6.16061732      0.10011223   -747.46307773     -6.41595816     -0.07360451      0.46391121     -0.67520143      0.11999501      0.00000000      0.00000000      0.00000000      0.46391121     -0.67520143      0.11999501      0.00000000       0    6566    6566       0       0       0      25      0.64414414 term
H        5.57958805     -4.62134401      1.98471136       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.01709260      0.00000000     -4.14494205   2910.85773629       1       1      0.07576912     -0.00455629      0.00097368     -0.00294810      0.00012877     -0.00026364      0.00034159   -119.59616566     -4.02677538      2.18626690   -119.45686344     -4.49274213      2.55757118      0.37484246     -0.76743296      0.99432507      0.37484246     -0.76743296      0.99432507      0.35575564     -1.40617171      1.52014871      0.00000000       0    6747    6747       0       0       0      25      0.64414414 term
H        7.58040902     -9.08757712      4.96603885       1      -6       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.04634150      0.00000000     -4.30220865   2910.85773629       1       1      0.00860020      0.00089563      0.00249568      0.00220497      0.00002958      0.00006852     -0.00005235   -743.51618063     -9.79907184      5.34432953   -743.45333110     -9.89178913      5.49053987      0.08610596      0.19946529     -0.15237367     -0.05808771      0.10589692     -0.10132644      0.08610596      0.19946529     -0.15237367      0.00000000       0   31402   31402       0       0       0      26      0.64414414 term
H        7.37039631     -6.40965691      4.97535403       1      -6       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.02639675      0.00000000     -4.18304790   2910.85773629       1       1      0.03760822     -0.00353654      0.00725853     -0.00198377     -0.00003649     -0.00002397     -0.00000884   -743.75117178     -5.79018627      5.42633898   -743.69659696     -5.84468243      5.43480427     -0.10621937     -0.06977546     -0.02574137     -0.10621937     -0.06977546     -0.02574137      0.09252634      0.39324738     -0.17314227      0.00000000       0   31406   31406       0       0       0      26      0.64414414 term
H        9.05799924     -4.95160711      1.53233042       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.02574407      0.00000000     -4.14126568   2910.85773629       1       1      0.07481067     -0.00052947     -0.00131955      0.00686207     -0.00019972     -0.00010076      0.00024500   -115.52572288     -5.01038913      1.86110663   -115.78213228     -5.39614883      1.66025709     -0.58136440     -0.29329706      0.71316950     -0.58136440     -0.29329706      0.71316950     -0.47681556      0.03184410      0.84937135      0.00000000       0    6928    6928       0       0       0      27      0.64414414 term
H        9.60576454     -5.71219411      2.41382518       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       1    F    F      0.00000000     -0.02574407      0.00000000     -4.14126568   2910.85773629       1       1      0.07481067     -0.00052947     -0.00131955      0.00686207     -0.00019972     -0.00010076      0.00024500   -115.52572288     -5.01038913      1.86110663   -115.78213228     -5.39614883      1.66025709     -0.58136440     -0.29329706      0.71316950     -0.58136440     -0.29329706      0.71316950     -0.47681556      0.03184410      0.84937135      0.00000000       0    6928    6928       0       0       0      28      0.64414414 term
H       10.74001970     -7.43396972      3.96964584       1      -5       0       0       0       0       0       1       0       0       1       4       4       0       0    F    F      0.00000000     -0.03509485      0.00000000     -4.23902095   2910.85773629       1       1      0.00811107      0.00029159      0.00042004      0.00321761      0.00003282      0.00002160      0.00011166   -113.70496336     -7.45205442      4.38119149   -113.64174799     -7.35790124      4.32227808      0.09552644      0.06287875      0.32503090      0.09552644      0.06287875      0.32503090     -0.61519270      0.42485813      0.28644278      0.00000000       0   31761   31761       0       0       0      28      0.64414414 term      """)

   def test_read_int_lattice(self):
      cio = CInOutput()
      at = cio.read(str="""6
EMP2=2.67539000e-03  ELONG=1.46920613e-03  COM=3.16174712e+00 Lattice="20 0 0 0 20 0 0 0 20" Properties=species:S:1:pos:R:3
O	8.433042	7.254127	4.000330
H	9.329070	6.902494	4.019458
H	8.505653	8.201400	4.069496
O	5.939224	8.979049	3.011011
H	5.602166	8.307183	3.670115
H	6.859570	9.215634	3.262673""")
      self.assertArrayAlmostEqual(at.lattice, np.diag([20.0,20.0,20.0]))

   def test_read_empty_file(self):
      os.system('touch empty.xyz')
      self.assertRaises(RuntimeError, Atoms, 'empty.xyz')

   def test_missing_value(self):
      cio = CInOutput()
      self.assertRaises(RuntimeError, cio.read, str="""1
Properties="species:S:1:pos:R:3" Lattice="10. 0 0 0 10.0 0 0 0 10"  State=
H 0. 0. 0.
""")

   def test_missing_value_quoted(self):
      cio = CInOutput()
      self.assertRaises(RuntimeError, cio.read, str="""1
Properties="species:S:1:pos:R:3" Lattice="10. 0 0 0 10.0 0 0 0 10"  State=""
H 0. 0. 0.
""")

   def test_one_frame_per_file_low_level(self):
      cio = CInOutput('test00000.xyz', action=OUTPUT, one_frame_per_file=True)
      self.al.write(cio)
      cio = CInOutput('test00000.xyz', action=INPUT, one_frame_per_file=True)
      al = [cio.read() for i in range(cio.n_frame)]
      self.assertEqual(al, list(self.al))

   def test_one_frame_per_file_high_level(self):
      self.al.write('test00000.xyz', one_frame_per_file=True)
      al = AtomsList('test00000.xyz', one_frame_per_file=True)
      al2 = AtomsList('test0000*.xyz', no_compute_index=False)
      self.assertEqual(al, self.al)
      self.assertEqual(al, al2)

   def test_read_xyz_indices_all(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', indices=frange(self.at.n))
      self.assertEqual(at, self.at)

   def test_read_xyz_indices_subset(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', indices=[1,5,10])
      sub = self.at.select(list=[1,5,10], orig_index=False)
      self.assertEqual(at, sub)

   def test_read_xyz_indices_empty(self):
      self.at.write('test.xyz')
      at = Atoms('test.xyz', indices=[])
      sub = self.at.select(list=[], orig_index=False)
      self.assertEqual(at, sub)

   def test_read_xyz_bad_indices(self):
      self.at.write('test.xyz')
      self.assertRaises(RuntimeError, Atoms, 'test.xyz', indices=[self.at.n+1])

   def test_read_xyz_bare_keys(self):
      self.at.params['bare_key'] = None
      self.at.write('test.xyz')
      at = Atoms('test.xyz')
      self.assertEqual(at, self.at)
      self.assert_(at.params['bare_key'] is None)

   if 'netcdf' in available_modules:

      def test_read_nc_indices(self):
         self.at.write('test.nc')
         self.assertRaises(RuntimeError, Atoms, 'test.nc', indices=[])


class TestPythonNetCDF(QuippyTestCase):
   def setUp(self):
      self.at = supercell(diamond(5.44,14), 2,2,2)
      self.al = AtomsList([ supercell(diamond(5.44+0.01*x,14),2,2,2) for x in range(5) ])
      self.al.write('test3.nc', netcdf4=False)

   def tearDown(self):
      if os.path.exists('test3.nc'): os.remove('test3.nc')
      if os.path.exists('dataset.nc'): os.remove('dataset.nc')

   if 'netcdf' in quippy.available_modules:

      def testpupynere_read(self):
         from quippy.pupynere import netcdf_file
         nc = netcdf_file('test3.nc', 'r')
         al = AtomsList(nc, format=quippy.netcdf.netcdf_file)
         self.assertEqual(list(self.al), list(al))
         nc.close()

   if 'netCDF4' in quippy.available_modules:

      def testnetcdf4_read(self):
         from netCDF4 import Dataset
         nc = Dataset('test3.nc','r')
         al = AtomsList(nc)
         for a, b in zip(self.al, al):
            self.assertEqual(a, b)
         nc.close()

      def testnetcdf4_write(self):
         from netCDF4 import Dataset
         nc = Dataset('dataset.nc','w')
         al2 = AtomsList(self.al)
         al2.write(nc)
         nc.close()
         al = AtomsList('dataset.nc')
         self.assertEqual(list(self.al), list(al))


if __name__ == '__main__':
   unittest.main()
