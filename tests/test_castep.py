# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   libAtoms+QUIP: atomistic simulation library
# HND X
# HND X   Portions of this code were written by
# HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HND X
# HND X   Copyright 2006-2010.
# HND X
# HND X   Not for distribution
# HND X
# HND X   Portions of this code were written by Noam Bernstein as part of
# HND X   his employment for the U.S. Government, and are not subject
# HND X   to copyright in the USA.
# HND X
# HND X   When using this software, please cite the following reference:
# HND X
# HND X   http://www.libatoms.org
# HND X
# HND X  Additional contributions by
# HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HND X
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from quippy import *
from quippy import castep

from StringIO import StringIO
from quippytest import *
import unittest
import xml.dom.minidom

class TestCastepCell(QuippyTestCase):
   def setUp(self):
      self.cell_lines = """fix_all_cell  TRUE
            
         %BLOCK SPECIES_POT
         H  H_00PBE.usp
         C  C_00PBE.usp
         %ENDBLOCK SPECIES_POT

         %BLOCK LATTICE_CART
         10.000000 0.000000 0.000000
         0.000000 10.000000 0.000000
         0.000000 0.000000 10.000000
         %ENDBLOCK LATTICE_CART

         %BLOCK POSITIONS_ABS
         H 0.000000 0.000000 0.000000
         %ENDBLOCK POSITIONS_ABS""".split("\n")

      self.cell = castep.CastepCell(self.cell_lines)


   def testkeys(self):
      self.assertEqual(self.cell.keys(),['fix_all_cell','SPECIES_POT', 'LATTICE_CART', 'POSITIONS_ABS'])

   def testrepr(self):
      self.assertEqual(repr(self.cell), "CastepCell({'fix_all_cell': 'TRUE', 'POSITIONS_ABS': ['H 0.000000 0.000000 0.000000'], 'SPECIES_POT': ['H  H_00PBE.usp', 'C  C_00PBE.usp'], 'LATTICE_CART': ['10.000000 0.000000 0.000000', '0.000000 10.000000 0.000000', '0.000000 0.000000 10.000000']})")

   def testvalues(self):
      self.assertEqual(self.cell.values(),['TRUE',
                                           ['H  H_00PBE.usp', 'C  C_00PBE.usp'],
                                           ['10.000000 0.000000 0.000000',
                                            '0.000000 10.000000 0.000000',
                                            '0.000000 0.000000 10.000000'],
                                           ['H 0.000000 0.000000 0.000000']])

   def testbadparseblock(self):
      bad_lines = """%BLOCK SPECIES_POT
      %BLOCK LATTICE_CART
      H  H_00PBE.usp
      C  C_00PBE.usp
      %ENDBLOCK SPECIES_POT
      """.split("\n")
      self.assertRaises(ValueError, castep.CastepCell, bad_lines)

   def testbadparseblock2(self):
      bad_lines = """%BLOCK SPECIES_POT
      %ENDBLOCK LATTICE_CART""".split("\n")
      self.assertRaises(ValueError, castep.CastepCell, bad_lines)

   def testbadkeyword(self):
      bad_lines = "BADKEYWORD = 1".split("\n")
      self.assertRaises(ValueError, castep.CastepCell, bad_lines)
      
   def testcopy(self):
      cell_copy = self.cell.copy()

      self.assertEqual(cell_copy.keys(), self.cell.keys())
      self.assertEqual(cell_copy.values(), self.cell.values())

   def testread_xml(self):
      cell_xml = xml.dom.minidom.parseString("""<params><castep_cell>
      %s
      </castep_cell></params>""" % '\n'.join(self.cell_lines))

      new_cell = castep.CastepCell(xml=cell_xml)
      self.assertEqual(new_cell.keys(), self.cell.keys())
      self.assertEqual(new_cell.values(), self.cell.values())

   def testread_xml_badxml(self):
      cell_xml = xml.dom.minidom.parseString("""<params><bad_castep_cell>
      </bad_castep_cell></params>""")
      self.assertRaises(ValueError, castep.CastepCell, xml=cell_xml)

   def testread_xml_multiple_stanzas(self):
      cell_xml = xml.dom.minidom.parseString("""<params><castep_cell>
      </castep_cell><castep_cell></castep_cell></params>""")
      self.assertRaises(ValueError, castep.CastepCell, xml=cell_xml)
      
   def testread_xml_badchild(self):
      cell_xml = xml.dom.minidom.parseString("""<params><castep_cell>
      <child/></castep_cell></params>""")
      self.assertRaises(ValueError, castep.CastepCell, xml=cell_xml)      
      

   def testwrite(self):
      f = StringIO()
      
      self.cell.write(f)
      f.seek(0)
      lines = f.readlines()
      f.close()

      self.assertEqual([line.strip() for line in lines if line.strip() != ''],
                       [line.strip() for line in self.cell_lines if line.strip() != ''])

   def testto_atoms(self):
      at = self.cell.to_atoms()
      self.assertEqual(at.species[1].stripstrings(), 'H')
      self.assertEqual(at.n, 1)
      latt_error = at.lattice - farray([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
      self.assertEqual((latt_error < 1e-8).all(), True)

   def testatomsread(self):
      at = Atoms(self.cell_lines, format='cell')
      self.assertEqual(at.species[1].stripstrings(), 'H')
      self.assertEqual(at.n, 1)
      latt_error = at.lattice - farray([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
      self.assertEqual((latt_error < 1e-8).all(), True)

   def testatomswrite(self):
      f = StringIO()
      at = Atoms(self.cell_lines, format='cell')
      at.write(f, format='cell')
      f.seek(0)
      lines = f.readlines()
      f.close()

      self.assertEqual([line.strip() for line in lines if line.strip() != ''],
                       [line.strip() for line in self.cell_lines if line.strip() != ''][5:])
      

   def test_to_atoms_without_lattice(self):
      del self.cell['LATTICE_CART']
      self.assertRaises(ValueError, self.cell.to_atoms)

   def test_to_atoms_without_pos(self):
      del self.cell['POSITIONS_ABS']
      self.assertRaises(ValueError, self.cell.to_atoms)

   def test_to_atoms_frac_pos(self):
      def convertStr(s):
         try:
            ret = float(s)
         except:
            ret = s
         return ret

      c = castep.CastepCell("""%block lattice_cart
ang
2.4866946  -4.3070816   0.0000000
2.4866946   4.3070816   0.0000000
0.0000000   0.0000000   5.4621192
%endblock lattice_cart

%block positions_frac
O     0.414382   0.259186   0.792179
O    -0.259186   0.155196   1.458845
O    -0.155196  -0.414382   1.125512
O     0.259186   0.414382  -0.792179
O    -0.414382  -0.155196  -0.125512
O     0.155196  -0.259186  -0.458845
Si    0.474452   0.000000   0.666667
Si    0.000000   0.474452   0.333333
Si   -0.474452  -0.474452   0.000000
%endblock positions_frac""".split('\n'))

      a = c.to_atoms()
      self.assertEqual(a.n, 9)
      self.assertEqual(list(a.z), [ 8,  8,  8,  8,  8,  8, 14, 14, 14])
      self.assertArrayAlmostEqual(a.pos, FortranArray([[ 1.67495791, -0.66844184,  4.32697613],
                                                       [-0.25859137,  1.78477709,  7.96838528],
                                                       [-1.41636654, -1.11633525,  6.14768071],
                                                       [ 1.67495791,  0.66844184, -4.32697613],
                                                       [-1.41636654,  1.11633525, -0.68556151],
                                                       [-0.25859137, -1.78477709, -2.50626608],
                                                       [ 1.17981723, -2.04350348,  3.64141462],
                                                       [ 1.17981723,  2.04350348,  1.82070458],
                                                       [-2.35963445,  0.        ,  0.        ]]).T)
      c2 = castep.CastepCell()
      c2.update_from_atoms(a, frac_pos=True)
      self.assertEqual([[convertStr(r) for r in x.split()] for x in c['POSITIONS_FRAC']], [[convertStr(r) for r in x.split()] for x in c2['POSITIONS_FRAC']])

   def test_to_atoms_unit_bohr(self):
      c = castep.CastepCell("""%block lattice_cart
bohr
 4.69917141  -8.13920403   0.       
 4.69917141   8.13920403   0.       
 0.           0.          10.3219086
%endblock lattice_cart

%block positions_abs
bohr
O    3.16521149  -1.26317191   8.17679924
O   -0.48866684   3.37273965  15.05806476
O   -2.67654465  -2.10956774  11.617432  
O    3.16521149   1.26317191  -8.17679924
O   -2.67654465   2.10956774  -1.29552339
O   -0.48866684  -3.37273965  -4.73615615
Si   2.22953127  -3.86166163   6.88127584
Si   2.22953127   3.86166163   3.44063276
Si  -4.45906255   0.           0.        
%endblock positions_abs""".split('\n'))

      a = c.to_atoms()
      self.assertEqual(a.n, 9)
      self.assertEqual(list(a.z), [ 8,  8,  8,  8,  8,  8, 14, 14, 14])
      self.assertArrayAlmostEqual(a.pos, FortranArray([[ 1.67495791, -0.66844184,  4.32697613],
                                                       [-0.25859137,  1.78477709,  7.96838528],
                                                       [-1.41636654, -1.11633525,  6.14768071],
                                                       [ 1.67495791,  0.66844184, -4.32697613],
                                                       [-1.41636654,  1.11633525, -0.68556151],
                                                       [-0.25859137, -1.78477709, -2.50626608],
                                                       [ 1.17981723, -2.04350348,  3.64141462],
                                                       [ 1.17981723,  2.04350348,  1.82070458],
                                                       [-2.35963445,  0.        ,  0.        ]]).T)


   def testupdate_from_atoms(self):
      a = diamond(5.44, 14)
      c = castep.CastepCell()
      c.update_from_atoms(a)
      self.assertEqual(c['POSITIONS_ABS'], ['Si 0.000000 0.000000 0.000000',
                                            'Si 1.360000 1.360000 1.360000',
                                            'Si 2.720000 2.720000 0.000000',
                                            'Si 4.080000 4.080000 1.360000',
                                            'Si 2.720000 0.000000 2.720000',
                                            'Si 4.080000 1.360000 4.080000',
                                            'Si 0.000000 2.720000 2.720000',
                                            'Si 1.360000 4.080000 4.080000'])
      self.assertEqual(c['LATTICE_CART'], ['5.440000 0.000000 0.000000',
                                           '0.000000 5.440000 0.000000',
                                           '0.000000 0.000000 5.440000'])

   def testnonorthorhombic(self):
      import os
      
      shear = fidentity(3)
      shear[1,2] = 0.05
      a = transform(diamond(5.44, 14), shear)
      c = castep.CastepCell(atoms=a)
      a2 = c.to_atoms()
      self.assertEqual(a, a2)

      c.write('test.cell')

      d = castep.CastepCell('test.cell')
      a3 = d.to_atoms()
      self.assertEqual(a, a3)
      
      os.remove('test.cell')


class TestCastepParam(QuippyTestCase):

   def setUp(self):
      self.param_str = """#comment line
task                    : SinglePoint
xc_functional           : PBE
opt_strategy            : speed
fix_occupancy           : false
popn_calculate          : false
finite_basis_corr       : none
spin_polarised          : false
max_scf_cycles          : 60
num_dump_cycles         : -1
iprint                  : 3
metals_method           : DM
elec_energy_tol		: 1e-4
cut_off_energy          : 200
devel_code              : EWALD:PREC=25.0:OLD=F:R2R=0.004:ENDEWALD
mix_charge_amp		: 0.25"""

      self.castep_output = """ +-------------------------------------------------+
 |                                                 |
 |      CCC   AA    SSS  TTTTT  EEEEE  PPPP        |
 |     C     A  A  S       T    E      P   P       |
 |     C     AAAA   SS     T    EEE    PPPP        |
 |     C     A  A     S    T    E      P           |
 |      CCC  A  A  SSS     T    EEEEE  P           |
 |                                                 |
 +-------------------------------------------------+
 |                                                 |
 | Welcome to Academic Release CASTEP version 4.4  |
 | Ab Initio Total Energy Program                  |
 |                                                 |
 | Authors:                                        |
 | M. Segall, M. Probert, C. Pickard, P. Hasnip,   |
 | S. Clark, K. Refson, M. Payne                   |
 |                                                 |
 | Contributors:                                   |
 | P. Lindan, P. Haynes, J. White, V. Milman,      |
 | N. Govind, M. Gibson, P. Tulip, V. Cocula,      |
 | B. Montanari, D. Quigley, M. Glover,            |
 | L. Bernasconi, A. Perlov, M. Plummer            |
 |                                                 |
 | Copyright (c) 2000 - 2008                       |
 |                                                 |
 |     Distributed under the terms of an           |
 |     Agreement between the United Kingdom        |
 |     Car-Parrinello (UKCP) Consortium,           |
 |     Daresbury Laboratory and Accelrys, Inc.     |
 |                                                 |
 | Please cite                                     |
 |                                                 |
 |     "First principles methods using CASTEP"     |
 |                                                 |
 |         Zeitschrift fuer Kristallographie       |
 |           220(5-6) pp. 567-570 (2005)           |
 |                                                 |
 | S. J. Clark, M. D. Segall, C. J. Pickard,       |
 | P. J. Hasnip, M. J. Probert, K. Refson,         |
 | M. C. Payne                                     |
 |                                                 |
 |       in all publications arising from          |
 |              your use of CASTEP                 |
 |                                                 |
 +-------------------------------------------------+
 

 ******************************* User Parameters *******************************
iprint          : 3                                       # verbosity control (I)
reuse           : NULL                                    # reuse filename (S)
task            : SinglePoint                             # type of calculation (S)
optstrategy     : speed                                   # optimization strategy (S)
xcfunctional    : PBE                                     # exchange-correlation functional (S)
cutoffenergy    : 200.0000000 ev                          # maximum energy of planewaves in basis set (P)
finitebasiscorr : none                                    # finite basis set correction (S)
spinpolarised   : F                                       # spin polarised (L)
metalsmethod    : DM                                      # treatment of metals or finite temperature insulator (S)
elecenergytol   : 0.0001000 ev                            # total energy per atom convergence tolerance (P)
maxscfcycles    : 60                                      # maximum SCF cycles (I)
fixoccupancy    : F                                       # treat system as an insulator (L)
numdumpcycles   : -1                                      # frequency of wavefunction dumps (I)
mixchargeamp    : 0.2500000                               # charge density mixing amplitude (R)
popncalculate   : F                                       # population analysis on/off (L)
develcode       : EWALD:PREC=25.0:OLD=F:R2R=0.004:ENDEWALD # developers code (S)

 *******************************************************************************
  
"""
      self.castep_lines = [s+"\n" for s in self.castep_output.split("\n")]


      self.param_lines = self.param_str.split("\n")
      self.param_xml = '<p><castep_param label="DM" order="0">%s</castep_param></p>' % self.param_str

      self.param = castep.CastepParam(self.param_lines)

   def testkeys(self):
      self.assertEqual(self.param.keys(), ['task',
                                           'xc_functional',           
                                           'opt_strategy',            
                                           'fix_occupancy',           
                                           'popn_calculate',          
                                           'finite_basis_corr',       
                                           'spin_polarised',          
                                           'max_scf_cycles',         
                                           'num_dump_cycles',         
                                           'iprint',
                                           'metals_method',
                                           'elec_energy_tol',         
                                           'cut_off_energy',          
                                           'devel_code',              
                                           'mix_charge_amp'])

   def testvalues(self):
      self.assertEqual(self.param.values(), ['SinglePoint',
                                             'PBE',
                                             'speed',
                                             False,
                                             False,
                                             'none',
                                             False,
                                             '60',
                                             '-1',
                                             '3',
                                             'DM',
                                             '1e-4',
                                             '200',
                                             'EWALD:PREC=25.0:OLD=F:R2R=0.004:ENDEWALD',
                                             '0.25'])

   def testreadxml(self):
      xmlparam = castep.CastepParam(xml=xml.dom.minidom.parseString(self.param_xml))
      self.assertEqual(self.param.keys(), xmlparam.keys())
      self.assertEqual(self.param.values(), xmlparam.values())

   def testcopy(self):
      cp = self.param.copy()
      self.assertEqual(self.param.keys(), cp.keys())
      self.assertEqual(self.param.values(), cp.values())

   def testreadbadparam(self):
      self.assertRaises(ValueError, castep.CastepParam, ['badparam: false'])

   def testread_xml_badxml(self):
      param_xml = xml.dom.minidom.parseString("""<params><bad_castep_param>
      </bad_castep_param></params>""")
      self.assertRaises(ValueError, castep.CastepParam, xml=param_xml)

   def testread_xml_multiple_stanzas(self):
      param_xml = xml.dom.minidom.parseString("""<params><castep_param>
      </castep_param><castep_param></castep_param></params>""")
      self.assertRaises(ValueError, castep.CastepParam, xml=param_xml)
      
   def testread_xml_badchild(self):
      param_xml = xml.dom.minidom.parseString("""<params><castep_param>
      <child/></castep_param></params>""")
      self.assertRaises(ValueError, castep.CastepParam, xml=param_xml)      

   def testread_castep_output(self):
      p = castep.CastepParam()
      p.read_from_castep_output(self.castep_lines)

      rep = [('ev','')]
      
      for k,v in self.param.iteritems():
         v2 = p[k]
         for a,b in rep:
            if isinstance(v2, str):
               v2 = v2.replace(a,b)
         try:
            v = eval(v)
            v2 = eval(v2)
         except:
            pass
         self.assertEqual(v, v2)

   def testwrite(self):
      f = StringIO()
      
      self.param.write(f)
      f.seek(0)
      lines = f.readlines()
      f.close()

      for a,b in zip([[x for x in line.strip().split() if x not in (':','=')] for line in lines if line.strip() != ''],
                     [[x for x in line.strip().split() if x not in (':','=')] for line in self.param_lines
                      if (line.strip() != '' and not line.startswith('#'))]):
         self.assertEqual(a,b)

   def testparsernospaces1(self):
      p = castep.CastepParam(['#comment',
                       'elec_energy_tol:value'])
      self.assertEqual(p.keys(), ['elec_energy_tol'])
      self.assertEqual(p.values(), ['value'])

   def testparsernospaces2(self):
      p = castep.CastepParam(['#comment',
                       'elec_energy_tol: value'])
      self.assertEqual(p.keys(), ['elec_energy_tol'])
      self.assertEqual(p.values(), ['value'])

   def testparsernospaces3(self):
      p = castep.CastepParam(['#comment',
                       'elec_energy_tol :value'])
      self.assertEqual(p.keys(), ['elec_energy_tol'])
      self.assertEqual(p.values(), ['value'])

   def testfortrancomment(self):
      p = castep.CastepParam(['!comment',
                              'elec_energy_tol :value'])
      self.assertEqual(p.keys(), ['elec_energy_tol'])
      self.assertEqual(p.values(), ['value'])
      
      

   
class TestReadGeom(QuippyTestCase):

   def setUp(self):
      self.geom_string = """ BEGIN header
  
 END header
  
                                      0
                     -1.1024193863342781E+002   -1.1024193863342781E+002                             <-- E
                      1.2314768852609946E+001    0.0000000000000000E+000    0.0000000000000000E+000  <-- h
                      0.0000000000000000E+000    7.1099339546491658E+000    0.0000000000000000E+000  <-- h
                      0.0000000000000000E+000    0.0000000000000000E+000    7.8391993010819832E+001  <-- h
                     -1.0690207662406373E-004    1.7101111405881319E-006   -1.5920944508643492E-006  <-- S
                      1.7101111405881319E-006   -1.0105502034559073E-004   -1.6011442565378357E-006  <-- S
                     -1.5920944508643492E-006   -1.6011442565378357E-006   -8.3612739601073363E-005  <-- S
 Si              1    1.3120122888016088E+000    1.7780470994819315E+000    3.9027321440363131E+001  <-- R
 Si              2    7.1395440734519555E+000    5.3315315866539335E+000    3.9081817362632258E+001  <-- R
 Si              3    9.3186742066424362E+000    1.7774839610938253E+000    3.3351375927324533E+001  <-- R
 Si              4    3.1361139039094752E+000    5.3324499935553424E+000    3.3020241106487596E+001  <-- R
 Si              5    5.1311544759733669E+000    1.7774839610938253E+000    2.7300503080006163E+001  <-- R
 Si              6    1.1288536067689135E+001    5.3324499935553424E+000    2.7300503080006163E+001  <-- R
 Si              7    1.0262290054685388E+000    1.7774839610938253E+000    2.1495264394565503E+001  <-- R
 Si              8    7.1836162663627121E+000    5.3324499935553424E+000    2.1495264394565503E+001  <-- R
 Si              9    9.2360742772997888E+000    1.7774839610938253E+000    1.5690025709124834E+001  <-- R
 Si             10    3.0786907958578849E+000    5.3324499935553424E+000    1.5690025709124834E+001  <-- R
 Si             11    5.1705212508090099E+000    1.7774839610938253E+000    9.7185856101475459E+000  <-- R
 Si             12    1.1061952234974164E+001    5.3324499935553424E+000    9.7705436302179223E+000  <-- R
 Si             13    6.8242149220046722E+000    1.7774839610938253E+000    5.5487820619901962E+000  <-- R
 Si             14    7.5995903420117816E-001    5.3324499935553424E+000    3.1957045278374334E+000  <-- R
 Si             15    1.1230978789081989E+001    5.3548621455115031E+000    4.0661815494031430E+001  <-- R
 Si             16    5.2447780392674721E+000    1.6986086819624426E+000    4.0940393251613322E+001  <-- R
 Si             17    1.2528128381835097E+000    1.7774839610938253E+000    3.4671292940547119E+001  <-- R
 Si             18    7.1442457120747980E+000    5.3324499935553424E+000    3.4723250960617477E+001  <-- R
 Si             19    9.2360761670259244E+000    1.7774839610938253E+000    2.8751810861640198E+001  <-- R
 Si             20    3.0786907958578849E+000    5.3324499935553424E+000    2.8751810861640198E+001  <-- R
 Si             21    5.1311544759733669E+000    1.7774839610938253E+000    2.2946574065925663E+001  <-- R
 Si             22    1.1288536067689135E+001    5.3324499935553424E+000    2.2946574065925663E+001  <-- R
 Si             23    1.0262290054685388E+000    1.7774839610938253E+000    1.7141335380485007E+001  <-- R
 Si             24    7.1836162663627121E+000    5.3324499935553424E+000    1.7141335380485007E+001  <-- R
 Si             25    9.1786530589743336E+000    1.7774839610938253E+000    1.1421595464277434E+001  <-- R
 Si             26    2.9960946459675073E+000    5.3324499935553424E+000    1.1090458753714360E+001  <-- R
 Si             27    2.8719283005685607E+000    1.7774839610938253E+000    3.8215005653877174E+000  <-- R
 Si             28    9.4175238910163692E+000    5.3324499935553424E+000    5.6061976111372491E+000  <-- R
 Si              1   -1.9922292705419436E-002    1.9425978802557920E-003    2.4401933737649269E-002  <-- F
 Si              2   -3.5912638452564364E-003   -6.7281387656061382E-003    1.5729240482904947E-002  <-- F
 Si              3    6.1364358378183060E-003    2.5415370102216579E-004   -1.2227725215869616E-002  <-- F
 Si              4    3.4919714580724588E-003    3.8681930969282965E-005    1.4727990616849021E-002  <-- F
 Si              5    1.0338080606808078E-003    1.4625704828549854E-005    2.4225852081862388E-003  <-- F
 Si              6   -9.1953565038260668E-004    3.8332046980592605E-005    3.3028156826770813E-003  <-- F
 Si              7    7.5191454715354406E-004    2.2226930229417235E-005    1.6122847768421588E-003  <-- F
 Si              8   -9.4627137166587759E-004    1.1524229696366736E-005    1.8173648640441426E-003  <-- F
 Si              9    3.2989247820703563E-004    1.4053441191491703E-005    9.4894550143690241E-003  <-- F
 Si             10   -1.7494169864933607E-003    1.6441911692346736E-005   -1.6934412343769852E-002  <-- F
 Si             11   -9.4530449772623750E-004    1.5141427378280058E-005   -6.1919951959722545E-003  <-- F
 Si             12    6.0637167600436316E-003    1.4746712188161918E-005   -4.1089237257320062E-003  <-- F
 Si             13    2.5173427682640055E-003    1.3645883749411781E-005   -4.7108218178698249E-004  <-- F
 Si             14   -4.4897940039217961E-003    1.3779538413977732E-005   -1.0368589068478037E-002  <-- F
 Si             15    1.7174887439929384E-002   -1.4257093082497632E-003   -1.9514538822384336E-002  <-- F
 Si             16    4.1600525645794648E-003    6.1523543653281458E-003   -8.2208149845951119E-003  <-- F
 Si             17   -1.6447435119013587E-002    8.4332560678299380E-005   -3.9178462959600967E-003  <-- F
 Si             18    7.3632183367870639E-003   -5.5208702634572805E-004   -7.7689232456340223E-003  <-- F
 Si             19   -3.4484085320728331E-005   -8.1431572223518989E-005    1.9948006093133203E-002  <-- F
 Si             20   -3.5618633659583413E-005    4.3616057589002283E-005   -8.2290932871435761E-003  <-- F
 Si             21   -1.1284725095408576E-004    2.7728980060445610E-005    1.5981878117899849E-003  <-- F
 Si             22    2.1001112111821265E-004   -1.4094085546831510E-005    8.2234929590574782E-004  <-- F
 Si             23   -4.9499111276618156E-004    6.4870353026348692E-006    3.8717056356716081E-004  <-- F
 Si             24    6.7711774919394134E-004    1.8714992637719303E-005    1.0373900853276562E-003  <-- F
 Si             25   -1.5834140913376965E-003    1.5311961244285903E-005   -5.5808627400879629E-004  <-- F
 Si             26   -3.6389508575841423E-003    1.2275931418144633E-005    1.1849755833669704E-003  <-- F
 Si             27    5.6713209529111281E-003    1.5271199611201598E-005    7.6509393332644532E-004  <-- F
 Si             28   -6.7006986325723319E-004    1.5416335506263813E-005   -7.3481310860435203E-004  <-- F
  
"""
      self.geom_lines = self.geom_string.split("\n")

      self.at, = castep.CastepGeomReader(self.geom_lines)

   def testspecies(self):
      self.assert_(all(self.at.species.stripstrings() == 'Si'))
      self.assertEqual((self.at.species.stripstrings() == 'Si').sum(), 28)
      self.assertEqual(self.at.n, 28)

   def testatomsread(self):
      a = Atoms(self.geom_lines, format='geom')
      self.assertEqual(self.at, a)

   def testatomsreader(self):
      a, = AtomsList(self.geom_lines, format='geom')
      self.assertEqual(self.at, a)

class TestReadMD(QuippyTestCase):
   def setUp(self):
      self.md_lines = """ BEGIN header
  
 END header
  
                      0.0000000000000000E+000
                     -9.2387338525880364E+002   -9.2377220551495782E+002    1.0117974384583524E-001  <-- E
                      9.5004454315338255E-004                                                        <-- T
                      2.0035758980959645E+001    0.0000000000000000E+000    0.0000000000000000E+000  <-- h
                      0.0000000000000000E+000    1.8428996658316365E+001    0.0000000000000000E+000  <-- h
                      0.0000000000000000E+000    0.0000000000000000E+000    1.9027293401156189E+001  <-- h
 O               1    1.5574319736307456E+000    2.9859895100366080E-001    2.9825629648594973E+000  <-- R
 O               2    4.7824599161643224E+000   -1.7744144388914538E+000   -9.5366449833963802E+000  <-- R
 O               3   -3.5230945219665903E+000   -5.3660880626537981E+000   -6.5994562468008695E+000  <-- R
 O               4   -3.7411421156673561E+000   -7.6750681235605764E+000    1.1196954281342777E-001  <-- R
 O               5    2.2458878183413491E+000   -1.9920710124505006E+000   -1.3099078907535422E+000  <-- R
 O               6   -6.1235181930072553E+000    1.9013793905244636E+000    5.5717337489802576E-001  <-- R
 O               7    8.4644827104715823E+000   -7.8892850037450044E+000   -5.5141075359747109E+000  <-- R
 O               8    3.3491767406007504E+000    7.1813666407462087E+000   -9.2804031883700038E+000  <-- R
 O               9    1.1581775359812359E+000   -3.4889593441920135E+000   -6.1211068449794661E+000  <-- R
 O              10   -4.6210994126094551E-001   -9.3891327168805514E+000    6.0536319763583020E+000  <-- R
 O              11    9.3637750623251286E+000    8.5528349905035448E+000    9.2546120369924996E+000  <-- R
 O              12    5.2401203843403312E+000   -5.1439545284027135E+000   -3.2639486021436928E+000  <-- R
 O              13   -6.2276999873606593E+000   -1.2288400542280480E+000   -6.7389266738096092E+000  <-- R
 O              14    2.3459263441963225E+000    2.9873529815273407E+000   -1.1271803969520151E+000  <-- R
 O              15   -7.6120994308100265E+000    3.1275260988818854E+000   -9.0525187843438673E+000  <-- R
 O              16   -9.4966839982173106E+000   -5.6925226183367394E-001    2.8666853315581289E+000  <-- R
 O              17   -1.5515517347121865E+000   -4.9386017566912832E+000    5.4617334768002452E+000  <-- R
 O              18    6.2383609321826139E+000    5.6308553766364575E+000    6.3236913224823486E+000  <-- R
 O              19    7.1191910597965578E+000   -8.9461015176570733E+000   -6.6170207144462534E-001  <-- R
 O              20    2.4211845029538659E+000    2.6350204547875742E+000    7.3622379454056253E+000  <-- R
 O              21    3.3509246391740084E+000   -5.6152970516182972E+000    5.9829587380454825E+000  <-- R
 O              22   -8.8186764133811379E-001    4.3750771558866051E+000    4.4331962112215439E+000  <-- R
 O              23   -4.4428630321826788E+000    7.2458932062331192E+000    3.1981743385036547E+000  <-- R
 O              24    5.6263890405939421E+000    3.0549402987976948E+000   -8.3015404497264207E+000  <-- R
 O              25   -8.8268133109266884E+000    5.8934639116416747E+000    5.9155107363509787E+000  <-- R
 O              26    3.2975018723935795E-001   -8.3106836585631765E+000   -7.3450833569773719E+000  <-- R
 O              27   -5.3354004618459729E+000   -3.1488263497292630E+000    8.1252781990521559E+000  <-- R
 O              28    9.3754203016646209E-001   -7.2576764946557066E+000    1.9543256461283378E+000  <-- R
 O              29   -3.2544323373035772E-001   -6.6705576192148843E+000   -2.9754556655483140E+000  <-- R
 O              30   -6.1668661877643816E+000    1.9506517693592795E+000    5.4766864901287047E+000  <-- R
 O              31    8.1014679350792065E+000   -8.9257380373177426E+000    4.3351466784873658E+000  <-- R
 O              32    9.2357879994708920E+000    4.6252010720406016E-002   -7.1892911635502044E+000  <-- R
 O              33    8.2243381091960419E+000   -4.6757391982083458E+000    1.7701395929585584E+000  <-- R
 O              34   -7.3827178809038072E+000    2.6376798580609933E+000   -4.0841871246479835E+000  <-- R
 O              35    4.9877856443461379E+000   -5.6302287335079526E-001   -4.8887885293713031E+000  <-- R
 O              36    5.7251409453372371E-001    5.4535517822105550E+000   -5.2934442205214109E+000  <-- R
 O              37    7.2772147864171011E+000   -5.6727512476142632E+000    9.1136735166152878E+000  <-- R
 O              38   -2.5776463636856892E+000    3.0943523531903452E+000   -2.5908180749765863E+000  <-- R
 O              39   -6.0176254799737379E+000    6.4542357272158455E+000   -1.3971095326321388E+000  <-- R
 O              40   -4.1978664913404344E+000    5.9345708648887507E+000    7.7497050310066475E+000  <-- R
 O              41   -8.4785367146789721E+000   -8.6209292200150429E+000    1.0391816614145393E+000  <-- R
 O              42   -1.3903373445254497E+000   -2.8473373722480544E-001    6.7226667426017839E+000  <-- R
 O              43    6.1333670419540089E+000    9.8054315624444194E-002    9.3552358822475035E-001  <-- R
 O              44   -4.6354514949892627E+000   -2.1502817371158796E+000    3.2262715190676756E+000  <-- R
 O              45   -8.1334866483685850E+000   -5.6784227057297034E+000   -8.0120126329134376E+000  <-- R
 O              46    4.0367914167311891E+000    8.7613735898776337E+000   -4.3449818123145016E+000  <-- R
 O              47    6.6147808854916210E+000   -2.0163912842525833E+000    5.2558693055039400E+000  <-- R
 O              48   -9.7520433633192863E-003    7.3278599762991865E+000   -8.8443218567972748E-001  <-- R
 Si              1    4.3653557596138617E+000    4.5922966602807147E+000    8.5158138239868855E+000  <-- R
 Si              2    8.8431783751272697E+000   -7.5707715792868484E+000    1.6266371129691324E+000  <-- R
 Si              3    9.3043048782187032E+000   -7.2117231677377811E+000   -8.3011212614478218E+000  <-- R
 Si              4    6.2421067323369703E+000    1.4058729762385963E-001   -7.6751751541906907E+000  <-- R
 Si              5   -8.0081187839741403E+000    1.1629797140131717E+000   -6.7458399638165512E+000  <-- R
 Si              6   -7.8131615660019682E-001   -8.1509893769760264E+000   -4.0110037698785550E-001  <-- R
 Si              7   -6.6389745022265982E+000    2.9828595663141988E-001    2.9229552383475359E+000  <-- R
 Si              8   -6.5711957315462499E-001   -5.8680322204165014E+000   -5.9610570455425753E+000  <-- R
 Si              9   -3.2416320307553490E+000   -2.6369297922767263E+000    5.9392065321761072E+000  <-- R
 Si             10    3.3344132182170201E+000   -2.8803979308791519E+000   -3.9621140957242420E+000  <-- R
 Si             11    3.4424638628125843E-002    4.6820629329020536E+000   -2.4048055552003662E+000  <-- R
 Si             12    8.7092286723552785E+000    7.4515515463249047E+000    6.4661180016160111E+000  <-- R
 Si             13   -2.4468332121734249E+000    6.7545622870888344E+000    5.3608960407328219E+000  <-- R
 Si             14   -5.7539856393718312E+000    8.5648265262909344E+000    7.9516381595950980E-001  <-- R
 Si             15    3.9259128363081092E-001    1.7833306238754252E+000    5.3225633335248199E+000  <-- R
 Si             16    7.6055843712536741E-001   -6.7278316145594967E+000    4.9463516433194661E+000  <-- R
 Si             17    7.8550172956638056E+000   -1.8287883502630864E+000    2.6297198787657630E+000  <-- R
 Si             18   -6.6231745327796645E+000    4.2872192363561243E+000    7.3253176773246684E+000  <-- R
 Si             19   -5.5708592240809507E+000    3.5056144811171572E+000   -1.9809310656645285E+000  <-- R
 Si             20    6.3228031992450981E+000   -7.9787080734941780E+000   -3.3161373523105691E+000  <-- R
 Si             21    5.4637468265967648E+000   -3.9704546050154925E+000    7.4400637322093957E+000  <-- R
 Si             22    3.0649436495326516E+000    5.2083865018656361E-001    3.5108434057025178E-001  <-- R
 Si             23    2.1694399579319446E+000    7.8235735101072397E+000   -6.5909767006852125E+000  <-- R
 Si             24   -5.7389701648415690E+000   -3.9142266372981123E+000   -7.9432893578058570E+000  <-- R
 O               1    9.9962249904713300E-005    1.9908363788661700E-004   -2.5466050427574749E-004  <-- V
 O               2   -2.1824345521718096E-005    1.4583543841869385E-004   -2.2755027745507153E-005  <-- V
 O               3   -1.9980125623920753E-004    1.2283508379539552E-004   -2.4855573502788282E-004  <-- V
 O               4    1.0854063701026119E-004   -4.4569550085162140E-005    1.3240809923333860E-004  <-- V
 O               5   -2.2319204610550921E-004   -3.4307825049854766E-005    2.7909386737235759E-004  <-- V
 O               6    2.0966023447212707E-004    2.3752031215641572E-004    8.5751773958941437E-005  <-- V
 O               7   -1.5506760426890263E-005   -2.3824941115798472E-004    1.3266222340499618E-005  <-- V
 O               8   -7.9602763615327086E-005   -2.9719252541528560E-004    1.1014581782090857E-004  <-- V
 O               9    3.3273100838750657E-004   -2.1983271691999019E-004    2.5976294828710135E-005  <-- V
 O              10    3.5967431932418966E-005   -3.1142624953605340E-004   -4.2918110962020101E-005  <-- V
 O              11    5.2923281063510330E-005    1.9166603913700090E-005    1.3638293402614015E-004  <-- V
 O              12    2.9917924149102420E-005    1.9303700722332471E-004   -7.1135944229446289E-005  <-- V
 O              13    6.4230825381313033E-005    3.3654514510506674E-004   -4.8598083675892879E-005  <-- V
 O              14   -1.0906847400203236E-004    1.3164568186547411E-004   -2.9443177250829393E-004  <-- V
 O              15   -3.5384151685324286E-005   -7.5928205056264119E-005    2.4654034688807858E-005  <-- V
 O              16   -9.9557362924056756E-006   -1.1271683392736717E-004    1.7229679776409164E-004  <-- V
 O              17   -1.4498759458163974E-004   -8.3721499003015104E-005    3.1607803990166316E-004  <-- V
 O              18   -2.2564995327493828E-004    3.2292523840146115E-005   -3.0780320803772104E-005  <-- V
 O              19   -3.4477055758051506E-004    6.3221692919363919E-005    2.1055563904925033E-004  <-- V
 O              20    3.5230807778033904E-004   -6.3042280365797470E-005   -2.7027646292502043E-004  <-- V
 O              21   -4.8455288687729656E-006   -1.6294685270802259E-004   -1.4415830573173192E-004  <-- V
 O              22   -1.8745945984083539E-004   -1.2348989532472702E-004   -2.5061967053003926E-004  <-- V
 O              23    1.6416949149988298E-004   -6.3355480328348092E-005   -2.1804229670270022E-004  <-- V
 O              24    3.1677110010819866E-006    4.5403412140162119E-005   -1.9658493714780075E-005  <-- V
 O              25   -9.2050441026262630E-005    2.6204814848507560E-004   -3.3959448084580621E-005  <-- V
 O              26   -4.4111379587714727E-005    2.1124546264082627E-004   -3.6642518146272324E-004  <-- V
 O              27    2.8758864025608037E-005   -1.1180933709428613E-004    8.2166333426821417E-005  <-- V
 O              28    7.1377295246943177E-005   -2.7344365276240674E-004   -2.5658003062321855E-004  <-- V
 O              29   -1.5263534507642790E-005   -8.2358715602887385E-005    2.5213012569330476E-004  <-- V
 O              30   -3.9213030777441304E-005   -1.6012569334212677E-004    2.5929897709299672E-004  <-- V
 O              31    1.3793974230205249E-004    8.4219646193880822E-005   -1.4784993727755762E-004  <-- V
 O              32    1.2304756525339882E-004    2.8832862708024105E-004   -1.2613176233395868E-004  <-- V
 O              33    1.7473467984898941E-005   -4.3489112527343347E-004   -1.2552927423215506E-004  <-- V
 O              34    2.0047564524773654E-004   -5.5529803039165773E-005   -1.0068255967656893E-004  <-- V
 O              35    6.9978873588445994E-005   -7.4032275722065905E-005    3.9734458032696809E-005  <-- V
 O              36    1.0689224700476891E-004   -1.1846160657975216E-004   -7.6189238930930405E-005  <-- V
 O              37    7.5108585392896885E-004   -1.1160734632400855E-004   -9.5004001503921865E-005  <-- V
 O              38    2.4648768250456918E-005    5.9321399407196901E-005   -1.6195778200994791E-004  <-- V
 O              39   -2.4208689913698614E-004   -5.6142874801187699E-005   -2.1235307615505265E-005  <-- V
 O              40    3.0836892635161895E-004   -3.3712434851189426E-004    5.1518763584847409E-004  <-- V
 O              41   -2.6690233131961873E-004   -2.3255201613587001E-005    2.0335741746145477E-004  <-- V
 O              42   -2.3185160186798821E-004   -3.1372423907358421E-004   -1.8844489532320909E-004  <-- V
 O              43   -4.0067160672332162E-006    2.1217645204934459E-005    2.4945713584303573E-004  <-- V
 O              44    1.3353251498462277E-004   -2.8184821756718976E-005   -2.3753238431767679E-004  <-- V
 O              45   -1.0018505066216647E-004   -1.2020023120846364E-004    6.0156382297556667E-005  <-- V
 O              46   -1.3510882159784181E-004   -2.8144543095536537E-004   -2.9846439032653847E-004  <-- V
 O              47   -1.2314929948024458E-004   -1.1020991872588233E-004   -1.4564182065715590E-004  <-- V
 O              48   -2.5919014819453918E-004   -5.4853927001870073E-005    1.7521606553566564E-004  <-- V
 Si              1    1.6110966850658522E-005   -1.3234684081284541E-004    1.1964106533583512E-004  <-- V
 Si              2   -2.3332784954351447E-004    1.7698569307532624E-004    2.1508791873507910E-004  <-- V
 Si              3    4.7291023585504364E-006    1.5969346677479864E-004   -4.7458875485476323E-005  <-- V
 Si              4    1.3267108773227178E-005    7.6077162633864556E-005    1.4653387393265663E-005  <-- V
 Si              5   -1.4784627442788545E-004    7.7724035845779852E-005    1.1889446262136406E-004  <-- V
 Si              6   -1.2225151989325684E-004   -8.4899431995521144E-005    7.4418150923932577E-005  <-- V
 Si              7    3.4687847520197701E-005   -5.0689627694374290E-005    2.6787528968495838E-007  <-- V
 Si              8    5.8301237251150769E-005    1.2160064862273995E-004   -2.4505936464514864E-005  <-- V
 Si              9   -1.6258391891218826E-005   -2.5562977030092407E-004    2.2589432604191013E-004  <-- V
 Si             10    1.0229694754060861E-005   -9.7699032087193282E-005    5.6320158858574600E-005  <-- V
 Si             11    1.4833673546201843E-004    6.8492158615930487E-006    9.1952652546084703E-005  <-- V
 Si             12    1.1983847549174861E-004   -2.9587591578156683E-005   -7.5337123178116836E-005  <-- V
 Si             13    1.0142431536220775E-005    8.0669949417164877E-005    8.6838187957358510E-005  <-- V
 Si             14   -1.8521947423455498E-004   -6.6444666090773863E-005    9.1540054799528096E-005  <-- V
 Si             15   -7.7121758274447534E-006    2.2254882190160530E-004    4.4352314182428282E-005  <-- V
 Si             16    1.1390583276529631E-004    3.9504549613616886E-005    1.9940315139994617E-004  <-- V
 Si             17    4.2429820162277302E-005    2.0367798985418200E-004   -1.1840356495160194E-004  <-- V
 Si             18    4.3250997981509499E-005    1.3948328860229062E-004    9.9740140204449576E-006  <-- V
 Si             19    4.7908577554786383E-007    1.3870688482255225E-004   -1.2715921977459916E-004  <-- V
 Si             20    1.9640147344598350E-005    4.5960437230606945E-005   -1.2969544600819705E-004  <-- V
 Si             21    2.9309351065819148E-004    2.4821463975887892E-004    3.2264338293282975E-005  <-- V
 Si             22   -1.0094408714709732E-004    2.4550191299613343E-007   -1.2461342824654671E-004  <-- V
 Si             23   -1.6222278054761326E-004   -1.0208795173918208E-005   -3.3892062915773204E-005  <-- V
 Si             24   -1.0760476751198925E-004    2.0022816754643486E-004   -1.5645811674112192E-004  <-- V
 O               1    6.5930020453396397E-003    5.4981386843022317E-003   -8.1497942728569195E-003  <-- F
 O               2    4.2494292456798132E-003   -3.5310219778806129E-003   -5.7248812880128230E-003  <-- F
 O               3    4.8754992277283181E-003    2.5502207589686776E-003   -3.6367121645393551E-003  <-- F
 O               4   -8.4141108483867470E-003   -1.0769122560083257E-002    1.2574598653774243E-003  <-- F
 O               5   -1.1219681908743917E-002    2.3695990705585666E-002    2.0575593953132590E-002  <-- F
 O               6    1.0606298064220466E-002    2.1130670979717818E-002   -3.7879066927440266E-002  <-- F
 O               7   -5.2553004984099662E-003   -2.5594648557930666E-003    2.1529515626146509E-002  <-- F
 O               8   -6.9415921108810107E-004   -1.2776391676787825E-002    4.7730089424537075E-003  <-- F
 O               9    1.6573647303801457E-002    8.7028169397755581E-003    6.6714368011531328E-003  <-- F
 O              10   -2.8425502619309167E-003   -1.1697470556033103E-002   -5.5313994873713995E-003  <-- F
 O              11   -9.2784137852657976E-004    5.9884837413488800E-003   -8.2363531293972309E-003  <-- F
 O              12   -3.9767825021024506E-003    1.1295419807093941E-004   -4.0009915317457817E-003  <-- F
 O              13    8.1051416358948521E-003    6.6118926210009084E-006    1.7508400507173852E-002  <-- F
 O              14   -2.4386263587114678E-002    2.9997985030331485E-002   -1.7839514332069054E-002  <-- F
 O              15    8.2271300796915910E-003    6.6430980060985041E-004    5.3707353331751764E-003  <-- F
 O              16    7.2313242018327466E-003   -2.4584483647550415E-003   -2.0416903774508614E-003  <-- F
 O              17    2.1036975812322060E-002   -1.9186811397383548E-002    3.8726692170235745E-003  <-- F
 O              18    7.0978664222543274E-003   -4.3438398283939606E-003    1.7010836820289379E-002  <-- F
 O              19    2.6497318702203763E-002    1.5248697475942895E-002    4.3230264913404995E-002  <-- F
 O              20    9.5245524128963823E-005    4.8922882980980886E-003   -2.2648381792833437E-003  <-- F
 O              21    3.3363601428497464E-003    6.6149496849471985E-004    1.7850702771195257E-002  <-- F
 O              22    1.2538248052681902E-003   -1.4350310378773858E-002   -6.5610930592049782E-003  <-- F
 O              23   -1.7715984425213906E-002    1.1102497119267161E-002   -1.0579560572261521E-002  <-- F
 O              24    1.6584723888281003E-002   -1.6445733888983384E-002    1.2899483942818715E-002  <-- F
 O              25    1.0899278444798674E-002   -9.5378688240657615E-003   -1.7862078848148589E-003  <-- F
 O              26   -4.0981298008257598E-004   -4.1779411520873827E-003   -1.6445049640114048E-002  <-- F
 O              27    2.4180106245768915E-002   -8.0401980002791987E-003   -5.9028956398248095E-003  <-- F
 O              28    3.2685594085503961E-003    6.6682732651498480E-003   -7.3831323541844046E-003  <-- F
 O              29   -8.6256284685520422E-004    1.0234083832583972E-002   -5.8747998356754334E-003  <-- F
 O              30    4.4916902087720522E-003   -1.0499210436492747E-002   -1.3849001323926011E-002  <-- F
 O              31    6.7642889274642003E-003    2.2917239409616062E-002   -1.7470372492327871E-002  <-- F
 O              32   -2.8746567074570783E-003   -4.0381520769494439E-003    8.4783348733196166E-004  <-- F
 O              33    7.6744986854640805E-003   -1.2135253532220107E-002   -2.0398605724804009E-002  <-- F
 O              34   -3.1040617019928076E-002   -2.2374529653222570E-002   -5.0949992320946061E-002  <-- F
 O              35    1.6499437125825665E-002    1.8024650525411766E-002   -3.5480758271273213E-002  <-- F
 O              36    2.2626151886861377E-002    1.6786190783966391E-002   -2.1291275726579392E-002  <-- F
 O              37   -6.8078213526916469E-004    4.2248456709047455E-003   -1.0339541741013217E-002  <-- F
 O              38   -1.8066598536408088E-003    7.7593698463168303E-003    8.9447543687981995E-003  <-- F
 O              39    5.0103870534854279E-004    3.8640545460843700E-003    6.2175962337550766E-003  <-- F
 O              40    2.8395238544897009E-002    1.6610032364819054E-002    1.9948882361602439E-003  <-- F
 O              41    2.3944117365478888E-002    9.1884619641719525E-003    4.5305174497634465E-003  <-- F
 O              42   -7.2822723262276510E-003    5.0165984629668070E-003   -6.7258329107944131E-003  <-- F
 O              43   -1.7473460898964923E-002   -7.5728197840993850E-003   -1.2163887950149422E-003  <-- F
 O              44   -2.4688134223142529E-002    3.2070481503702029E-002    9.0559497102861292E-003  <-- F
 O              45   -1.9487187567159687E-002   -1.5909653179285062E-002   -9.6482533654433710E-004  <-- F
 O              46   -3.1890132877300764E-003   -1.7550236330879300E-002   -1.9546064277451244E-002  <-- F
 O              47   -2.4660745887879736E-002   -3.3292233023975361E-002    6.2990802320872050E-002  <-- F
 O              48   -1.7364686112650431E-002   -1.1692135539947336E-004    1.0176186876172125E-002  <-- F
 Si              1    6.3340361497085159E-003    1.3550726326852872E-002   -2.0294961611372127E-002  <-- F
 Si              2   -3.6456898663280439E-002   -5.3847702623518352E-002    6.9387203412286432E-003  <-- F
 Si              3    5.1404557292643576E-003   -1.2167601819822507E-002    2.8883697270241848E-003  <-- F
 Si              4   -2.0388255014711275E-002    1.1519350725893533E-003    4.1059660076609147E-002  <-- F
 Si              5   -5.2376149473554786E-003    2.4918769808366707E-002   -5.0467775600885334E-003  <-- F
 Si              6    1.2824096806328482E-002   -5.7446232825309725E-003   -1.2114513345148070E-002  <-- F
 Si              7    1.6883822412106113E-002   -2.6031957148542860E-002    4.4781633091028547E-002  <-- F
 Si              8    6.0101195882645477E-003   -1.7965004930508309E-002    5.4922843469471433E-002  <-- F
 Si              9   -2.6489767099950275E-002    4.7305052798796314E-002    1.8324130904623086E-003  <-- F
 Si             10    1.7377596955562755E-003   -8.3721254301958262E-003    3.2080477485138875E-003  <-- F
 Si             11    1.4701863602160145E-002   -4.3054405102822244E-004    3.1102025277891512E-003  <-- F
 Si             12   -6.6635293112559616E-003   -5.9724464221757078E-003    2.8992090004938556E-003  <-- F
 Si             13    2.1928230442088768E-002    2.0512624995464988E-002    1.8352923284335777E-002  <-- F
 Si             14    1.0419267335537446E-002   -1.9415155330882738E-003   -1.0648708607262445E-002  <-- F
 Si             15   -9.6771103125269081E-003   -1.3905404204310103E-002    1.3214863207351497E-002  <-- F
 Si             16   -2.7170149627454503E-002   -3.9016543450222145E-002   -1.3823660488818015E-002  <-- F
 Si             17   -2.9804624019516743E-003    3.4302035978563038E-002   -3.2181059607815796E-002  <-- F
 Si             18   -4.6184124672672244E-002   -2.4119959452677246E-002    2.5973701922395792E-003  <-- F
 Si             19    3.4116338104415236E-002    7.9619465444897199E-003    4.0731290324822686E-002  <-- F
 Si             20   -1.9952462544296674E-002    2.0457197989827015E-002   -5.0617561911962851E-002  <-- F
 Si             21    8.9651333845522495E-003    4.4968277909280686E-002   -1.1604766150653950E-002  <-- F
 Si             22    2.1464738952520053E-002   -5.1463427296423010E-002   -8.1478783181049698E-004  <-- F
 Si             23   -1.4674827579124526E-002   -7.6464070760818093E-003    3.2463512378119661E-002  <-- F
 Si             24   -9.0055862101436398E-003    1.3242885361889266E-002   -4.1092259124119776E-002  <-- F
""".split("\n")
      self.at = Atoms(self.md_lines, format='md')

   def test_pos(self):
      self.assertArrayAlmostEqual(self.at.pos.T, [[ 0.82415757,  0.15801177,  1.57830446],
                                                [ 2.53076898, -0.93897975, -5.04657556],
                                                [-1.86434147, -2.83961172, -3.4922821 ],
                                                [-1.97972729, -4.06147144,  0.05925173],
                                                [ 1.18847274, -1.05415866, -0.69317345],
                                                [-3.24042651,  1.00616672,  0.29484347],
                                                [ 4.47921167, -4.17483013, -2.91794026],
                                                [ 1.77230813,  3.80021584, -4.91097823],
                                                [ 0.6128812 , -1.84627791, -3.23915048],
                                                [-0.24453807, -4.96851542,  3.20344432],
                                                [ 4.95509673,  4.52596569,  4.89733014],
                                                [ 2.77295249, -2.72206371, -1.72720734],
                                                [-3.29555715, -0.6502742 , -3.56608668],
                                                [ 1.24141085,  1.58083923, -0.59647822],
                                                [-4.02814984,  1.65501566, -4.79038699],
                                                [-5.02542911, -0.30123535,  1.51698466],
                                                [-0.82104588, -2.61339569,  2.8902251 ],
                                                [ 3.30119868,  2.97972056,  3.34635358],
                                                [ 3.76731394, -4.73407339, -0.35015768],
                                                [ 1.28123575,  1.39439288,  3.89592882],
                                                [ 1.77323308, -2.97148745,  3.16604565],
                                                [-0.46666429,  2.31519129,  2.34594658],
                                                [-2.35106204,  3.83436183,  1.6924011 ],
                                                [ 2.97735707,  1.6166049 , -4.39298634],
                                                [-4.67094879,  3.11868702,  3.1303537 ],
                                                [ 0.1744963 , -4.39782472, -3.886851  ],
                                                [-2.82337254, -1.66628727,  4.29971236],
                                                [ 0.49612591, -3.84059728,  1.03418467],
                                                [-0.17221716, -3.52990733, -1.57454344],
                                                [-3.26336528,  1.03224054,  2.89813789],
                                                [ 4.28711251, -4.7232975 ,  2.29406099],
                                                [ 4.88736889,  0.02447551, -3.80440932],
                                                [ 4.35213262, -2.47429481,  0.9367176 ],
                                                [-3.90676634,  1.39580017, -2.16125891],
                                                [ 2.63942269, -0.2979389 , -2.58703566],
                                                [ 0.30296143,  2.88589553, -2.80117025],
                                                [ 3.8509365 , -3.0018909 ,  4.82274868],
                                                [-1.36403181,  1.63746087, -1.37100198],
                                                [-3.1843905 ,  3.41543471, -0.73931858],
                                                [-2.22141544,  3.14043988,  4.10096759],
                                                [-4.48664873, -4.56199961,  0.54991129],
                                                [-0.73573489, -0.15067462,  3.55748229],
                                                [ 3.2456383 ,  0.05188811,  0.4950578 ],
                                                [-2.45297547, -1.13788017,  1.70726949],
                                                [-4.30405609, -3.00489211, -4.2397748 ],
                                                [ 2.13617818,  4.63631957, -2.29926552],
                                                [ 3.50039155, -1.06702839,  2.78128646],
                                                [-0.00516056,  3.87773678, -0.46802139],
                                                [ 2.31004695,  2.43013891,  4.50637493],
                                                [ 4.6796088 , -4.00628008,  0.86077935],
                                                [ 4.92362646, -3.81627983, -4.39276451],
                                                [ 3.30318087,  0.0743956 , -4.06152807],
                                                [-4.23771427,  0.61542241, -3.56974503],
                                                [-0.41345473, -4.31331814, -0.21225319],
                                                [-3.51319426,  0.15784614,  1.54676141],
                                                [-0.34773273, -3.10522915, -3.15445577],
                                                [-1.71539792, -1.39540325,  3.14289297],
                                                [ 1.76449561, -1.52424105, -2.09666064],
                                                [ 0.01821674,  2.47764118, -1.27256839],
                                                [ 4.60872567,  3.94319155,  3.42172254],
                                                [-1.29480847,  3.57436069,  2.83686422],
                                                [-3.04487829,  4.53231134,  0.4207826 ],
                                                [ 0.20775038,  0.94369799,  2.81657942],
                                                [ 0.40247022, -3.56021543,  2.61749676],
                                                [ 4.15669644, -0.96775319,  1.39158793],
                                                [-3.50483328,  2.26869888,  3.87639146],
                                                [-2.94797196,  1.85509143, -1.04826365],
                                                [ 3.3458836 , -4.22215079, -1.75482444],
                                                [ 2.89129051, -2.10107425,  3.93711246],
                                                [ 1.62189845,  0.27561596,  0.18578585],
                                                [ 1.14801827,  4.14005711, -3.48779492],
                                                [-3.03693244, -2.07131968, -4.20340801]])

   def test_vel(self):
      self.assertArrayAlmostEqual(self.at.velo.T, [[  2.18686555e-03,   4.35533564e-03,  -5.57118597e-03],
                                                 [ -4.77449331e-04,   3.19042935e-03,  -4.97809786e-04],
                                                 [ -4.37103491e-03,   2.68725257e-03,  -5.43763245e-03],
                                                 [  2.37453419e-03,  -9.75044217e-04,   2.89668061e-03],
                                                 [ -4.88275321e-03,  -7.50549340e-04,   6.10571255e-03],
                                                 [  4.58671893e-03,   5.19621146e-03,   1.87598419e-03],
                                                 [ -3.39240065e-04,  -5.21216189e-03,   2.90224006e-04],
                                                 [ -1.74146282e-03,  -6.50165534e-03,   2.40965059e-03],
                                                 [  7.27912767e-03,  -4.80926145e-03,   5.68281169e-04],
                                                 [  7.86856417e-04,  -6.81304530e-03,  -9.38915825e-04],
                                                 [  1.15779807e-03,   4.19306147e-04,   2.98363773e-03],
                                                 [  6.54511855e-04,   4.22305402e-03,  -1.55623494e-03],
                                                 [  1.40517225e-03,   7.36256922e-03,  -1.06317610e-03],
                                                 [ -2.38608163e-03,   2.88000127e-03,  -6.44125858e-03],
                                                 [ -7.74096045e-04,  -1.66107482e-03,   5.39354198e-04],
                                                 [ -2.17800787e-04,  -2.46589649e-03,   3.76932224e-03],
                                                 [ -3.17188115e-03,  -1.83156804e-03,   6.91481211e-03],
                                                 [ -4.93652463e-03,   7.06460769e-04,  -6.73378433e-04],
                                                 [ -7.54251586e-03,   1.38309554e-03,   4.60630762e-03],
                                                 [  7.70741354e-03,  -1.37917055e-03,  -5.91281495e-03],
                                                 [ -1.06005219e-04,  -3.56477429e-03,  -3.15373886e-03],
                                                 [ -4.10103449e-03,  -2.70157783e-03,  -5.48278499e-03],
                                                 [  3.59152186e-03,  -1.38602240e-03,  -4.77009259e-03],
                                                 [  6.92997413e-05,   9.93286545e-04,  -4.30067177e-04],
                                                 [ -2.01377959e-03,   5.73280482e-03,  -7.42927927e-04],
                                                 [ -9.65020860e-04,   4.62139883e-03,  -8.01625220e-03],
                                                 [  6.29155197e-04,  -2.44604326e-03,   1.79754581e-03],
                                                 [  1.56151495e-03,  -5.98210329e-03,  -5.61317928e-03],
                                                 [ -3.33919033e-04,  -1.80175454e-03,   5.51582909e-03],
                                                 [ -8.57860104e-04,  -3.50305603e-03,   5.67266144e-03],
                                                 [  3.01769589e-03,   1.84246596e-03,  -3.23450037e-03],
                                                 [  2.69190101e-03,   6.30774059e-03,  -2.75937372e-03],
                                                 [  3.82265557e-04,  -9.51407577e-03,  -2.74619314e-03],
                                                 [  4.38578846e-03,  -1.21482073e-03,  -2.20262370e-03],
                                                 [  1.53092180e-03,  -1.61959773e-03,   8.69267323e-04],
                                                 [  2.33847250e-03,  -2.59157439e-03,  -1.66678543e-03],
                                                 [  1.64314407e-02,  -2.44162432e-03,  -2.07839438e-03],
                                                 [  5.39238985e-04,   1.29776916e-03,  -3.54313648e-03],
                                                 [ -5.29611429e-03,  -1.22823285e-03,  -4.64563000e-04],
                                                 [  6.74616050e-03,  -7.37524040e-03,   1.12707156e-02],
                                                 [ -5.83899936e-03,  -5.08752047e-04,   4.44883274e-03],
                                                 [ -5.07219757e-03,  -6.86331822e-03,  -4.12259278e-03],
                                                 [ -8.76545830e-05,   4.64176601e-04,   5.45735232e-03],
                                                 [  2.92127935e-03,  -6.16596923e-04,  -5.19647556e-03],
                                                 [ -2.19173974e-03,  -2.62961013e-03,   1.31603601e-03],
                                                 [ -2.95576408e-03,  -6.15715750e-03,  -6.52947982e-03],
                                                 [ -2.69412664e-03,  -2.41105292e-03,  -3.18619359e-03],
                                                 [ -5.67028060e-03,  -1.20003465e-03,   3.83318681e-03],
                                                 [  3.52458237e-04,  -2.89534046e-03,   2.61737730e-03],
                                                 [ -5.10449331e-03,   3.87190080e-03,   4.70545991e-03],
                                                 [  1.03458166e-04,   3.49360025e-03,  -1.03825374e-03],
                                                 [  2.90243398e-04,   1.66433355e-03,   3.20570897e-04],
                                                 [ -3.23442024e-03,   1.70036205e-03,   2.60104394e-03],
                                                 [ -2.67448599e-03,  -1.85733758e-03,   1.62803949e-03],
                                                 [  7.58863059e-04,  -1.10893263e-03,   5.86028469e-06],
                                                 [  1.27545116e-03,   2.66024694e-03,  -5.36114266e-04],
                                                 [ -3.55683442e-04,  -5.59239051e-03,   4.94187076e-03],
                                                 [  2.23794153e-04,  -2.13735333e-03,   1.23211127e-03],
                                                 [  3.24515001e-03,   1.49839707e-04,   2.01164028e-03],
                                                 [  2.62169603e-03,  -6.47285198e-04,  -1.64814377e-03],
                                                 [  2.21885103e-04,   1.76480955e-03,   1.89975157e-03],
                                                 [ -4.05203052e-03,  -1.45360425e-03,   2.00261391e-03],
                                                 [ -1.68718608e-04,   4.86868145e-03,   9.70291766e-04],
                                                 [  2.49190811e-03,   8.64237637e-04,   4.36232561e-03],
                                                 [  9.28233529e-04,   4.45584588e-03,  -2.59030461e-03],
                                                 [  9.46198366e-04,   3.05146392e-03,   2.18200648e-04],
                                                 [  1.04809183e-05,   3.03447860e-03,  -2.78185132e-03],
                                                 [  4.29665815e-04,   1.00547254e-03,  -2.83733613e-03],
                                                 [  6.41198154e-03,   5.43017034e-03,   7.05844156e-04],
                                                 [ -2.20834512e-03,   5.37082425e-06,  -2.72615726e-03],
                                                 [ -3.54893383e-03,  -2.23336935e-04,  -7.41453748e-04],
                                                 [ -2.35406025e-03,   4.38037441e-03,  -3.42282077e-03]])

   def test_force(self):
      self.assertArrayAlmostEqual(self.at.force.T, [[  3.39025894e-01,   2.82725740e-01,  -4.19079393e-01],
                                                  [  2.18514501e-01,  -1.81572503e-01,  -2.94385317e-01],
                                                  [  2.50708323e-01,   1.31137662e-01,  -1.87007313e-01],
                                                  [ -4.32671101e-01,  -5.53770708e-01,   6.46612048e-02],
                                                  [ -5.76939408e-01,   1.21849719e+00,   1.05803989e+00],
                                                  [  5.45397933e-01,   1.08658310e+00,  -1.94782050e+00],
                                                  [ -2.70238495e-01,  -1.31613013e-01,   1.10709253e+00],
                                                  [ -3.56951121e-02,  -6.56988666e-01,   2.45438059e-01],
                                                  [  8.52251457e-01,   4.47516970e-01,   3.43059173e-01],
                                                  [ -1.46169854e-01,  -6.01508295e-01,  -2.84436080e-01],
                                                  [ -4.77115358e-02,   3.07940304e-01,  -4.23530429e-01],
                                                  [ -2.04494437e-01,   5.80834008e-03,  -2.05739316e-01],
                                                  [  4.16783261e-01,   3.39997287e-04,   9.00318413e-01],
                                                  [ -1.25399246e+00,   1.54255886e+00,  -9.17344976e-01],
                                                  [  4.23056161e-01,   3.41601933e-02,   2.76174395e-01],
                                                  [  3.71849749e-01,  -1.26418534e-01,  -1.04987971e-01],
                                                  [  1.08176510e+00,  -9.86625797e-01,   1.99140716e-01],
                                                  [  3.64987072e-01,  -2.23369290e-01,   8.74732653e-01],
                                                  [  1.36254731e+00,   7.84119778e-01,   2.22299024e+00],
                                                  [  4.89772319e-03,   2.51571652e-01,  -1.16462695e-01],
                                                  [  1.71562586e-01,   3.40154488e-02,   9.17920309e-01],
                                                  [  6.44742824e-02,  -7.37922843e-01,  -3.37385068e-01],
                                                  [ -9.10992811e-01,   5.70913522e-01,  -5.44023035e-01],
                                                  [  8.52821038e-01,  -8.45673883e-01,   6.63318326e-01],
                                                  [  5.60463594e-01,  -4.90457077e-01,  -9.18505290e-02],
                                                  [ -2.10734368e-02,  -2.14838434e-01,  -8.45638698e-01],
                                                  [  1.24339142e+00,  -4.13443724e-01,  -3.03539186e-01],
                                                  [  1.68076131e-01,   3.42896497e-01,  -3.79656040e-01],
                                                  [ -4.43547778e-02,   5.26257902e-01,  -3.02094441e-01],
                                                  [  2.30972064e-01,  -5.39891264e-01,  -7.12144487e-01],
                                                  [  3.47833823e-01,   1.17845217e+00,  -8.98362934e-01],
                                                  [ -1.47820834e-01,  -2.07650189e-01,   4.35973634e-02],
                                                  [  3.94638704e-01,  -6.24020007e-01,  -1.04893878e+00],
                                                  [ -1.59617317e+00,  -1.15054491e+00,  -2.61995470e+00],
                                                  [  8.48435415e-01,   9.26865064e-01,  -1.82449448e+00],
                                                  [  1.16348385e+00,   8.63180886e-01,  -1.09484173e+00],
                                                  [ -3.50072351e-02,   2.17250362e-01,  -5.31680767e-01],
                                                  [ -9.29022119e-02,   3.99002956e-01,   4.59957896e-01],
                                                  [  2.57644536e-02,   1.98697731e-01,   3.19721746e-01],
                                                  [  1.46014230e+00,   8.54122453e-01,   1.02581307e-01],
                                                  [  1.23125638e+00,   4.72489848e-01,   2.32968642e-01],
                                                  [ -3.74469608e-01,   2.57963940e-01,  -3.45856334e-01],
                                                  [ -8.98521746e-01,  -3.89410163e-01,  -6.25492448e-02],
                                                  [ -1.26951527e+00,   1.64913094e+00,   4.65675792e-01],
                                                  [ -1.00207176e+00,  -8.18107497e-01,  -4.96133280e-02],
                                                  [ -1.63985704e-01,  -9.02469699e-01,  -1.00509933e+00],
                                                  [ -1.26810691e+00,  -1.71195595e+00,   3.23911823e+00],
                                                  [ -8.92928320e-01,  -6.01233957e-03,   5.23280720e-01],
                                                  [  3.25709329e-01,   6.96806566e-01,  -1.04360919e+00],
                                                  [ -1.87468965e+00,  -2.76896100e+00,   3.56803449e-01],
                                                  [  2.64332938e-01,  -6.25683423e-01,   1.48525986e-01],
                                                  [ -1.04840653e+00,   5.92349002e-02,   2.11137322e+00],
                                                  [ -2.69329067e-01,   1.28137503e+00,  -2.59515812e-01],
                                                  [  6.59441763e-01,  -2.95400492e-01,  -6.22953503e-01],
                                                  [  8.68201307e-01,  -1.33861745e+00,   2.30276483e+00],
                                                  [  3.09052865e-01,  -9.23797964e-01,   2.82424698e+00],
                                                  [ -1.36215899e+00,   2.43252432e+00,   9.42264969e-02],
                                                  [  8.93592222e-02,  -4.30512124e-01,   1.64964496e-01],
                                                  [  7.56000442e-01,  -2.21394717e-02,   1.59933091e-01],
                                                  [ -3.42652553e-01,  -3.07115632e-01,   1.49083364e-01],
                                                  [  1.12759527e+00,   1.05480189e+00,   9.43745533e-01],
                                                  [  5.35780423e-01,  -9.98367717e-02,  -5.47578771e-01],
                                                  [ -4.97617164e-01,  -7.15044841e-01,   6.79535785e-01],
                                                  [ -1.39714567e+00,  -2.00631191e+00,  -7.10841409e-01],
                                                  [ -1.53261583e-01,   1.76388212e+00,  -1.65481710e+00],
                                                  [ -2.37488386e+00,  -1.24029854e+00,   1.33562184e-01],
                                                  [  1.75433315e+00,   4.09419871e-01,   2.09448777e+00],
                                                  [ -1.02599717e+00,   1.05195172e+00,  -2.60286044e+00],
                                                  [  4.61005828e-01,   2.31236249e+00,  -5.96741241e-01],
                                                  [  1.10376158e+00,  -2.64635660e+00,  -4.18980870e-02],
                                                  [ -7.54610193e-01,  -3.93194175e-01,   1.66934141e+00],
                                                  [ -4.63085996e-01,   6.80976742e-01,  -2.11304953e+00]])
      
      
class TestCastepOutput(QuippyTestCase):

   def setUp(self):
      castep_output = """+-------------------------------------------------+
 |                                                 |
 |      CCC   AA    SSS  TTTTT  EEEEE  PPPP        |
 |     C     A  A  S       T    E      P   P       |
 |     C     AAAA   SS     T    EEE    PPPP        |
 |     C     A  A     S    T    E      P           |
 |      CCC  A  A  SSS     T    EEEEE  P           |
 |                                                 |
 +-------------------------------------------------+
 |                                                 |
 | Welcome to Academic Release CASTEP version 4.3  |
 | Ab Initio Total Energy Program                  |
 |                                                 |
 | Authors:                                        |
 | M. Segall, M. Probert, C. Pickard, P. Hasnip,   |
 | S. Clark, K. Refson, M. Payne                   |
 |                                                 |
 | Contributors:                                   |
 | P. Lindan, P. Haynes, J. White, V. Milman,      |
 | N. Govind, M. Gibson, P. Tulip, V. Cocula,      |
 | B. Montanari, D. Quigley, M. Glover,            |
 | L. Bernasconi, A. Perlov, M. Plummer            |
 |                                                 |
 | Copyright (c) 2000 - 2008                       |
 |                                                 |
 |          Internal TCM version                   |
 |                                                 |
 |             11th April 2008                     |
 |                                                 |
 | Please cite                                     |
 |                                                 |
 |     "First principles methods using CASTEP"     |
 |                                                 |
 |         Zeitschrift fuer Krystallographie       |
 |           220(5-6) pp. 567-570 (2005)           |
 |                                                 |
 | S. J. Clark, M. D. Segall, C. J. Pickard,       |
 | P. J. Hasnip, M. J. Probert, K. Refson,         |
 | M. C. Payne                                     |
 |                                                 |
 |       in all publications arising from          |
 |              your use of CASTEP                 |
 |                                                 |
 +-------------------------------------------------+
 
 
 Pseudo atomic calculation performed for Si 3s2 3p2
 
 Converged in 15 iterations to a total energy of -101.8974 eV
 
Calculation parallelised over   32 nodes.
K-points are distributed over    8 groups, each containing    4 nodes.

 ************************************ Title ************************************
 

 ***************************** General Parameters ******************************
  
 output verbosity                               : normal  (1)
 write checkpoint data to                       : si_unrec.check
 type of calculation                            : geometry optimization
 stress calculation                             : off
 density difference calculation                 : off
 electron localisation func (ELF) calculation   : off
 unlimited duration calculation
 timing information                             : on
 memory usage estimate                          : on
 write final potential to formatted file        : off
 write final density to formatted file          : off
  
 output         length unit                     : A
 output           mass unit                     : amu
 output           time unit                     : ps
 output         charge unit                     : e
 output         energy unit                     : eV
 output          force unit                     : eV/A
 output       velocity unit                     : A/ps
 output       pressure unit                     : GPa
 output     inv_length unit                     : 1/A
 output      frequency unit                     : cm-1
 output force constant unit                     : eV/A**2
 output         volume unit                     : A**3
 output   IR intensity unit                     : (D/A)**2/amu
 output         dipole unit                     : D
 output         efield unit                     : eV/A/e
  
 wavefunctions paging                           : none
 data distribution                              : optimal for this architecture
 optimization strategy                          : maximize speed(+++)

 *********************** Exchange-Correlation Parameters ***********************
  
 using functional                               : Perdew Burke Ernzerhof
 Divergence correction                          : off

 ************************* Pseudopotential Parameters **************************
  
 pseudopotential representation                 : reciprocal space
 <beta|phi> representation                      : reciprocal space

 **************************** Basis Set Parameters *****************************
  
 plane wave basis set cut-off                   :   250.0000   eV
 finite basis set correction                    : none

 **************************** Electronic Parameters ****************************
  
 number of  electrons                           :  112.0    
 net charge of system                           :  0.000    
 net spin   of system                           :  0.000    
 number of  up  spins                           :  56.00    
 number of down spins                           :  56.00    
 treating system as non-spin-polarized
 number of bands                                :         67
 electron temperature                           :     0.0000   K

 ********************* Electronic Minimization Parameters **********************
  
 Method: Treating system as metallic with ensemble DFT treatment of electrons,
         and number of  SD  steps               :          1
         and number of  CG  steps               :          4
  
 total energy / atom convergence tol.           : 0.1000E-04   eV
 max force / atom convergence tol.              : ignored
 convergence tolerance window                   :          3   cycles
 max. number of SCF cycles                      :        200
 number of fixed-spin iterations                :         10
 smearing scheme                                : Gaussian
 smearing width                                 : 0.2000       eV
 Fermi energy convergence tolerance             : 0.1000E-06   eV
 number of occupancy cycles                     :          6

 ********************** Geometry Optimization Parameters ***********************
  
 optimization method                            : BFGS
 variable cell method                           : fixed basis quality
 max. number of steps                           :         30
 estimated bulk modulus                         :  500.0       GPa
 estimated <frequency>                          :  1668.       cm-1
 line minimiser tolerance                       :     0.4000
 total energy convergence tolerance             : 0.2000E-04   eV/atom
        max ionic |force| tolerance             : 0.1000E-01   eV/A
 max ionic |displacement| tolerance             : 0.1000E-02   A
   max |stress component| tolerance             : 0.1000       GPa
 convergence tolerance window                   :          2   steps
 backup results every                           :          5   steps

 *******************************************************************************
  
 
                           -------------------------------
                                      Unit Cell
                           -------------------------------
        Real Lattice(A)                      Reciprocal Lattice(1/A)
   6.5166950   0.0000000   0.0000000        0.9641675   0.0000000   0.0000000
   0.0000000   3.7624150   0.0000000        0.0000000   1.6699873   0.0000000
   0.0000000   0.0000000  41.4832560        0.0000000   0.0000000   0.1514632
 
                       Lattice parameters(A)       Cell Angles
                    a =    6.516695          alpha =   90.000000
                    b =    3.762415          beta  =   90.000000
                    c =   41.483256          gamma =   90.000000
 
                       Current cell volume = 1017.107669       A**3
 
                           -------------------------------
                                     Cell Contents
                           -------------------------------
                         Total number of ions in cell =   28
                      Total number of species in cell =    1
                        Max number of any one species =   28
 
            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            x  Element    Atom        Fractional coordinates of atoms  x
            x            Number           u          v          w      x
            x----------------------------------------------------------x
            x  Si           1         0.106540   0.250079   0.497848   x
            x  Si           2         0.579755   0.749871   0.498543   x
            x  Si           3         0.756707   0.250000   0.425444   x
            x  Si           4         0.254663   0.750000   0.421220   x
            x  Si           5         0.416667   0.250000   0.348256   x
            x  Si           6         0.916667   0.750000   0.348256   x
            x  Si           7         0.083333   0.250000   0.274202   x
            x  Si           8         0.583333   0.750000   0.274202   x
            x  Si           9         0.750000   0.250000   0.200148   x
            x  Si          10         0.250000   0.750000   0.200148   x
            x  Si          11         0.419863   0.250000   0.123974   x
            x  Si          12         0.898267   0.750000   0.124637   x
            x  Si          13         0.584289   0.727002   0.045298   x
            x  Si          14         0.090225   0.248917   0.049714   x
            x  Si          15         0.911993   0.753152   0.518699   x
            x  Si          16         0.425893   0.238906   0.522252   x
            x  Si          17         0.101733   0.250000   0.442281   x
            x  Si          18         0.580136   0.750000   0.442944   x
            x  Si          19         0.750000   0.250000   0.366770   x
            x  Si          20         0.250000   0.750000   0.366770   x
            x  Si          21         0.416667   0.250000   0.292716   x
            x  Si          22         0.916667   0.750000   0.292716   x
            x  Si          23         0.083333   0.250000   0.218662   x
            x  Si          24         0.583333   0.750000   0.218662   x
            x  Si          25         0.745337   0.250000   0.145698   x
            x  Si          26         0.243293   0.750000   0.141474   x
            x  Si          27         0.426070   0.253751   0.068644   x
            x  Si          28         0.900053   0.747332   0.069504   x
            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 
                         No user defined ionic velocities
 
                           -------------------------------
                                   Details of Species
                           -------------------------------
 
                               Mass of species in AMU
                                    Si   28.0855000
 
                          Electric Quadrupole Moment (Barn)
                                    Si    1.0000000 No Isotope Defined
 
                          Files used for pseudopotentials:
                                    Si Si_00PBE.usp
 
                           -------------------------------
                              k-Points For BZ Sampling
                           -------------------------------
                       MP grid size for SCF calculation is  4  4  1
                         Number of kpoints used =           8
 
             +++++++++++++++++++++++++++++++++++++++++++++++++++++++
             +  Number       Fractional coordinates        Weight  +
             +-----------------------------------------------------+
             +     1  -0.375000   0.375000   0.000000   0.1250000  +
             +     2   0.375000   0.375000   0.000000   0.1250000  +
             +     3   0.375000   0.125000   0.000000   0.1250000  +
             +     4  -0.375000   0.125000   0.000000   0.1250000  +
             +     5   0.125000   0.375000   0.000000   0.1250000  +
             +     6  -0.125000   0.375000   0.000000   0.1250000  +
             +     7  -0.125000   0.125000   0.000000   0.1250000  +
             +     8   0.125000   0.125000   0.000000   0.1250000  +
             +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
                           -------------------------------
                               Symmetry and Constraints
                           -------------------------------
 
         There are no symmetry operations specified or generated for this cell
                      Number of ionic constraints     =           3
 
                          Centre of mass is constrained
             Set iprint > 1 for details of linear ionic constraints
 
                         Number of cell constraints= 0
                         Cell constraints are:  1 2 3 4 5 6
 
                         External pressure/stress (GPa)
                          0.00000   0.00000   0.00000
                                    0.00000   0.00000
                                              0.00000
  
+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER NODE -----------------+
|                                                     Memory          Disk    |
| Model and support data                               34.2 MB         0.0 MB |
| Electronic energy minimisation requirements          17.6 MB         0.0 MB |
|                                               ----------------------------- |
| Approx. total storage required per node              51.8 MB         0.0 MB |
|                                                                             |
| Requirements will fluctuate during execution and may exceed these estimates |
+-----------------------------------------------------------------------------+
------------------------------------------------------------------------ <-- SCF
SCF loop      Energy           Fermi           Energy gain       Timer   <-- SCF
                               energy          per atom          (sec)   <-- SCF
------------------------------------------------------------------------ <-- SCF
Initial   3.54031177E+004  2.07624515E+002                         9.81  <-- SCF
      1   1.24048433E+004  5.15608121E+001   8.21366942E+002      78.47  <-- SCF
      2   3.89718601E+002  6.29371518E+000   4.29111596E+002     150.50  <-- SCF
      3  -2.50760290E+003  3.86196194E+000   1.03475768E+002     255.60  <-- SCF
      4  -2.97887188E+003  2.42654645E+000   1.68310352E+001     357.72  <-- SCF
      5  -2.99820130E+003  1.43462681E+000   6.90336528E-001     445.54  <-- SCF
      6  -2.99875468E+003  1.30448416E+000   1.97634816E-002     557.25  <-- SCF
      7  -2.99888455E+003  1.21689643E+000   4.63819249E-003     662.47  <-- SCF
      8  -2.99893812E+003  1.15021819E+000   1.91303282E-003     767.59  <-- SCF
      9  -2.99894517E+003  1.13675715E+000   2.52067323E-004     869.06  <-- SCF
     10  -2.99894766E+003  1.12484872E+000   8.88423503E-005     974.20  <-- SCF
     11  -2.99894868E+003  1.11660756E+000   3.64655975E-005    1079.95  <-- SCF
     12  -2.99894908E+003  1.11115264E+000   1.40302298E-005    1186.48  <-- SCF
     13  -2.99894925E+003  1.10459834E+000   6.35042893E-006    1270.42  <-- SCF
     14  -2.99894926E+003  1.10453910E+000   3.04612412E-007    1328.27  <-- SCF
------------------------------------------------------------------------ <-- SCF
 
Final energy, E             =  -2998.588093584     eV
Final free energy (E-TS)    =  -2998.949262410     eV
(energies not corrected for finite basis set)
 
NB est. 0K energy (E-0.5TS)      =  -2998.768677997     eV
 

WARNING - doing cell minimisation without finite basis correction:
        => now searching for a stationary and not minimum energy
 
+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER NODE -----------------+
|                                                     Memory          Disk    |
| Model and support data                               34.2 MB         0.0 MB |
| Geometry minimisation requirements                   20.4 MB         0.0 MB |
|                                               ----------------------------- |
| Approx. total storage required per node              54.6 MB         0.0 MB |
|                                                                             |
| Requirements will fluctuate during execution and may exceed these estimates |
+-----------------------------------------------------------------------------+
 
 
 ******************************** Forces *********************************
 *                                                                       *
 *                      Cartesian components (eV/A)                      *
 * --------------------------------------------------------------------- *
 *              x                    y                    z              *
 *                                                                       *
 * Si   1     -0.94623              0.09706              1.35555         *
 * Si   2     -0.07382             -0.34057              0.94944         *
 * Si   3      0.12301              0.01579             -0.50095         *
 * Si   4      0.12532              0.00476              0.78360         *
 * Si   5      0.06472              0.00354              0.01463         *
 * Si   6     -0.03845              0.00409              0.07130         *
 * Si   7      0.00513              0.00056             -0.06598         *
 * Si   8      0.03838              0.00244             -0.02917         *
 * Si   9      0.01983              0.00106              0.52932         *
 * Si  10      0.01505             -0.00817             -0.81851         *
 * Si  11     -0.52546             -0.06765              0.37344         *
 * Si  12      0.87906             -0.00659              0.40109         *
 * Si  13     -0.29319              0.86885              0.33083         *
 * Si  14     -0.82913             -0.15815              0.89809         *
 * Si  15      1.06981             -0.13557             -1.01460         *
 * Si  16      0.28528              0.35565             -0.48995         *
 * Si  17     -0.85544              0.00706             -0.08505         *
 * Si  18      0.37325             -0.02508             -0.25474         *
 * Si  19     -0.01886             -0.00121              0.97372         *
 * Si  20     -0.03053              0.00493             -0.39012         *
 * Si  21     -0.02792              0.00423             -0.07803         *
 * Si  22     -0.00190              0.00328             -0.06052         *
 * Si  23      0.08379              0.00762              0.01907         *
 * Si  24     -0.02017              0.00623              0.08109         *
 * Si  25     -0.25471              0.02076             -0.69736         *
 * Si  26     -0.25003              0.04108              0.42654         *
 * Si  27      0.20315             -0.78916             -1.24898         *
 * Si  28      0.88007              0.08315             -1.47375         *
 *                                                                       *
 *************************************************************************
 
 ***************** Stress Tensor *****************
 *                                               *
 *          Cartesian components (GPa)           *
 * --------------------------------------------- *
 *             x             y             z     *
 *                                               *
 *  x     -3.164402     -0.105182     -0.217689  *
 *  y     -0.105182     -2.411314      0.075166  *
 *  z     -0.217689      0.075166     -2.604180  *
 *                                               *
 *  Pressure:    2.7266                          *
 *                                               *
 *************************************************
 BFGS: finished iteration     0 with enthalpy= -2.99894926E+003 eV
  
 +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
 | Parameter |      value      |    tolerance    |    units   | OK? | <-- BFGS
 +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
 |  dE/ion   |   0.000000E+000 |   2.000000E-005 |         eV | No  | <-- BFGS
 |  |F|max   |   1.718538E+000 |   1.000000E-002 |       eV/A | No  | <-- BFGS
 |  |dR|max  |   0.000000E+000 |   1.000000E-003 |          A | No  | <-- BFGS
 |   Smax    |   3.164402E+000 |   1.000000E-001 |        GPa | No  | <-- BFGS
 +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
  

================================================================================
 Starting BFGS iteration          1 ...
================================================================================
  
 +------------+-------------+-------------+-----------------+ <-- min BFGS
 |    Step    |   lambda    |   F.delta   |    enthalpy     | <-- min BFGS
 +------------+-------------+-------------+-----------------+ <-- min BFGS
 |  previous  |    0.000000 |    0.006113 |    -2998.949262 | <-- min BFGS
 +------------+-------------+-------------+-----------------+ <-- min BFGS
  

--------------------------------------------------------------------------------
 BFGS: starting iteration         1 with trial guess (lambda=  1.000000)
--------------------------------------------------------------------------------
 
+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER NODE -----------------+
|                                                     Memory          Disk    |
| Model and support data                               34.3 MB         0.0 MB |
| Geometry minimisation requirements                   20.6 MB         0.0 MB |
|                                               ----------------------------- |
| Approx. total storage required per node              54.9 MB         0.0 MB |
|                                                                             |
| Requirements will fluctuate during execution and may exceed these estimates |
+-----------------------------------------------------------------------------+
 
 
                           -------------------------------
                                      Unit Cell
                           -------------------------------
        Real Lattice(A)              Reciprocal Lattice(1/A)
   6.5304426   0.0004570   0.0009457        0.9621378  -0.0000674  -0.0001394
   0.0002638   3.7684632  -0.0001885       -0.0001167   1.6673071   0.0000834
   0.0060203  -0.0020788  41.5552759       -0.0000219   0.0000076   0.1512007
 
                       Lattice parameters(A)       Cell Angles
                    a =    6.530443          alpha =   90.005732
                    b =    3.768463          beta  =   89.983402
                    c =   41.555276          gamma =   89.991980
 
                Current cell volume = 1022.664217 A**3
 
                           -------------------------------
                                     Cell Contents
                           -------------------------------
 
            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            x  Element    Atom        Fractional coordinates of atoms  x
            x            Number           u          v          w      x
            x----------------------------------------------------------x
            x  Si           1         0.106034   0.250169   0.497962   x
            x  Si           2         0.579715   0.749556   0.498623   x
            x  Si           3         0.756773   0.250015   0.425402   x
            x  Si           4         0.254730   0.750004   0.421285   x
            x  Si           5         0.416701   0.250003   0.348257   x
            x  Si           6         0.916646   0.750004   0.348262   x
            x  Si           7         0.083336   0.250001   0.274197   x
            x  Si           8         0.583354   0.750002   0.274200   x
            x  Si           9         0.750010   0.250001   0.200193   x
            x  Si          10         0.250008   0.749992   0.200080   x
            x  Si          11         0.419583   0.249937   0.124006   x
            x  Si          12         0.898737   0.749994   0.124671   x
            x  Si          13         0.584132   0.727805   0.045326   x
            x  Si          14         0.089782   0.248770   0.049789   x
            x  Si          15         0.912564   0.753027   0.518613   x
            x  Si          16         0.426046   0.239235   0.522211   x
            x  Si          17         0.101276   0.250007   0.442274   x
            x  Si          18         0.580336   0.749977   0.442922   x
            x  Si          19         0.749990   0.249999   0.366851   x
            x  Si          20         0.249984   0.750004   0.366737   x
            x  Si          21         0.416652   0.250004   0.292709   x
            x  Si          22         0.916665   0.750003   0.292711   x
            x  Si          23         0.083378   0.250007   0.218663   x
            x  Si          24         0.583323   0.750006   0.218669   x
            x  Si          25         0.745201   0.250019   0.145640   x
            x  Si          26         0.243159   0.750038   0.141510   x
            x  Si          27         0.426179   0.253021   0.068540   x
            x  Si          28         0.900523   0.747409   0.069380   x
            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 
------------------------------------------------------------------------ <-- SCF
SCF loop      Energy           Fermi           Energy gain       Timer   <-- SCF
                               energy          per atom          (sec)   <-- SCF
------------------------------------------------------------------------ <-- SCF
Initial  -2.99538389E+003  1.02887254E+000                      1358.14  <-- SCF
      1  -2.99904364E+003  1.05023026E+000   1.30705257E-001    1451.33  <-- SCF
      2  -2.99911098E+003  1.07864734E+000   2.40506967E-003    1558.74  <-- SCF
      3  -2.99911383E+003  1.08246895E+000   1.01730832E-004    1667.56  <-- SCF
      4  -2.99911392E+003  1.08362147E+000   3.09905433E-006    1744.87  <-- SCF
      5  -2.99911392E+003  1.08362147E+000   4.48215836E-008    1793.75  <-- SCF
------------------------------------------------------------------------ <-- SCF
 
Final energy, E             =  -2998.752694423     eV
Final free energy (E-TS)    =  -2999.113917987     eV
(energies not corrected for finite basis set)
 
NB est. 0K energy (E-0.5TS)      =  -2998.933306205     eV
 
 
 ******************************** Forces *********************************
 *                                                                       *
 *                      Cartesian components (eV/A)                      *
 * --------------------------------------------------------------------- *
 *              x                    y                    z              *
 *                                                                       *
 * Si   1     -0.91129              0.08933              1.33909         *
 * Si   2     -0.04461             -0.32135              0.92970         *
 * Si   3      0.10609              0.01266             -0.46687         *
 * Si   4      0.14473              0.00379              0.74776         *
 * Si   5      0.07651              0.00260              0.01429         *
 * Si   6     -0.04214              0.00260              0.07650         *
 * Si   7      0.01449             -0.00007             -0.06360         *
 * Si   8      0.04083              0.00158             -0.02571         *
 * Si   9      0.02826              0.00013              0.46589         *
 * Si  10      0.01954             -0.00810             -0.74676         *
 * Si  11     -0.50692             -0.05839              0.31912         *
 * Si  12      0.83695             -0.00573              0.33713         *
 * Si  13     -0.27067              0.82061              0.38956         *
 * Si  14     -0.78760             -0.13883              0.92323         *
 * Si  15      1.02752             -0.14112             -1.05804         *
 * Si  16      0.25516              0.33139             -0.53599         *
 * Si  17     -0.82092              0.00953             -0.05532         *
 * Si  18      0.36652             -0.01923             -0.21166         *
 * Si  19     -0.01887              0.00255              0.89650         *
 * Si  20     -0.03440              0.00783             -0.33759         *
 * Si  21     -0.02753              0.00702             -0.07381         *
 * Si  22     -0.00720              0.00626             -0.05727         *
 * Si  23      0.08791              0.01012              0.00511         *
 * Si  24     -0.03706              0.00928              0.07246         *
 * Si  25     -0.26411              0.02162             -0.64766         *
 * Si  26     -0.22117              0.03740              0.39430         *
 * Si  27      0.18498             -0.74557             -1.19508         *
 * Si  28      0.80501              0.06209             -1.43527         *
 *                                                                       *
 *************************************************************************
 
 ***************** Stress Tensor *****************
 *                                               *
 *          Cartesian components (GPa)           *
 * --------------------------------------------- *
 *             x             y             z     *
 *                                               *
 *  x     -2.867761     -0.089309     -0.177955  *
 *  y     -0.089309     -2.165759      0.066595  *
 *  z     -0.177955      0.066595     -2.252184  *
 *                                               *
 *  Pressure:    2.4286                          *
 *                                               *
 *************************************************
  
 +------------+-------------+-------------+-----------------+ <-- min BFGS
 |    Step    |   lambda    |   F.delta   |    enthalpy     | <-- min BFGS
 +------------+-------------+-------------+-----------------+ <-- min BFGS
 |  previous  |    0.000000 |    0.006113 |    -2998.949262 | <-- min BFGS
 | trial step |    1.000000 |    0.005645 |    -2999.113918 | <-- min BFGS
 +------------+-------------+-------------+-----------------+ <-- min BFGS
  

--------------------------------------------------------------------------------
 BFGS: improving iteration         1 with line minimization (lambda=  3.659543)
--------------------------------------------------------------------------------
 
+---------------- MEMORY AND SCRATCH DISK ESTIMATES PER NODE -----------------+
|                                                     Memory          Disk    |
| Model and support data                               34.4 MB         0.0 MB |
| Geometry minimisation requirements                   20.7 MB         0.0 MB |
|                                               ----------------------------- |
| Approx. total storage required per node              55.1 MB         0.0 MB |
|                                                                             |
| Requirements will fluctuate during execution and may exceed these estimates |
+-----------------------------------------------------------------------------+
 
 
                           -------------------------------
                                      Unit Cell
                           -------------------------------
        Real Lattice(A)              Reciprocal Lattice(1/A)
   6.5670050   0.0016723   0.0034610        0.9567813  -0.0002442  -0.0005050
   0.0009655   3.7845488  -0.0006900       -0.0004229   1.6602206   0.0003028
   0.0220315  -0.0076073  41.7468160       -0.0000793   0.0000275   0.1505070
 
                       Lattice parameters(A)       Cell Angles
                    a =    6.567006          alpha =   90.020879
                    b =    3.784549          beta  =   89.939569
                    c =   41.746822          gamma =   89.970799
 
                Current cell volume = 1037.539534 A**3
 
                           -------------------------------
                                     Cell Contents
                           -------------------------------
 
            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            x  Element    Atom        Fractional coordinates of atoms  x
            x            Number           u          v          w      x
            x----------------------------------------------------------x
            x  Si           1         0.104690   0.250408   0.498265   x
            x  Si           2         0.579610   0.748718   0.498835   x
            x  Si           3         0.756948   0.250054   0.425290   x
            x  Si           4         0.254908   0.750016   0.421460   x
            x  Si           5         0.416793   0.250012   0.348261   x
            x  Si           6         0.916591   0.750014   0.348278   x
            x  Si           7         0.083343   0.250002   0.274182   x
            x  Si           8         0.583408   0.750008   0.274193   x
            x  Si           9         0.750039   0.250004   0.200311   x
            x  Si          10         0.250029   0.749972   0.199897   x
            x  Si          11         0.418836   0.249771   0.124089   x
            x  Si          12         0.899985   0.749978   0.124760   x
            x  Si          13         0.583715   0.729943   0.045400   x
            x  Si          14         0.088604   0.248381   0.049990   x
            x  Si          15         0.914084   0.752693   0.518387   x
            x  Si          16         0.426451   0.240110   0.522102   x
            x  Si          17         0.100060   0.250024   0.442255   x
            x  Si          18         0.580866   0.749915   0.442866   x
            x  Si          19         0.749963   0.249996   0.367069   x
            x  Si          20         0.249940   0.750017   0.366650   x
            x  Si          21         0.416612   0.250014   0.292692   x
            x  Si          22         0.916663   0.750011   0.292697   x
            x  Si          23         0.083497   0.250026   0.218668   x
            x  Si          24         0.583294   0.750021   0.218687   x
            x  Si          25         0.744839   0.250070   0.145484   x
            x  Si          26         0.242804   0.750139   0.141605   x
            x  Si          27         0.426467   0.251079   0.068261   x
            x  Si          28         0.901773   0.747613   0.069051   x
            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 
------------------------------------------------------------------------ <-- SCF
SCF loop      Energy           Fermi           Energy gain       Timer   <-- SCF
                               energy          per atom          (sec)   <-- SCF
------------------------------------------------------------------------ <-- SCF
Initial  -2.98304086E+003  9.45338714E-001                      1823.31  <-- SCF
      1  -2.99913583E+003  9.79560691E-001   5.74820046E-001    1920.02  <-- SCF
      2  -2.99946713E+003  1.01446046E+000   1.18324096E-002    2030.02  <-- SCF
      3  -2.99948509E+003  1.02524541E+000   6.41308081E-004    2134.97  <-- SCF
      4  -2.99948740E+003  1.02907530E+000   8.26229187E-005    2240.54  <-- SCF
      5  -2.99948893E+003  1.03170213E+000   5.45154577E-005    2350.29  <-- SCF
      6  -2.99948901E+003  1.03175388E+000   2.94871450E-006    2416.19  <-- SCF
      7  -2.99948901E+003  1.03175388E+000   8.52130140E-008    2472.20  <-- SCF
------------------------------------------------------------------------ <-- SCF
 
Final energy, E             =  -2999.131245702     eV
Final free energy (E-TS)    =  -2999.489014446     eV
(energies not corrected for finite basis set)
 
NB est. 0K energy (E-0.5TS)      =  -2999.310130074     eV
 
 
 ******************************** Forces *********************************
 *                                                                       *
 *                      Cartesian components (eV/A)                      *
 * --------------------------------------------------------------------- *
 *              x                    y                    z              *
 *                                                                       *
 * Si   1     -0.80356              0.06560              1.27814         *
 * Si   2      0.01936             -0.26641              0.85116         *
 * Si   3      0.06695              0.00278             -0.38119         *
 * Si   4      0.18659             -0.00044              0.63665         *
 * Si   5      0.10700             -0.00150              0.01117         *
 * Si   6     -0.05267             -0.00208              0.08686         *
 * Si   7      0.03972             -0.00381             -0.06135         *
 * Si   8      0.04766             -0.00294             -0.01895         *
 * Si   9      0.05360             -0.00392              0.29987         *
 * Si  10      0.03077             -0.00902             -0.55495         *
 * Si  11     -0.45765             -0.03437              0.19802         *
 * Si  12      0.71394             -0.00452              0.19531         *
 * Si  13     -0.20774              0.69223              0.54657         *
 * Si  14     -0.63433             -0.08206              0.97571         *
 * Si  15      0.85346             -0.14734             -1.11827         *
 * Si  16      0.17753              0.27023             -0.64980         *
 * Si  17     -0.70790              0.01546              0.02647         *
 * Si  18      0.34927             -0.00603             -0.09317         *
 * Si  19     -0.01639              0.00954              0.68137         *
 * Si  20     -0.04318              0.01355             -0.20711         *
 * Si  21     -0.02371              0.01260             -0.05983         *
 * Si  22     -0.01899              0.01167             -0.04383         *
 * Si  23      0.09795              0.01418             -0.02979         *
 * Si  24     -0.07887              0.01502              0.04988         *
 * Si  25     -0.28027              0.02132             -0.51913         *
 * Si  26     -0.14250              0.02429              0.31049         *
 * Si  27      0.11787             -0.63761             -1.05937         *
 * Si  28      0.60608              0.03359             -1.35094         *
 *                                                                       *
 *************************************************************************
 
 ***************** Stress Tensor *****************
 *                                               *
 *          Cartesian components (GPa)           *
 * --------------------------------------------- *
 *             x             y             z     *
 *                                               *
 *  x     -2.105375     -0.048259     -0.077313  *
 *  y     -0.048259     -1.532628      0.045671  *
 *  z     -0.077313      0.045671     -1.371095  *
 *                                               *
 *  Pressure:    1.6697                          *
 *                                               *
 *************************************************
  
 +------------+-------------+-------------+-----------------+ <-- min BFGS
 |    Step    |   lambda    |   F.delta   |    enthalpy     | <-- min BFGS
 +------------+-------------+-------------+-----------------+ <-- min BFGS
 |  previous  |    0.000000 |    0.006113 |    -2998.949262 | <-- min BFGS
 | trial step |    1.000000 |    0.005645 |    -2999.113918 | <-- min BFGS
 |  line step |    3.659543 |    0.004401 |    -2999.489014 | <-- min BFGS
 +------------+-------------+-------------+-----------------+ <-- min BFGS
  
 BFGS: finished iteration     1 with enthalpy= -2.99948901E+003 eV
  
 +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
 | Parameter |      value      |    tolerance    |    units   | OK? | <-- BFGS
 +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
 |  dE/ion   |   1.927686E-002 |   2.000000E-005 |         eV | No  | <-- BFGS
 |  |F|max   |   1.511174E+000 |   1.000000E-002 |       eV/A | No  | <-- BFGS
 |  |dR|max  |   2.202829E-002 |   1.000000E-003 |          A | No  | <-- BFGS
 |   Smax    |   2.105375E+000 |   1.000000E-001 |        GPa | No  | <-- BFGS
 +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
 """
      self.lines = [x + '\n' for x in castep_output.split("\n")]
      self.al = AtomsList(self.lines, format='castep', abort=False)
      #self.al.loadall(progress=False)
     
   def testread(self):
      self.assertEqual(len(self.al), 3)
      self.assertEqual(self.al[0].n, 28)
      self.assertEqual(self.al[2].n, 28)
      self.assertAlmostEqual(self.al[0].energy, -2998.9492624099998)
      self.assertAlmostEqual(self.al[2].energy, -2999.4890144460001)
      self.assertEqual(self.al[0].force.shape, (3,28))
      self.assertEqual(self.al[2].force.shape, (3,28))


   def testabort(self):
      self.assertRaises(ValueError, AtomsList, self.lines, format='castep', abort=True)

   def testatoms_ref(self):
      a = Atoms.read(self.lines, format='castep', atoms_ref=self.al[0], abort=False)
      self.assertEqual(a.n, self.al[0].n)
      self.assertAlmostEqual(a.energy, self.al[0].energy)
      self.assertArrayAlmostEqual(a.force, self.al[0].force)


if __name__ == '__main__':
   unittest.main()


