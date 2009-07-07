from quippy import *
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

   def testupdate_from_atoms(self):
      a = diamond(5.44, 14)
      c = castep.CastepCell()
      a.write(c)
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
                                             'false',
                                             'false',
                                             'none',
                                             'false',
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
      self.assertEqual((self.at.species.stripstrings() == 'Si').count(), 28)
      self.assertEqual(self.at.n, 28)

   def testatomsread(self):
      a = Atoms(self.geom_lines, format='geom')
      self.assertEqual(self.at.data, a.data)

   def testatomsreader(self):
      a, = AtomsList(self.geom_lines, format='geom')
      self.assertEqual(self.at.data, a.data)


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
      self.al.loadall(progress=False)
     
   def testread(self):
      self.assertEqual(len(self.al), 2)
      self.assertEqual(self.al[0].n, 28)
      self.assertEqual(self.al[1].n, 28)
      self.assertAlmostEqual(self.al[0].energy, -2998.58809358)
      self.assertAlmostEqual(self.al[1].energy, -2999.13124570)
      self.assertEqual(self.al[0].force.shape, (3,28))
      self.assertEqual(self.al[1].force.shape, (3,28))


   def testabort(self):
      al =  AtomsList(self.lines, format='castep')
      self.assertRaises(ValueError, al.loadall)

   def testfromcluster(self):
      a = Atoms.read(self.lines, format='castep', cluster=self.al[0], abort=False)
      self.assertEqual(a.n, self.al[0].n)
      self.assertAlmostEqual(a.energy, self.al[0].energy)
      self.assertArrayAlmostEqual(a.force, self.al[0].force)
      

def getTestSuite():
   tl = unittest.TestLoader()
   return unittest.TestSuite([tl.loadTestsFromTestCase(TestCastepCell),
                              tl.loadTestsFromTestCase(TestCastepParam),
                              tl.loadTestsFromTestCase(TestReadGeom),
                              tl.loadTestsFromTestCase(TestCastepOutput)])

if __name__ == '__main__':
   suite = getTestSuite()
   unittest.TextTestRunner(verbosity=2).run(suite)


