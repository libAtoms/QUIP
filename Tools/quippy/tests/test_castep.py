from quippy import *
from StringIO import StringIO
import unittest
import xml.dom.minidom

class TestCastepCell(unittest.TestCase):
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
         H 0.0 0.0 0.0
         %ENDBLOCK POSITIONS_ABS""".split("\n")

      self.cell = castep.CastepCell(self.cell_lines)


   def testkeys(self):
      self.assertEqual(self.cell.keys(),['fix_all_cell','SPECIES_POT', 'LATTICE_CART', 'POSITIONS_ABS'])

   def testrepr(self):
      self.assertEqual(repr(self.cell), "CastepCell({'fix_all_cell': 'TRUE', 'POSITIONS_ABS': ['H 0.0 0.0 0.0'], 'SPECIES_POT': ['H  H_00PBE.usp', 'C  C_00PBE.usp'], 'LATTICE_CART': ['10.000000 0.000000 0.000000', '0.000000 10.000000 0.000000', '0.000000 0.000000 10.000000']})")

   def testvalues(self):
      self.assertEqual(self.cell.values(),['TRUE',
                                           ['H  H_00PBE.usp', 'C  C_00PBE.usp'],
                                           ['10.000000 0.000000 0.000000',
                                            '0.000000 10.000000 0.000000',
                                            '0.000000 0.000000 10.000000'],
                                           ['H 0.0 0.0 0.0']])

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

   def test_to_atoms_without_lattice(self):
      del self.cell['LATTICE_CART']
      self.assertRaises(ValueError, self.cell.to_atoms)

   def test_to_atoms_without_pos(self):
      del self.cell['POSITIONS_ABS']
      self.assertRaises(ValueError, self.cell.to_atoms)

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


class TestCastepParam(unittest.TestCase):

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

      self.castep_lines = [s+"\n" for s in castep_output.split("\n")]


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

   
class TestReadGeom(unittest.TestCase):

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

      self.at = castep.read_geom(self.geom_lines)

   def testspecies(self):
      self.assert_(all(self.at.species.stripstrings() == 'Si'))
      

def getTestSuite():
   tl = unittest.TestLoader()
   return unittest.TestSuite([tl.loadTestsFromTestCase(TestCastepCell),
                              tl.loadTestsFromTestCase(TestCastepParam),
                              tl.loadTestsFromTestCase(TestReadGeom)])

if __name__ == '__main__':
   suite = getTestSuite()
   unittest.TextTestRunner(verbosity=2).run(suite)


castep_output=""" +-------------------------------------------------+
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
