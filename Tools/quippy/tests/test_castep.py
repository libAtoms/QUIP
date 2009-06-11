from quippy import *
from StringIO import StringIO
import unittest
import xml.dom.minidom

class TestCastep(unittest.TestCase):
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

      

if __name__ == '__main__':
   #unittest.main()
   suite = unittest.TestLoader().loadTestsFromTestCase(TestCastep)
   unittest.TextTestRunner(verbosity=2).run(suite)
