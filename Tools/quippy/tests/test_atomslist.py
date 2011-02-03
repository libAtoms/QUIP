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
import unittest, itertools, sys, quippy
from quippytest import *

class TestAtomsList(QuippyTestCase):

   def setUp(self):
      self.listsrc = [diamond(5.44+0.01*x, 14) for x in range(5)]
      self.listal = AtomsList(self.listsrc)
      
      self.gensrc = (diamond(5.44+0.01*x,14) for x in range(5))
      self.genal = AtomsList(self.gensrc)
      
   def testgetitem(self):
      a1 = self.listal[0]
      a5 = self.listal[4]
      self.assertEqual(a1, diamond(5.44,14))
      self.assertEqual(a5, diamond(5.44+0.04,14))

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
      a1 = self.genal[0]
      a5 = self.genal[4]
      self.assertEqual(a1, diamond(5.44,14))
      self.assertEqual(a5, diamond(5.44+0.04,14))

   def testgengetslice(self):
      sl = self.genal[0:2]
      self.assert_(isinstance(sl, AtomsList))
      self.assertEqual(list(sl), list(self.genal)[0:2])

   def testgengetsliceopen(self):
      sl = self.genal[0:]
      self.assert_(isinstance(sl, AtomsList))
      self.assertEqual(list(sl), list(self.genal)[0:])

   def testgengetitemoutofrange(self):
      self.assertRaises(IndexError, self.genal.__getitem__, 6)

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
      class testwriter(object):
         def write(self, at):
            return at.n

      g = testwriter()
      self.assertEqual(self.listal.write(g), [8, 8, 8, 8, 8])

      g = testwriter()
      self.assertEqual(self.genal.write(g),  [8, 8, 8, 8, 8])

   def testwriteproperties(self):
      class testwriter_properties(object):
         def write(self, at, properties):
            return properties
      g = testwriter_properties()
      self.assertEqual(self.listal.write(g, properties=['pos']), [['pos']]*5)

   def testwriteproperties_exception(self):
      class testwriter(object):
         def write(self, at):
            return at.n
      g = testwriter()
      self.assertRaises(ValueError, self.listal.write, g, properties=['pos'])

   def testfrom_at(self):
      nl = AtomsList(self.listal[0])
      self.assertEqual(list(nl), [self.listal[0]])

   def testatomsreader_list(self):
      ar = AtomsReader(self.listal)
      self.assertEqual(list(ar), list(self.listal))

   def testatomsreader_gen(self):
      ar = AtomsReader(self.genal)
      self.assertEqual(list(ar), list(self.genal))

      
if __name__ == '__main__':
   unittest.main()
