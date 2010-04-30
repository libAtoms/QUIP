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

from quippy.ordereddict import *
import unittest
from quippytest import *

class TestOrderedDict(QuippyTestCase):

   def setUp(self):
      self.odict = OrderedDict({1:'a',
                                2:'b',
                                3:'c'})

   def testkeys(self):
      self.assertEqual(self.odict.keys(), [1,2,3])

   def testvalues(self):
      self.assertEqual(self.odict.values(), ['a','b','c'])

   def testiter(self):
      self.assertEqual([k for k in self.odict], [1,2,3])

   def testcopy(self):
      odict_copy = self.odict.copy()
      self.assertEqual(odict_copy.keys(), [1,2,3])
      self.assertEqual(odict_copy.values(), ['a','b','c'])

      from copy import copy
      odict_copy2 = copy(self.odict)
      self.assertEqual(odict_copy2.keys(), [1,2,3])
      self.assertEqual(odict_copy2.values(), ['a','b','c'])
      

   def testsetitem(self):
      odict_copy = self.odict.copy()
      odict_copy[4] = 'd'
      self.assertEqual(odict_copy[4], 'd')
      self.assertEqual(odict_copy.keys(), [1,2,3,4])
      
   def testclear(self):
      odict_copy = self.odict.copy()
      odict_copy.clear()
      self.assertEqual(odict_copy.keys(), [])
      self.assertEqual(odict_copy.values(), [])

   def testiterkeys(self):
      self.assertEqual(self.odict.keys(), [k for k in self.odict.iterkeys()])

   def testitervalues(self):
      self.assertEqual(self.odict.values(), [k for k in self.odict.itervalues()])

   def testiteritems(self):
      self.assertEqual(zip(self.odict.keys(),self.odict.values()),
                       [kv for kv in self.odict.iteritems()])

   def testitems(self):
      self.assertEqual(zip(self.odict.keys(),self.odict.values()),
                       self.odict.items())

   def testpopitem(self):
      odict_copy = self.odict.copy()
      key, val = odict_copy.popitem()
      self.assertEqual(key, 3)
      self.assertEqual(val, 'c')
      self.assertEqual(odict_copy.keys(), [1,2])
      self.assertEqual(odict_copy.values(), ['a','b'])
      key, val = odict_copy.popitem()
      key, val = odict_copy.popitem()
      self.assertRaises(KeyError, odict_copy.popitem)

   def testsetdefault(self):
      odict_copy = self.odict.copy()
      self.assertEqual(odict_copy.setdefault(-1, 'default'), 'default')
      self.assertEqual(odict_copy[-1],'default')

   def testrename(self):
      odict_copy = self.odict.copy()

      odict_copy.rename(1,4)
      self.assertEqual(odict_copy.keys(), [4,2,3])
      self.assertEqual(odict_copy.values(), ['a','b','c'])

      odict_copy.rename(4,4)
      self.assertEqual(odict_copy.keys(), [4,2,3])

      self.assertRaises(ValueError, odict_copy.rename, 4, 2)

   def testupdate(self):
      odict_copy = self.odict.copy()

      odict_copy.update({4:'d'})
      self.assertEqual(odict_copy.keys(), [1,2,3,4])
      self.assertEqual(odict_copy.values(), ['a','b','c','d'])
                                          

if __name__ == '__main__':
   unittest.main()
