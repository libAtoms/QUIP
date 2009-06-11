from quippy.ordereddict import *
import unittest

class TestOrderedDict(unittest.TestCase):

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

   def testpopitem(self):
      odict_copy = self.odict.copy()
      key, val = odict_copy.popitem()
      self.assertEqual(key, 3)
      self.assertEqual(val, 'c')
      self.assertEqual(odict_copy.keys(), [1,2])
      self.assertEqual(odict_copy.values(), ['a','b'])
      

   def testsetdefault(self):
      odict_copy = self.odict.copy()
      self.assertEqual(odict_copy.setdefault(-1, 'default'), 'default')
      self.assertEqual(odict_copy[-1],'default')

   def testrename(self):
      odict_copy = self.odict.copy()

      odict_copy.rename(1,4)
      self.assertEqual(odict_copy.keys(), [4,2,3])
      self.assertEqual(odict_copy.values(), ['a','b','c'])
      

if __name__ == '__main__':
   #unittest.main()
   suite = unittest.TestLoader().loadTestsFromTestCase(TestOrderedDict)
   unittest.TextTestRunner(verbosity=2).run(suite)
