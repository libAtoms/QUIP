from quippy import args_str, farray
import unittest
from quippytest import *

from quippy import Dictionary

class TestDictionary(QuippyTestCase):

   def testinitstring(self):
      params = Dictionary("a=1 b=2 c=3")
      self.assertEqual(params['a'], 1)
      self.assertEqual(params['b'], 2)
      self.assertEqual(params['c'], 3)
      self.assertEqual(params.keys(), ["a","b","c"])
      self.assertEqual(params.values(), [1,2,3])

   def testinitdict(self):
      params = Dictionary({'a':1, 'b':2})
      self.assertEqual(sorted(params.keys()), ['a','b'])

   def testinitdictionary(self):
      p1 = Dictionary("a=1 b=2 c=3")
      p2 = Dictionary(p1)

      self.assertEqual(p2['a'], 1)
      self.assertEqual(p2['b'], 2)
      self.assertEqual(p2['c'], 3)
      self.assertEqual(p2.keys(), ["a","b","c"])
      self.assertEqual(p2.values(), [1,2,3])

   def testinitlines(self):
      lines = ['a=1\n', 'b=2\n']
      p = Dictionary(lines)
      self.assertEqual(p.keys(), ['a','b'])
      self.assertEqual(p.values(), [1,2])

   def testiniterror(self):
      self.assertRaises(TypeError, Dictionary, 0)
                        
   def testquotedstring(self):
      params = Dictionary('a="one two" b="2.943"')
      self.assertEqual(params.keys(), ["a", "b"])
      self.assertEqual(params['a'], "one two")
      self.assertAlmostEqual(params['b'], 2.943)

   def testcopy(self):
      params = Dictionary('a="one two" b="2.943"')
      params2 = params.copy()
      from copy import copy
      params3 = copy(params)
      self.assertEqual(params.keys(), params2.keys())
      self.assertEqual(params.values(), params2.values())
      self.assertEqual(params.values(), params3.values())

   def testvaluelists(self):
      params = Dictionary('xyz={1.0 2.0 3.0} abc="1 2 3"')
      self.assertAlmostEqual(params['xyz'][1], 1.0)
      self.assertAlmostEqual(params['xyz'][2], 2.0)
      self.assertAlmostEqual(params['xyz'][3], 3.0)
      self.assertArrayAlmostEqual(params['abc'], farray([1, 2, 3]))

   def testrepr(self):
      params = Dictionary('a="one two" b="2.943"')
      self.assertEqual(repr(params), """Dictionary('a="one two" b=2.943')""")

   def testwrite(self):
      params = Dictionary('a="one two" b="2.943" c=44')
      from StringIO import StringIO
      s = StringIO()
      params.write(s)
      self.assertEqual(s.getvalue(), """a="one two"
b=2.943
c=44""")

   def testasstring(self):
      params = Dictionary('a="one two" b="2.943" c="44 45"')
      from StringIO import StringIO
      s = StringIO()
      s.write(params.asstring(' '))
      self.assertEqual(s.getvalue(), 'a="one two" b=2.943 c="44 45"')

   def testargs_str(self):
      arg = args_str({'a':1, 'b':True, 'c':False})
      self.assertEqual(' '.join(sorted(arg.split())), 'a=1 b=T c=F')
      

if __name__ == '__main__':
   unittest.main()
