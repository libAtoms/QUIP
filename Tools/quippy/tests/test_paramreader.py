from quippy.paramreader import *
import unittest

class TestParamReader(unittest.TestCase):

   def testinitstring(self):
      params = ParamReader("a=1 b=2 c=3")
      self.assertEqual(params['a'], 1)
      self.assertEqual(params['b'], 2)
      self.assertEqual(params['c'], 3)
      self.assertEqual(params.keys(), ["a","b","c"])
      self.assertEqual(params.values(), [1,2,3])

   def testinitdict(self):
      params = ParamReader({'a':1, 'b':2})
      self.assertEqual(sorted(params.keys()), ['a','b'])

   def testinitparamreader(self):
      p1 = ParamReader("a=1 b=2 c=3")
      p2 = ParamReader(p1)

      self.assertEqual(p2['a'], 1)
      self.assertEqual(p2['b'], 2)
      self.assertEqual(p2['c'], 3)
      self.assertEqual(p2.keys(), ["a","b","c"])
      self.assertEqual(p2.values(), [1,2,3])

   def testinitlines(self):
      lines = ['a=1\n', 'b=2\n']
      p = ParamReader(lines)
      self.assertEqual(p.keys(), ['a','b'])
      self.assertEqual(p.values(), [1,2])

   def testiniterror(self):
      self.assertRaises(TypeError, ParamReader, 0)
                        
   def testquotedstring(self):
      params = ParamReader('a="one two" b="2.943"')
      self.assertEqual(params.keys(), ["a", "b"])
      self.assertEqual(params['a'], "one two")
      self.assertAlmostEqual(params['b'], 2.943)

   def testcopy(self):
      params = ParamReader('a="one two" b="2.943"')
      params2 = params.copy()
      from copy import copy
      params3 = copy(params)
      self.assertEqual(params.keys(), params2.keys())
      self.assertEqual(params.values(), params2.values())
      self.assertEqual(params.values(), params3.values())

   def testvaluelists(self):
      params = ParamReader('xyz={1.0 2.0 3.0} abc="1 2 3"')
      self.assertAlmostEqual(params['xyz'][0], 1.0)
      self.assertAlmostEqual(params['xyz'][1], 2.0)
      self.assertAlmostEqual(params['xyz'][2], 3.0)
      self.assertEqual(params['abc'], [1, 2, 3])

   def testrepr(self):
      params = ParamReader('a="one two" b="2.943"')
      self.assertEqual(repr(params), """ParamReader('a="one two" b=2.943')""")

   def testwrite(self):
      params = ParamReader('a="one two" b="2.943" c=44')
      from StringIO import StringIO
      s = StringIO()
      params.write(s)
      self.assertEqual(s.getvalue(), """a="one two"
b=2.943
c=44""")

   def testasstring(self):
      params = ParamReader('a="one two" b="2.943" c="44 45"')
      from StringIO import StringIO
      s = StringIO()
      s.write(params.asstring(' '))
      self.assertEqual(s.getvalue(), 'a="one two" b=2.943 c="44 45"')

   def testargs_str(self):
      arg = args_str({'a':1, 'b':True, 'c':False})
      self.assertEqual(' '.join(sorted(arg.split())), 'a=1 b=T c=F')
      

def getTestSuite():
   return unittest.TestLoader().loadTestsFromTestCase(TestParamReader)

if __name__ == '__main__':
   suite = getTestSuite()
   unittest.TextTestRunner(verbosity=2).run(suite)
