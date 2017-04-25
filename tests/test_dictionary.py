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

from quippy.util import args_str
from quippy.farray import FortranArray, farray, fzeros, s2a, a2s
from numpy import dtype
import unittest
from quippytest import *

from quippy.dictionary import Dictionary

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
      self.assertEqual(' '.join(sorted(arg.split())), 'a=1 b c=F')

   def test_int_array(self):
      d = Dictionary()
      d['int'] = [i for i in range(10)]
      self.assertArrayAlmostEqual(FortranArray(d.get_array('int')), d['int'])

   def test_real_array(self):
      d = Dictionary()
      d['real'] = [float(i) for i in range(10)]
      self.assertArrayAlmostEqual(FortranArray(d.get_array('real')), d['real'])

   def test_complex_array(self):
      d = Dictionary()
      d['complex'] = fzeros(10,dtype='D')
      a = FortranArray(d.get_array('complex'))
      a[:] = 1+2j
      self.assertEqual(a.dtype.kind, 'c')
      self.assertArrayAlmostEqual(a, d['complex'])

   def test_logical_array(self):
      d = Dictionary()
      d['logical'] = fzeros(5,dtype='bool')
      a = FortranArray(d.get_array('logical'))
      a[:] = [True,False,False,True,True]
      self.assertEqual(a.dtype, dtype('int32'))    # Fortran logical represented as int32 internally
      self.assertArrayAlmostEqual(a, d['logical'])

   def test_string_array(self):
      d = Dictionary()
      d['string'] = ['one', 'two', 'three']
      a = FortranArray(d.get_array('string'))
      a[:,1] = farray(list('hello'),'S1')
      self.assertEqual(list(a2s(d['string'])), list(a2s(a)))

   def test_int_array2(self):
      d = Dictionary()
      d.add_array('int', 0, (3, 3))
      d_int = FortranArray(d.get_array('int'))
      d_int[1,1] = 3
      d_int[2,2] = 1
      d_int[:,3] = (4,1,5)
      self.assert_(all(d_int == d['int']))

   def test_real_array2(self):
      d = Dictionary()
      d.add_array('real', 0.0, (3, 3))
      d_real = FortranArray(d.get_array('real'))
      d_real[1,1] = 3.
      d_real[2,2] = 1.
      d_real[:,3] = (4.,1.,5.)
      self.assertArrayAlmostEqual(d_real, d['real'])

   def test_get_value(self):
      d = Dictionary()
      d['i'] = 1
      d['s'] = 'string'
      d['ia'] = [1,2,3]
      d['sa'] = ['string', 'longer string']
      d['la'] = [True, False]
      d['ia2'] = [[1,2,3],
                  [4,5,6]]
      for k in d:
         v1 = d[k]
         v2 = d.get_value(k)
         if hasattr(v1, '__iter__'):
            self.assert_(all(v1 == v2))
         else:
            self.assertEqual(v1, v2)

   def test_deepcopy(self):
      d1 = Dictionary()
      d1['ia2'] = [[1,2,3],
                   [4,5,6]]
      d2 = Dictionary()
      d2.deepcopy(d1) # copy from d1 into d2
      self.assert_(all(d1['ia2'] == d2['ia2']))
      d1['ia2'][1,1] = 0
      self.assert_(not all(d1['ia2'] == d2['ia2']))
      self.assert_(not all(d1._fpointer == d2._fpointer))

   def test_subset(self):
      d = Dictionary()
      d['i'] = 1
      d['s'] = 'string'
      d['ia'] = [1,2,3]
      d['sa'] = ['string', 'longer string']
      d['la'] = [True, False]
      d['ia2'] = [[1,2,3],
                  [4,5,6]]
      d2 = d.subset(['i', 'ia', 'la'])
      self.assert_(d2['i'] == 1 and list(d2['ia']) == [1,2,3] and list(d2['la']) == [True, False])

   def test_subset_bad_key(self):
      d = Dictionary()
      self.assertRaises(RuntimeError, d.subset, ['bad_key'])

   def test_bcast(self):
      from quippy.mpi_context import MPI_context, bcast
      d = Dictionary()
      mpi = MPI_context()
      bcast(mpi, d)

   def test_none_value(self):
      d = Dictionary('a=1 b=two c')
      self.assert_(d['a'] == 1)
      self.assert_(d['b'] == 'two')
      self.assert_(d['c'] is None)

   @skip
   def test_dict_in_dict(self):
      d = Dictionary()
      d2 = Dictionary()
      d2['a'] = 1
      d['d2'] = d2
      print d
      print d['d2']


if __name__ == '__main__':
   unittest.main()
