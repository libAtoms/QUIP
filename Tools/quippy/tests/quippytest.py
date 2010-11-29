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

import unittest, logging
from numpy import all, unravel_index, loadtxt, isnan
from quippy import frange, farray, FortranArray
from StringIO import StringIO

def string_to_array(s):
   return loadtxt(StringIO(s)).T

class QuippyTestCase(unittest.TestCase):

   def assertDictionariesEqual(self, d1, d2):
      if d1 != d2:
         if sorted(d1.keys()) != sorted(d2.keys()):
             self.fail('Dictionaries differ: d1.keys() (%r) != d2.keys() (%r)'  % (d1.keys(), d2.keys()))
         for key in d1:
            v1, v2 = d1[key], d2[key]
            if isinstance(v1, FortranArray):
               try:
                  self.assertArrayAlmostEqual(v1, v2)
               except AssertionError:
                  print key, v1, v2
                  raise
            else:
               if v1 != v2:
                  self.fail('Dictionaries differ: key=%s value1=%r value2=%r' % (key, v1, v2))
      
   def assertEqual(self, a, b):
      if a == b: return
      # Repeat comparison with debug-level logging
      import logging
      level = logging.root.level
      logging.root.setLevel(logging.DEBUG)
      a == b
      logging.root.setLevel(level)
      self.fail('%s != %s' % (a,b))
      
   
   def assertArrayAlmostEqual(self, a, b, tol=1e-7):
      a = farray(a)
      b = farray(b)
      self.assertEqual(a.shape, b.shape)

      if isnan(a).any() or isnan(b).any():
         print 'a'
         print a
         print 'b'
         print b
         self.fail('Not a number (NaN) found in array')
      
      if a.dtype.kind != 'f':
         self.assert_((a == b).all())
      else:
         absdiff = abs(a-b)
         if absdiff.max() > tol:
            loc = [x+1 for x in unravel_index(absdiff.argmax()-1, absdiff.shape) ]
            print 'a'
            print a
            print
            print 'b'
            print b
            print
            print 'Absolute difference'
            print absdiff
            self.fail('Maximum abs difference between array elements is %e at location %r' % (absdiff.max(), loc))
   

def skip(f):
   def g(self):
      logging.warning('skipping test %s' % f.__name__)
   return g
