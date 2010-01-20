import unittest
from numpy import all, unravel_index, loadtxt
from quippy import frange, farray
from StringIO import StringIO

def string_to_array(s):
   return loadtxt(StringIO(s)).T


class QuippyTestCase(unittest.TestCase):
   
   def assertArrayAlmostEqual(self, a, b, tol=1e-7):
      a = farray(a)
      b = farray(b)
      self.assertEqual(a.shape, b.shape)
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
         if hasattr(a, 'transpose_on_print') and a.transpose_on_print:
            print absdiff.T
         else:
            print absdiff
         print
         
         self.fail('Maximum abs difference between array elements is %e at location %r' % (absdiff.max(), loc))
   
