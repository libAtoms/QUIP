import unittest
from numpy import all, unravel_index
from quippy import frange

class QuippyTestCase(unittest.TestCase):
   
   def assertArrayAlmostEqual(self, a, b, tol=1e-7):
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
      
