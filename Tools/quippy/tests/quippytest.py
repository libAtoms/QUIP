import unittest
from numpy import all, unravel_index
from quippy import frange

class QuippyTestCase(unittest.TestCase):
   
   def assertArrayAlmostEqual(self, a, b, tol=1e-7):
      self.assertEqual(a.shape, b.shape)
      absdiff = abs(a-b)
      if absdiff.max() > tol:
         loc = [x+1 for x in unravel_index(absdiff.argmax()-1, absdiff.shape) ]
         self.fail('Maximum abs difference between array elements is %e at location %r' % (absdiff.max(), loc))
      
