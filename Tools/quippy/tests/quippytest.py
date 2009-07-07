import unittest
from numpy import all

class QuippyTestCase(unittest.TestCase):
   
   def assertArrayAlmostEqual(self, a, b, tol=1e-8):
      self.assert_(abs(a - b).max() < tol)
      
