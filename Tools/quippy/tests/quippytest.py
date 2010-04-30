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
   
