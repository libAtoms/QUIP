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

from quippy import args_str, FortranArray, farray, fzeros, s2a, a2s
from numpy import dtype, int32
import unittest
from quippytest import *

from quippy import sort_array_i, sort_array_r

class TestLinearAlgebra(QuippyTestCase):

   def test_sort_array_i(self):
      a = farray([1,  4,  9,  7, -5,    3, -10, 11], dtype=int32)
      r = farray([1., 2., 3., 4., 5., 6.6, 7., 8.])
      sort_array_i(a, r_data=r)
      self.assertArrayAlmostEqual(a, [-10, -5,  1,   3,  4,  7,  9, 11])
      self.assertArrayAlmostEqual(r, [ 7., 5., 1., 6.6, 2., 4., 3., 8.])

   def test_sort_array_r(self):
      a = farray([1., 4., 9., 7., -5.,   3., -10., 11.])
      i = farray([1,  2,  3,  4,   5,    6,    7,   8 ], dtype=int32)
      r = farray([1., 2., 3., 4.,  5., 6.6,    7.,  8.])
      sort_array_r(a, i_data=i, r_data=r)
      self.assertArrayAlmostEqual(a, [-10., -5., 1.,  3., 4., 7., 9., 11.])
      self.assertArrayAlmostEqual(i, [  7.,  5., 1.,  6., 2., 4., 3.,  8.])
      self.assertArrayAlmostEqual(r, [  7.,  5., 1., 6.6, 2., 4., 3.,  8.])

      

if __name__ == '__main__':
   unittest.main()
