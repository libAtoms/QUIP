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

from numpy import dtype
import unittest
import numpy as np
from quippytest import *

from quippy import *

class TestParamReader(QuippyTestCase):

    # here we test only the array parameters (param_register_multiple_integer
    # and param_register_multiple_real) since then the Python arrays are valid
    # targets which can be set from Fortran.

    def setUp(self):
        self.params = Dictionary()
        self.a = fzeros(3, np.int32)
        self.b = fzeros(3)
        param_register(self.params, 'a', '0 0 0', self.a, help_string="Param a")
        param_register(self.params, 'b', '0.0 0.0 0.0', self.b, help_string="Param b")

    def test_initial_values(self):
        self.assertEqual(list(self.a), [0,0,0])
        self.assertArrayAlmostEqual(self.b, [0., 0., 0.])

    def test_param_read_line(self):
        result = param_read_line(self.params, 'a={1 2 3} b={3.5 0 0}')
        self.assertEqual(result, 1)
        self.assertEqual(list(self.a), [1,2,3])
        self.assertArrayAlmostEqual(self.b, [3.5, 0., 0.])
    
