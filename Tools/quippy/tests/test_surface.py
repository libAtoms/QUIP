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

from quippy import *
from quippy import surface
import unittest
from quippytest import *


class TestSurfaceRutile(QuippyTestCase):

    def setUp(self):
        self.rutile = Atoms("""6
Lattice="4.587000 0.000000 0.000000 0.000000 4.587000 0.000000 0.000000 0.000000 2.954000" Properties=species:S:1:pos:R:3:Z:I:1:primitive_index:I:1
Ti              0.00000000      0.00000000      0.00000000      22       1
Ti              2.29350000      2.29350000      1.47700000      22       2
O               1.39903500      1.39903500      0.00000000       8       3
O              -1.39903500     -1.39903500     -0.00000000       8       4
O               3.69253500      0.89446500      1.47700000       8       5
O               0.89446500      3.69253500      1.47700000       8       6""", format='string')

        self.rot = surface.crack_rotation_matrix(self.rutile, y=[1,-1,0],z=[1,1,0])
        self.slab = surface.orthorhombic_slab(self.rutile, rot=self.rot, periodicity=[0.0,0.0,0.0], max_nrep=20, verbose=False, graphics=False)

    def testrot(self):
        self.assertArrayAlmostEqual(self.rot.T, [[-0.,          0.,          1.        ],
                                                 [ 0.70710678, -0.70710678,  0.        ],
                                                 [ 0.70710678,  0.70710678,  0.        ]])

    def testlattice(self):
        self.assertArrayAlmostEqual(diag(self.slab.lattice), [2.954, 6.48699761, 6.48699761])


if __name__ == '__main__':
   unittest.main()
