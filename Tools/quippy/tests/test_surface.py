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
from numpy import *
from quippy.structure_tools import * 
import quippy.sio2
import unittest
from quippytest import *

class TestSurfaceSilicon(QuippyTestCase):

    def setUp(self):
        self.dia = diamond(5.44, 14)

    def test_111_11b0(self):
        y, z, shift = (1,1,1), (1,-1,0), None
        surf = orthorhombic_slab(self.dia, rot=rotation_matrix(self.dia, y, z), shift=shift, verbose=False)
        self.assertArrayAlmostEqual(surf.lattice, diag([6.6626121, 9.42235639, 3.84666089]))
       
    def test_110_001b(self):
        y, z, shift = (1,1,0), (0,0,-1), None
        surf = orthorhombic_slab(self.dia, rot=rotation_matrix(self.dia, y, z), shift=shift, verbose=False)
        self.assertArrayAlmostEqual(surf.lattice, diag([3.84666089, 3.84666089, 5.44]))

class TestSurfaceQuartz(QuippyTestCase):

    def setUp(self):
        self.aq = alpha_quartz(**quippy.sio2.quartz_params['ASAP_JRK'])
        self.h = 10.0
        self.d = 15.0

    def test_0001_1b10(self):
        y, z, shift = (0,0,0,1), (-1,1,0), None
        surf = orthorhombic_slab(self.aq, rot=rotation_matrix(self.aq, y, z), shift=shift, verbose=False)
        self.assertArrayAlmostEqual(surf.lattice, diag([4.84038097, 5.328524, 8.38378577]))

    def test_101b0_010(self):
        y, z, shift = (1,0,-1,0), (0,1,0), [0.0, 0.07, 0.0]
        surf = orthorhombic_slab(self.aq, rot=rotation_matrix(self.aq, y, z), shift=shift, verbose=False)
        self.assertArrayAlmostEqual(surf.lattice, diag([5.328524, 8.38378577, 4.84038097]))

    def test_101b1_010(self):
        y, z, shift = [1,0,-1,1], [0,1,0], None
        surf = orthorhombic_slab(self.aq, rot=rotation_matrix(self.aq, y, z), shift=shift,
                                         periodicity=[0.0, self.h, 0.0], vacuum=[0.0, self.d, 0.0], verbose=False)
        self.assertArrayAlmostEqual(surf.lattice, diag([13.55951828 , self.h+self.d, 4.84038097]))
        

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


    def test_11b0_110(self):
        slab = orthorhombic_slab(self.rutile, rot=rotation_matrix(self.rutile, y=[1,-1,0],z=[1,1,0]),
                                              max_nrep=20, verbose=False, graphics=False)
        self.assertArrayAlmostEqual(slab.lattice, diag([2.954, 6.48699761, 6.48699761]))


if __name__ == '__main__':
   unittest.main()
