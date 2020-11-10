# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2020
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

#from quippy.system_module import verbosity_push
#verbosity_push(1)

import unittest
import quippy
import numpy as np
import quippytest
from ase.build import bulk

from quippy.convert import ase_to_quip
from quippy.nye_tensor import nye_tensor

try:
    from matscipy.dislocation import (BCCScrew111Dislocation, 
                                      BCCEdge111Dislocation)
except ImportError:
    BCCScrew111Dislocation = None
    BCCEdge111Dislocation = None

@unittest.skipIf(BCCScrew111Dislocation is None, 'matscipy not available')
class TestNyeTensor(quippytest.QuippyTestCase):
    
    def build_disloc(self, cls):
        alat = 3.14339177996466
        C11 = 523.0266819809012
        C12 = 202.1786296941397
        C44 = 160.88179872237012
        
        d = cls(alat, C11, C12, C44)
        _, disloc = d.build_cylinder(radius=20.0)
        return d.unit_cell, disloc
    
    def test_screw(self):
        unit_cell, disloc = self.build_disloc(BCCScrew111Dislocation)
        alpha = nye_tensor(disloc, unit_cell, cutoff=3.0)
        asum = alpha.sum(axis=-1)        
        mask = np.ones((3, 3), dtype=bool)        
        mask[2, 2] = False # only alpha_33 should be significantly non-zero for screw
        assert np.all(np.abs(asum[mask]) < 1e-1)
        assert np.abs(asum[~mask]) > 1e-1
        
    
    def test_edge(self):
        unit_cell, disloc = self.build_disloc(BCCEdge111Dislocation)
        alpha = nye_tensor(disloc, unit_cell, cutoff=3.0)
        asum = alpha.sum(axis=-1)
        mask = np.ones((3, 3), dtype=bool)        
        mask[2, 0] = False # only alpha_31 should be significantly non-zero for edge
        assert np.all(np.abs(asum[mask]) < 1e-1)
        assert np.abs(asum[~mask]) > 1e-1


if __name__ == '__main__':
    unittest.main()
