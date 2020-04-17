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

import unittest

import ase
import numpy as np
import quippytest
from quippy.potential import Potential

xml_string = """
<SW_params n_types="1">
      <per_type_data type="1" atomic_num="14" />
      <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
      p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />
      <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.1675" />
</SW_params>
"""

class Test_NeighbourList(quippytest.QuippyTestCase):
    """
    Tests the NeighbourList with varying PBC settings as accessed by the wrapper.
    """

    def setUp(self):
        p0 = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    p0.append([i * 2 - 1,
                               j * 2 - 1,
                               k * 2 - 1])
        self.p0 = np.array(p0)
        self.sw_pot = Potential('IP SW', param_str=xml_string)

    def test_pbc(self):
        for cell in [tuple([4] * 3),
                     tuple([10] * 3),
                     tuple([20] * 3),
                     ((4, 4, 0), (4, 0, 4), (0, 4, 4)),
                     (4, 10, 10)]:
            for pbc in [tuple([True] * 3),
                        tuple([False] * 3),
                        (True, False, False),
                        (False, True, False)]:
                init_energy = None
                for offset_dir in [(1, 1, 1),
                                   (1, 0, 0)]:
                    for lat_offset in [-5.0, -1, -0.5, 0.0, 0.5, 1.0, 5.0]:
                        at = ase.Atoms(symbols=['Si'] * len(self.p0),
                                       positions=self.p0,
                                       cell=cell, pbc=pbc)
                        at.calc = self.sw_pot
                        offset_v = np.dot(offset_dir, at.get_cell()) * lat_offset
                        at.positions += offset_v
                        if init_energy is None:
                            init_energy = at.get_potential_energy()
                        self.assertAlmostEqual(init_energy,
                                               at.get_potential_energy())

    def test_large_offset_no_pbc(self):
        cell = [10] * 3
        pbc = [False] * 3
        offset = 5e6
        at = ase.Atoms(symbols=['Si'] * len(self.p0),
                       positions=self.p0 + offset,
                       cell=cell,
                       pbc=pbc)
        at.calc = self.sw_pot
        e = at.get_potential_energy()
        self.assertAlmostEqual(e, 20.863001205176083)


if __name__ == '__main__':
    unittest.main()
