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
from quippy.descriptors import Descriptor


class Test_NeighbourList(quippytest.QuippyTestCase):
    """
    Tests the NeighbourList with varying PBC settings as accessed by the wrapper.
    """

    def setUp(self):
        p0 = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    p0.append([i * 2 - 1, j * 2 - 1, k * 2 - 1])
        self.p0 = np.array(p0)

    def test_pbc(self):
        for cell in [[4] * 3,
                     [10] * 3,
                     [20] * 3,
                     [[4, 4, 0], [4, 0, 4], [0, 4, 4]],
                     [4, 10, 10]]:
            for pbc in [[True] * 3,
                        [False] * 3,
                        [True, False, False],
                        [False, True, False]]:
                init_size = None
                for offset_dir in [[1, 1, 1],
                                   [1, 0, 0]]:
                    for lat_offset in [-5.0, -1, -0.5, 0.0, 0.5, 1.0, 5.0]:
                        at = ase.Atoms(symbols=['Si'] * len(self.p0),
                                       positions=self.p0,
                                       cell=cell, pbc=pbc)
                        offset_v = np.dot(offset_dir, at.get_cell()) * lat_offset
                        at.positions += offset_v
                        desc = Descriptor('distance_2b cutoff=4.5')
                        sz = desc.sizes(at)
                        if init_size is None:
                            # print("cell", cell, "pbc", pbc, "offset dir", offset_dir,
                            #      "offset mag", lat_offset, "desc size", sz)
                            init_size = sz
                        self.assert_(init_size == sz)
                        #    print("cell", cell, "pbc", pbc, "offset dir", offset_dir,
                        #          "offset mag", lat_offset, "desc size", sz, "PROBLEM")

    def test_large_offset_no_pbc(self):
        cell = [10] * 3
        pbc = [False] * 3
        offset = 5e6
        at = ase.Atoms(symbols=['Si'] * len(self.p0),
                       positions=self.p0 + offset,
                       cell=cell,
                       pbc=pbc)
        desc = Descriptor('distance_2b cutoff=4.5')
        sz = desc.sizes(at)
        # print("cell", cell, "pbc", pbc, "desc size", sz)
        self.assert_(sz == (56, 112))


if __name__ == '__main__':
    unittest.main()
