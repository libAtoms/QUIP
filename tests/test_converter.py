# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright T. K. Stenczel 2019
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

import quippy
import quippytest
import ase
import numpy as np


class Test_Converter(quippytest.QuippyTestCase):
    def setUp(self):
        # this is a base atoms object, which is used afterwards
        self.at_base = ase.Atoms('HCOO', positions=[[0, 0, 0], [1, 0, 0], [2., 0, 0], [3., 0, 0]], cell=[3, 3, 3],
                                 pbc=True)
        self.at_base.set_momenta([[0, 0, 0], [1, 0, 0], [2., 0, 0], [5., 0, 0]])

        self.ref_real_1d = np.array([0., 0., 11., 24.])
        self.ref_real_2d = np.array([[0., 0., 11.],
                                     [0.1, 6.9, 0.],
                                     [0.2, 0., 0.],
                                     [0.03, 0., 24.]])

        self.at_base.arrays['dummy_real_1d'] = self.ref_real_1d
        self.at_base.arrays['dummy_real_2d'] = self.ref_real_2d

    def test_convert_add_properties(self):
        # default
        at_quip = quippy.convert.ase_to_quip(self.at_base)
        raw_arrays = quippy.convert.get_dict_arrays(at_quip.properties)
        self.assertNotIn('dummy_real_2d', raw_arrays.keys())
        self.assertNotIn('dummy_real_1d', raw_arrays.keys())

        # 2d only
        at_quip = quippy.convert.ase_to_quip(self.at_base, add_arrays='dummy_real_2d')
        raw_arrays = quippy.convert.get_dict_arrays(at_quip.properties)
        self.assertArrayAlmostEqual(raw_arrays['dummy_real_2d'], self.ref_real_2d.T, tol=1E-06)

        # 1d only
        at_quip = quippy.convert.ase_to_quip(self.at_base, add_arrays=['dummy_real_1d'])
        raw_arrays = quippy.convert.get_dict_arrays(at_quip.properties)
        self.assertArrayAlmostEqual(raw_arrays['dummy_real_1d'], self.ref_real_1d, tol=1E-06)


if __name__ == '__main__':
    unittest.main()
