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

import unittest
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

        self.ref_int_1d = np.array([0, 0, 11, 24])
        self.ref_int_2d = np.array([[0, 0, 11],
                                    [1, 9, 0],
                                    [2, 0, 0],
                                    [3, 0, 24]])

        self.ref_str_1d = np.array(['a', 'b', 'c', 'd'])
        self.ref_str_2d = np.array([['a', 'b'],
                                    ['c', 'd'],
                                    ['a', 'r'],
                                    ['f', 'd']])

        self.at_base.arrays['dummy_real_1d'] = self.ref_real_1d
        self.at_base.arrays['dummy_real_2d'] = self.ref_real_2d
        self.at_base.arrays['dummy_int_1d'] = self.ref_int_1d
        self.at_base.arrays['dummy_int_2d'] = self.ref_int_2d

        self.ref_bool_1d = np.array([True, False, True, True])
        self.ref_bool_2d = np.zeros((2, 4), dtype=bool)

        self.ref_info = dict(real=0.2, logical_T=True, logical_F=False, integer=2)
        self.at_base.info.update(self.ref_info)

    def test_convert_add_from_info(self):
        at = self.at_base.copy()

        # fixme this is a hack to get the first element of the array
        # real
        raw_dict_output = self._convert_at_add_arrays_and_info(at, add_info='real')
        self.assertAlmostEqual(raw_dict_output['real'], self.ref_info['real'], delta=1E-06)

        # logical
        raw_dict_output = self._convert_at_add_arrays_and_info(at, add_info=['logical_T', 'logical_F'])
        self.assertEqual(raw_dict_output['logical_T'], self.ref_info['logical_T'])
        self.assertEqual(raw_dict_output['logical_F'], self.ref_info['logical_F'])

        # integer
        raw_dict_output = self._convert_at_add_arrays_and_info(at, add_info=['integer'])
        self.assertEqual(raw_dict_output['integer'], self.ref_info['integer'])

    def test_convert_add_arguent_types(self):
        # testing the different ways of specifying arguments

        # None
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_info=None)
        for key in self.at_base.arrays.keys():
            self.assertNotIn(key, raw_dict_output.keys())
        for key in self.at_base.info.keys():
            self.assertNotIn(key, raw_dict_output.keys())

        # string
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_info='real')
        self.assertIn('real', raw_dict_output.keys())

        # list
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_info=['real'])
        self.assertIn('real', raw_dict_output.keys())

        # tuple
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_info=('real'))
        self.assertIn('real', raw_dict_output.keys())

        # np.array
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_info=np.array(['real']))
        self.assertIn('real', raw_dict_output.keys())

        # True
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_info=True, add_arrays=True)
        self.assertIn('real', raw_dict_output.keys())
        for key in self.at_base.arrays.keys():
            if key in ['numbers', 'positions', 'momenta']:
                continue
            self.assertIn(key, raw_dict_output.keys())
        for key in self.at_base.info.keys():
            self.assertIn(key, raw_dict_output.keys())

    def test_convert_add_wrong_shape(self):
        """ Passing an array of the wrong shape"""
        # setup of a separate atoms for this test
        at = self.at_base.copy()
        at.arrays['wrong_shape'] = np.zeros((1, 6), dtype=float)
        self.assertRaises(RuntimeError, self._convert_at_add_arrays_and_info, at, add_arrays='wrong_shape')

    def test_convert_add_real(self):
        # default
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base)
        self.assertNotIn('dummy_real_2d', raw_dict_output.keys())
        self.assertNotIn('dummy_real_1d', raw_dict_output.keys())

        # 2d only
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_arrays='dummy_real_2d')
        self.assertArrayAlmostEqual(raw_dict_output['dummy_real_2d'], self.ref_real_2d.T, tol=1E-06)

        # 1d only
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_arrays=['dummy_real_1d'])
        self.assertArrayAlmostEqual(raw_dict_output['dummy_real_1d'], self.ref_real_1d, tol=1E-06)

    def test_convert_add_int(self):
        # 2d only
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_arrays='dummy_int_2d')
        self.assertArrayIntEqual(raw_dict_output['dummy_int_2d'], self.ref_int_2d.T)

        # 1d only
        raw_dict_output = self._convert_at_add_arrays_and_info(self.at_base, add_arrays=['dummy_int_1d'])
        self.assertArrayIntEqual(raw_dict_output['dummy_int_1d'], self.ref_int_1d)

    def test_convert_add_arrays_bool(self):
        # setup of a separate atoms for this test
        at = self.at_base.copy()
        at.arrays['dummy_bool_1d'] = self.ref_bool_1d
        at.arrays['dummy_bool_2d'] = self.ref_bool_2d

        # 1d bool array
        raw_dict_output = self._convert_at_add_arrays_and_info(at, add_arrays='dummy_bool_1d')
        self.assertArrayIntEqual(raw_dict_output['dummy_bool_1d'], self.ref_bool_1d)

        # 2d bool array
        self.assertRaises(TypeError, self._convert_at_add_arrays_and_info, at, add_arrays='dummy_bool_2d')

    def test_convert_add_arrays_unsupported_types(self):
        # setup of a separate atoms for this test
        at = self.at_base.copy()
        at.arrays.update(dummy_obj_2d=np.zeros((3, 4), dtype=object),
                         dummy_cplx_2d=np.zeros((3, 4), dtype=complex),
                         dummy_void_2d=np.zeros((3, 4), dtype=np.void))

        # all should raise an error, the
        self.assertRaises(TypeError, self._convert_at_add_arrays_and_info, at, add_arrays='dummy_obj_2d')
        self.assertRaises(TypeError, self._convert_at_add_arrays_and_info, at, add_arrays='dummy_void_2d')
        self.assertRaises(TypeError, self._convert_at_add_arrays_and_info, at, add_arrays='dummy_cplx_2d')

    @staticmethod
    def _convert_at_add_arrays_and_info(at, add_info=None, add_arrays=None):
        """
        Method for calling the converter and returning both the arrays and info in one dict.
        Yes, clashes are possible between keys but the tests don't need that right now.
        """
        at_quip = quippy.convert.ase_to_quip(at, add_arrays=add_arrays, add_info=add_info)
        raw_arrays = quippy.convert.get_dict_arrays(at_quip.properties)
        raw_info = quippy.convert.get_dict_arrays(at_quip.params)
        raw_arrays.update(raw_info)
        return raw_arrays


if __name__ == '__main__':
    unittest.main()
