# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
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
import os

import quippy
import quippytest
import ase

import numpy as np
import ase.build

# ref data made in python2, with the old version of quippy, to reproduce the same behaviour
"""
import quippy
from quippy import descriptors
import ase
import ase.build
ase.Atoms(positions=[[0., 0., 0.], [0.875, 0.875, 0.875], [0.2, 0.2, 0.1]],
          symbols=['C', 'C', 'H'],
          cell=[[0.0, 1.75, 1.75], [1.75, 0.0, 1.75], [1.75, 1.75, 0.0]],
          pbc=True)
desc = quippy.descriptors.Descriptor("soap cutoff=1.3 l_max=4 n_max=4 atom_sigma=0.5 n_Z=2 Z={1 6}")
at = quippy.Atoms(at)
at.set_cutoff(3.)
at.calc_connect()
data = desc.calc(at, grad=True)
for key, val in data.iteritems():
    print("{} {}".format(key, val.shape))
------------------------
cutoff (3,)
grad_index_0based (119, 2)
cutoff_grad (119, 3)
descriptor_index_0based (3, 1)
descriptor (3, 51)
grad (119, 3, 51)
    
"""


@unittest.skipIf(os.environ['HAVE_GAP'] != '1', 'GAP support not enabled')
class Test_Descriptor(quippytest.QuippyTestCase):
    def setUp(self):
        at_h2 = ase.Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]])
        self.quip_at_h2 = quippy.convert.ase_to_quip(at_h2, None)

        # set up the atoms object for descriptor
        self.at_C2H = ase.Atoms(positions=[[0., 0., 0.], [0.875, 0.875, 0.875], [0.2, 0.2, 0.1]],
                                symbols=['C', 'C', 'H'],
                                cell=[[0.0, 1.75, 1.75], [1.75, 0.0, 1.75], [1.75, 1.75, 0.0]],
                                pbc=True)

        self.ref_shapes = {'cutoff': (3,),
                           'grad_index_0based': (7, 2),
                           'cutoff_grad': (7, 3),
                           'descriptor_index_0based': (3, 1),
                           'descriptor': (3, 51),
                           'grad': (7, 3, 51)}

        self.ref_grad_index_0based = np.array([[0, 0],
                                               [0, 2],
                                               [1, 1],
                                               [1, 2],
                                               [2, 2],
                                               [2, 1],
                                               [2, 0]], dtype="int32")

        self.ref_grad_array = np.array([[[1.88046072e-02, -7.69208358e-03, -1.26665811e-03,
                                          -8.34260643e-05, -2.94507806e-06, -2.78754894e-01,
                                          -2.60118673e-02, -4.28485490e-03, -2.82433778e-04,
                                          -9.98579772e-06, -1.79631756e-02, -4.39814001e-02,
                                          -7.24740987e-03, -4.78081025e-04, -1.69292884e-05,
                                          3.51431392e-02, 3.30989891e-03, 5.40034605e-04,
                                          3.48202791e-05, 1.17674588e-06, 3.21367088e-03,
                                          7.91457886e-03, 1.29176534e-03, 8.33551581e-05,
                                          2.82133185e-06, -2.87464795e-04, -7.12099923e-04,
                                          -1.15118476e-04, -7.26646328e-06, -2.35084115e-07,
                                          -3.44483503e-03, -3.24458191e-04, -5.29357344e-05,
                                          -3.41287931e-06, -1.15316199e-07, -3.15017852e-04,
                                          -7.75839390e-04, -1.26622528e-04, -8.16998322e-06,
                                          -2.76478780e-07, 3.98504667e-05, 9.87188049e-05,
                                          1.59583266e-05, 1.00722553e-06, 3.25796076e-08,
                                          -2.76218118e-06, -6.84272117e-06, -1.10611343e-06,
                                          -6.98072244e-08, -2.25755541e-09, 0.00000000e+00],
                                         [1.88046072e-02, -7.69208358e-03, -1.26665811e-03,
                                          -8.34260643e-05, -2.94507806e-06, -2.78754894e-01,
                                          -2.60118673e-02, -4.28485490e-03, -2.82433778e-04,
                                          -9.98579772e-06, -1.79631756e-02, -4.39814001e-02,
                                          -7.24740987e-03, -4.78081025e-04, -1.69292884e-05,
                                          3.51431392e-02, 3.30989891e-03, 5.40034605e-04,
                                          3.48202791e-05, 1.17674588e-06, 3.21367088e-03,
                                          7.91457886e-03, 1.29176534e-03, 8.33551581e-05,
                                          2.82133185e-06, -2.87464795e-04, -7.12099923e-04,
                                          -1.15118476e-04, -7.26646328e-06, -2.35084115e-07,
                                          -3.44483503e-03, -3.24458191e-04, -5.29357344e-05,
                                          -3.41287931e-06, -1.15316199e-07, -3.15017852e-04,
                                          -7.75839390e-04, -1.26622528e-04, -8.16998322e-06,
                                          -2.76478780e-07, 3.98504667e-05, 9.87188049e-05,
                                          1.59583266e-05, 1.00722553e-06, 3.25796076e-08,
                                          -2.76218118e-06, -6.84272117e-06, -1.10611343e-06,
                                          -6.98072244e-08, -2.25755541e-09, 0.00000000e+00],
                                         [9.40230362e-03, -3.84604179e-03, -6.33329055e-04,
                                          -4.17130322e-05, -1.47253903e-06, -1.39377447e-01,
                                          -1.30059337e-02, -2.14242745e-03, -1.41216889e-04,
                                          -4.99289886e-06, -8.98158780e-03, -2.19907000e-02,
                                          -3.62370494e-03, -2.39040513e-04, -8.46464422e-06,
                                          1.75715696e-02, 1.65494945e-03, 2.70017303e-04,
                                          1.74101395e-05, 5.88372939e-07, 1.60683544e-03,
                                          3.95728943e-03, 6.45882669e-04, 4.16775790e-05,
                                          1.41066592e-06, -1.43732397e-04, -3.56049961e-04,
                                          -5.75592382e-05, -3.63323164e-06, -1.17542058e-07,
                                          -1.72241751e-03, -1.62229095e-04, -2.64678672e-05,
                                          -1.70643966e-06, -5.76580993e-08, -1.57508926e-04,
                                          -3.87919695e-04, -6.33112642e-05, -4.08499161e-06,
                                          -1.38239390e-07, 1.99252333e-05, 4.93594025e-05,
                                          7.97916328e-06, 5.03612764e-07, 1.62898038e-08,
                                          -1.38109059e-06, -3.42136059e-06, -5.53056716e-07,
                                          -3.49036122e-08, -1.12877770e-09, 0.00000000e+00]],

                                        [[-1.88046072e-02, 7.69208358e-03, 1.26665811e-03,
                                          8.34260643e-05, 2.94507806e-06, 2.78754894e-01,
                                          2.60118673e-02, 4.28485490e-03, 2.82433778e-04,
                                          9.98579772e-06, 1.79631756e-02, 4.39814001e-02,
                                          7.24740987e-03, 4.78081025e-04, 1.69292884e-05,
                                          -3.51431392e-02, -3.30989891e-03, -5.40034605e-04,
                                          -3.48202791e-05, -1.17674588e-06, -3.21367088e-03,
                                          -7.91457886e-03, -1.29176534e-03, -8.33551581e-05,
                                          -2.82133185e-06, 2.87464795e-04, 7.12099923e-04,
                                          1.15118476e-04, 7.26646328e-06, 2.35084115e-07,
                                          3.44483503e-03, 3.24458191e-04, 5.29357344e-05,
                                          3.41287931e-06, 1.15316199e-07, 3.15017852e-04,
                                          7.75839390e-04, 1.26622528e-04, 8.16998322e-06,
                                          2.76478780e-07, -3.98504667e-05, -9.87188049e-05,
                                          -1.59583266e-05, -1.00722553e-06, -3.25796076e-08,
                                          2.76218118e-06, 6.84272117e-06, 1.10611343e-06,
                                          6.98072244e-08, 2.25755541e-09, 0.00000000e+00],
                                         [-1.88046072e-02, 7.69208358e-03, 1.26665811e-03,
                                          8.34260643e-05, 2.94507806e-06, 2.78754894e-01,
                                          2.60118673e-02, 4.28485490e-03, 2.82433778e-04,
                                          9.98579772e-06, 1.79631756e-02, 4.39814001e-02,
                                          7.24740987e-03, 4.78081025e-04, 1.69292884e-05,
                                          -3.51431392e-02, -3.30989891e-03, -5.40034605e-04,
                                          -3.48202791e-05, -1.17674588e-06, -3.21367088e-03,
                                          -7.91457886e-03, -1.29176534e-03, -8.33551581e-05,
                                          -2.82133185e-06, 2.87464795e-04, 7.12099923e-04,
                                          1.15118476e-04, 7.26646328e-06, 2.35084115e-07,
                                          3.44483503e-03, 3.24458191e-04, 5.29357344e-05,
                                          3.41287931e-06, 1.15316199e-07, 3.15017852e-04,
                                          7.75839390e-04, 1.26622528e-04, 8.16998322e-06,
                                          2.76478780e-07, -3.98504667e-05, -9.87188049e-05,
                                          -1.59583266e-05, -1.00722553e-06, -3.25796076e-08,
                                          2.76218118e-06, 6.84272117e-06, 1.10611343e-06,
                                          6.98072244e-08, 2.25755541e-09, 0.00000000e+00],
                                         [-9.40230362e-03, 3.84604179e-03, 6.33329055e-04,
                                          4.17130322e-05, 1.47253903e-06, 1.39377447e-01,
                                          1.30059337e-02, 2.14242745e-03, 1.41216889e-04,
                                          4.99289886e-06, 8.98158780e-03, 2.19907000e-02,
                                          3.62370494e-03, 2.39040513e-04, 8.46464422e-06,
                                          -1.75715696e-02, -1.65494945e-03, -2.70017303e-04,
                                          -1.74101395e-05, -5.88372939e-07, -1.60683544e-03,
                                          -3.95728943e-03, -6.45882669e-04, -4.16775790e-05,
                                          -1.41066592e-06, 1.43732397e-04, 3.56049961e-04,
                                          5.75592382e-05, 3.63323164e-06, 1.17542058e-07,
                                          1.72241751e-03, 1.62229095e-04, 2.64678672e-05,
                                          1.70643966e-06, 5.76580993e-08, 1.57508926e-04,
                                          3.87919695e-04, 6.33112642e-05, 4.08499161e-06,
                                          1.38239390e-07, -1.99252333e-05, -4.93594025e-05,
                                          -7.97916328e-06, -5.03612764e-07, -1.62898038e-08,
                                          1.38109059e-06, 3.42136059e-06, 5.53056716e-07,
                                          3.49036122e-08, 1.12877770e-09, 0.00000000e+00]]])

    def test_set_cutoff(self):
        # todo: move this to the test of the atoms objects
        self.quip_at_h2.set_cutoff(5)
        self.assertAlmostEqual(self.quip_at_h2.cutoff, 5.0, msg='cutoff not set to the given value', delta=0.01)

    def test_calc_connect(self):
        # todo: expand this to check for the result of the connectivity, not just the fact that the method exists
        with self.assertRaises(RuntimeError):
            self.quip_at_h2.calc_connect()

    def test_descriptor_type(self):
        generic_descriptor = quippy.descriptors_module.descriptor('distance_2b cutoff=4.0')
        self.assertIsInstance(generic_descriptor, quippy.descriptors_module.descriptor, 'Wrong type defined by generic '
                                                                                        'descriptor initialiser')
        # gives different if initialised directly, so let's check that too
        distance_2b = quippy.descriptors_module.distance_2b("distance_2b cutoff=4.0")
        self.assertIsInstance(distance_2b, quippy.descriptors_module.distance_2b, 'Wrong type defined by direct '
                                                                                  'descriptor initialiser')

    def test_descriptor_output(self):
        desc = quippy.descriptors.Descriptor("soap cutoff=1.3 l_max=4 n_max=4 atom_sigma=0.5 n_Z=2 Z={1 6}")
        data = desc.calc(self.at_C2H, grad=True)

        # has_grad_data, pos, has_data - not tested
        self.assertTupleEqual(data['data'].shape, self.ref_shapes['descriptor'])
        self.assertTupleEqual(data['grad_data'].shape, self.ref_shapes['grad'])
        self.assertTupleEqual(data['covariance_cutoff'].shape, self.ref_shapes['cutoff'])
        self.assertTupleEqual(data['grad_covariance_cutoff'].shape, self.ref_shapes['cutoff_grad'])
        self.assertTupleEqual(data['grad_index_0based'].shape, self.ref_shapes['grad_index_0based'])

        # test the indices
        self.assertArrayIntEqual(data["grad_index_0based"], self.ref_grad_index_0based)

        # test the gradient's values
        self.assertArrayAlmostEqual(data['grad_data'][:2], self.ref_grad_array)


if __name__ == '__main__':
    unittest.main()
