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
at = ase.build.bulk('C', 'diamond', 3.5)
at.append(ase.Atom('H', (0.2,0.2,0.1)))
desc = descriptors.Descriptor("soap cutoff=3 l_max=4 n_max=4 atom_sigma=0.5 n_Z=2 Z={1 6}")
at = quippy.Atoms(at)
at.set_cutoff(3.4)
at.calc_connect()
data = desc.calc(at, grad=True)
for key, val in data.items():
    print("{} {}".format(key, val.shape))
------------------------
cutoff (3,)
grad_index_0based (119, 2)
cutoff_grad (119, 3)
descriptor_index_0based (3, 1)
descriptor (3, 51)
grad (119, 3, 51)
    
"""


class Test_Descriptor(quippytest.QuippyTestCase):
    def setUp(self):
        at_h2 = ase.Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]])
        self.quip_at_h2 = quippy.convert.ase_to_quip(at_h2, None)

        # set up the atoms object for descriptor
        at_C2H = ase.build.bulk('C', 'diamond', 3.5)
        at_C2H.append(ase.Atom('H', (0.2, 0.2, 0.1)))
        self.at_C2H = at_C2H

        self.ref_shapes = {'cutoff': (3,),
                           'grad_index_0based': (119, 2),
                           'cutoff_grad': (119, 3),
                           'descriptor_index_0based': (3, 1),
                           'descriptor': (3, 51),
                           'grad': (119, 3, 51)}

        self.grad_array_ref = np.array([[[1.13442398e-01, -2.53763797e-04, 5.43671281e-06,
                                          -5.32119600e-04, -1.22629493e-04, -2.54495213e-02,
                                          -2.35215689e-03, 4.75523231e-05, -5.02800402e-03,
                                          -1.08561620e-03, -1.95363216e-01, -1.08721467e-02,
                                          1.99774343e-04, -2.37566525e-02, -4.73489941e-03,
                                          1.37877959e-01, 1.29099558e-03, -6.02464018e-05,
                                          1.60898786e-03, 1.20934568e-03, -6.20758739e-03,
                                          8.93220421e-03, -5.05352272e-04, 1.07208747e-02,
                                          8.71357869e-03, 8.31695444e-02, 5.34043119e-04,
                                          -7.42710259e-04, -2.67193323e-03, 3.30851049e-03,
                                          -1.26713676e-03, -1.17631427e-05, 5.17499647e-07,
                                          -1.41057631e-05, -1.10223791e-05, 5.22142751e-06,
                                          -8.09832380e-05, 4.52441394e-06, -9.39773456e-05,
                                          -7.95184362e-05, -1.07805976e-03, -2.24556230e-06,
                                          1.11276662e-05, 3.32533921e-05, -4.37921512e-05,
                                          6.98018370e-06, -1.28322312e-08, -8.24126629e-08,
                                          -2.07186884e-07, 2.89555893e-07, 0.00000000e+00],
                                         [1.13442398e-01, -2.53763797e-04, 5.43671281e-06,
                                          -5.32119600e-04, -1.22629493e-04, -2.54495213e-02,
                                          -2.35215689e-03, 4.75523231e-05, -5.02800402e-03,
                                          -1.08561620e-03, -1.95363216e-01, -1.08721467e-02,
                                          1.99774343e-04, -2.37566525e-02, -4.73489941e-03,
                                          1.37877959e-01, 1.29099558e-03, -6.02464018e-05,
                                          1.60898786e-03, 1.20934568e-03, -6.20758739e-03,
                                          8.93220421e-03, -5.05352272e-04, 1.07208747e-02,
                                          8.71357869e-03, 8.31695444e-02, 5.34043119e-04,
                                          -7.42710259e-04, -2.67193323e-03, 3.30851049e-03,
                                          -1.26713676e-03, -1.17631427e-05, 5.17499647e-07,
                                          -1.41057631e-05, -1.10223791e-05, 5.22142751e-06,
                                          -8.09832380e-05, 4.52441394e-06, -9.39773456e-05,
                                          -7.95184362e-05, -1.07805976e-03, -2.24556230e-06,
                                          1.11276662e-05, 3.32533921e-05, -4.37921512e-05,
                                          6.98018370e-06, -1.28322312e-08, -8.24126629e-08,
                                          -2.07186884e-07, 2.89555893e-07, 0.00000000e+00],
                                         [5.62089907e-02, -1.38457839e-04, 4.28026688e-05,
                                          -5.04141229e-04, 4.35011295e-06, -1.88932655e-02,
                                          -1.32007005e-03, 3.73714437e-04, -4.68615400e-03,
                                          7.28182168e-05, -1.07136020e-01, -6.25813743e-03,
                                          1.62136592e-03, -2.17410718e-02, 5.17941589e-04,
                                          7.80934164e-02, 2.83529951e-04, -4.81867697e-04,
                                          2.41304526e-03, 3.50601925e-04, 3.36108042e-03,
                                          2.47429899e-03, -3.13883417e-03, 1.64861857e-02,
                                          2.66528855e-03, 5.20708919e-02, 4.27468889e-03,
                                          1.38348816e-03, -6.99894308e-04, 2.08750362e-03,
                                          -7.13087093e-04, -2.60254817e-06, 4.33984399e-06,
                                          -2.18638435e-05, -3.08862541e-06, -5.27388207e-05,
                                          -2.24182986e-05, 2.84397065e-05, -1.49421137e-04,
                                          -2.36711867e-05, -6.69317942e-04, -5.21238036e-05,
                                          -1.56666324e-05, 8.45214919e-06, -2.82017984e-05,
                                          4.29940247e-06, 3.18081205e-07, 8.57747996e-08,
                                          -5.08881135e-08, 1.89987523e-07, 0.00000000e+00]],
                                        [[-2.14158073e-02, 3.79218050e-06, 7.23513681e-07,
                                          -7.84788979e-05, -3.50583379e-05, -2.89216230e-02,
                                          9.42399116e-05, 1.36864275e-05, -8.59430651e-04,
                                          -3.89273192e-04, -1.85999243e-02, 6.81736736e-04,
                                          8.43225045e-05, -4.63094307e-03, -2.14568411e-03,
                                          2.44466455e-02, 6.58481915e-04, 7.63834099e-05,
                                          -1.11484001e-03, -5.59330752e-04, 3.33924676e-02,
                                          3.63485980e-03, 2.89678116e-04, -7.41765752e-03,
                                          -4.14043902e-03, 4.03743618e-02, -7.19821888e-03,
                                          -1.90522488e-03, 1.93777581e-03, -2.42525548e-03,
                                          -1.95582360e-04, -5.54976344e-06, -6.34720740e-07,
                                          9.35790219e-06, 4.82053557e-06, -2.71507207e-04,
                                          -3.05150718e-05, -2.25446976e-06, 6.22363378e-05,
                                          3.57353505e-05, -4.90105993e-04, 8.71726430e-05,
                                          2.41405489e-05, -2.33132462e-05, 3.01491476e-05,
                                          2.97118096e-06, -5.27727787e-07, -1.52394728e-07,
                                          1.40208541e-07, -1.87271542e-07, 0.00000000e+00],
                                         [-7.13262538e-03, 1.21559754e-06, 2.15802281e-07,
                                          -2.57498259e-05, -1.21494351e-05, -9.63070849e-03,
                                          3.03991732e-05, 4.02761567e-06, -2.77321872e-04,
                                          -1.41451530e-04, -6.19188729e-03, 2.20220478e-04,
                                          2.46353751e-05, -1.47508784e-03, -8.10301233e-04,
                                          8.15552197e-03, 2.13260519e-04, 2.21561296e-05,
                                          -3.12256451e-04, -2.68959761e-04, 1.11385465e-02,
                                          1.17755131e-03, 8.23436392e-05, -2.08160223e-03,
                                          -2.00118852e-03, 1.34618090e-02, -2.32850330e-03,
                                          -5.66277047e-04, 5.10435977e-04, -1.24908308e-03,
                                          -6.52530791e-05, -1.79946963e-06, -1.85135523e-07,
                                          2.64035551e-06, 2.27137694e-06, -9.05709849e-05,
                                          -9.89870315e-06, -6.44420881e-07, 1.75971204e-05,
                                          1.69397481e-05, -1.63414980e-04, 2.82144027e-05,
                                          7.19229070e-06, -6.15369080e-06, 1.53724156e-05,
                                          9.90684444e-07, -1.70902096e-07, -4.55229079e-08,
                                          3.71006568e-08, -9.44582726e-08, 0.00000000e+00],
                                         [-7.13261775e-03, 1.25063479e-06, 2.30790120e-07,
                                          -2.57290744e-05, -1.21030381e-05, -9.63069593e-03,
                                          3.11208306e-05, 4.34616488e-06, -2.76819936e-04,
                                          -1.40330434e-04, -6.19187696e-03, 2.25093240e-04,
                                          2.67853396e-05, -1.47054829e-03, -8.02378204e-04,
                                          8.15553045e-03, 2.17634857e-04, 2.41402926e-05,
                                          -3.08810890e-04, -2.61269181e-04, 1.11385564e-02,
                                          1.19974567e-03, 9.21255151e-05, -2.04750454e-03,
                                          -1.96097972e-03, 1.34618137e-02, -2.39215608e-03,
                                          -5.97462401e-04, 5.95062268e-04, -1.35135072e-03,
                                          -6.52531544e-05, -1.83482509e-06, -2.01046527e-07,
                                          2.61227166e-06, 2.20985030e-06, -9.05710725e-05,
                                          -1.00752907e-05, -7.19750335e-07, 1.73150421e-05,
                                          1.66247639e-05, -1.63415039e-04, 2.89740666e-05,
                                          7.58164016e-06, -7.17685195e-06, 1.66062942e-05,
                                          9.90684816e-07, -1.75428638e-07, -4.79369204e-08,
                                          4.32667739e-08, -1.01876902e-07, 0.00000000e+00]]])

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
        desc = quippy.descriptors.Descriptor("soap cutoff=3 l_max=4 n_max=4 atom_sigma=0.5 n_Z=2 Z={1 6}")
        data = desc.calc(self.at_C2H, grad=True)

        # has_grad_data, pos, has_data - not tested
        self.assertTupleEqual(data['data'].shape, self.ref_shapes['descriptor'])
        self.assertTupleEqual(data['grad_data'].shape, self.ref_shapes['grad'])
        self.assertTupleEqual(data['covariance_cutoff'].shape, self.ref_shapes['cutoff'])
        self.assertTupleEqual(data['grad_covariance_cutoff'].shape, self.ref_shapes['cutoff_grad'])
        self.assertTupleEqual(data['ii'].shape, (self.ref_shapes['grad_index_0based'][0],))

        # test the gradient's values
        self.assertArrayAlmostEqual(data['grad_data'][:2], self.grad_array_ref)


if __name__ == '__main__':
    unittest.main()
