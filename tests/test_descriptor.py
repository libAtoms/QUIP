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


class Test_Descriptor(quippytest.QuippyTestCase):
    def setUp(self):
        at = ase.Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]])
        self.quip_at = quippy.convert.ase_to_quip(at, None)

    def test_set_cutoff(self):
        # todo: move this to the test of the atoms objects
        self.quip_at.set_cutoff(5)
        self.assertAlmostEqual(self.quip_at.cutoff, 5.0, msg='cutoff not set to the given value', delta=0.01)

    def test_calc_connect(self):
        # todo: expand this to check for the result of the connectivity, not just the fact that the method exists
        with self.assertRaises(RuntimeError):
            self.quip_at.calc_connect()

    def test_descriptor_type(self):
        generic_descriptor = quippy.descriptors_module.descriptor('distance_2b cutoff=4.0')
        self.assertIsInstance(generic_descriptor, quippy.descriptors_module.descriptor, 'Wrong type defined by generic '
                                                                                        'descriptor initialiser')
        # gives different if initialised directly, so let's check that too
        distance_2b = quippy.descriptors_module.distance_2b("distance_2b cutoff=4.0")
        self.assertIsInstance(distance_2b, quippy.descriptors_module.distance_2b, 'Wrong type defined by direct '
                                                                                  'descriptor initialiser')


if __name__ == '__main__':
    unittest.main()
