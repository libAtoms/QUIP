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
import quippy
import numpy as np
import quippytest
import ase
import ase.io

class TestCalculator_GAP_Potential(quippytest.QuippyTestCase):
    def setUp(self):
        self.pot_calculator = quippy.potential.potential("IP GAP", param_filename="GAP.xml")
        self.at_orig = ase.io.read('gap_sample.xyz')

        self.at = ase.Atoms(numbers=self.at_orig.arrays['numbers'], positions=self.at_orig.get_positions(), pbc=True,
                            cell=self.at_orig.get_cell())

        self.f = np.zeros((3, len(self.at)), order='F')

        self.energy_ref = self.at_orig.info['energy']
        self.forces_ref = self.at_orig.arrays['force']

        self.at.set_calculator(self.pot_calculator)

    def test_energy(self):
        self.assertAlmostEqual(self.at.get_potential_energy(), self.energy_ref, delta=1E-05)

    def test_forces(self):
        self.assertArrayAlmostEqual(self.at.get_forces(), self.forces_ref, tol=1E-06)

if __name__ == '__main__':
    unittest.main()
