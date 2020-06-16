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
import os

import quippy
import quippytest
import ase
import numpy as np


@unittest.skipIf(os.environ['HAVE_GAP'] != '1', 'GAP support not enabled')
class Test_Potential_Cell(quippytest.QuippyTestCase):
    def setUp(self):
        cell_sizes = np.linspace(2.5, 4.5, 5)
        self.at_list = []

        for c in cell_sizes:
            self.at_list.append(ase.Atoms('HH', positions=[[0., 0., 0.], [1., 1., 1.]], cell=[c, c, c], pbc=True))

        self.pot = quippy.potential.Potential('', param_filename='GAP.xml')

        # calculated in a notebook with teh correct code
        self.ref_energies = [0.36747083829015637,
                             2.8715032700273735,
                             4.10632306979403,
                             5.518256035535996,
                             5.885656871424537]

    def test_pot_calculation_with_changing_cell_size(self):
        # for bugfix on 13 Aug 2019, tks32
        ener_quippy = []
        for at in self.at_list:
            at.set_calculator(self.pot)
            ener_quippy.append(at.get_potential_energy())

        self.assertArrayAlmostEqual(ener_quippy, self.ref_energies, tol=1E-06)

if __name__ == '__main__':
    unittest.main()
