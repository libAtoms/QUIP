# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright Tamas K. Stenczel, 2020
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

import ase.build
import quippy
import quippytest


class TestCalculatorSumPotential(quippytest.QuippyTestCase):
    def setUp(self):
        self.pot1 = quippy.potential.Potential("IP glue", param_filename="glue.xml")
        self.pot2 = quippy.potential.Potential("IP GAP", param_filename="GAP.xml")

        self.sumpot = quippy.potential.Potential(args_str="Potential Sum", pot1=self.pot1, pot2=self.pot2)

        # atoms objects
        self.at = ase.build.bulk('HO', 'zincblende', 3., 3., 3.)

    def calcboth(self, attribute):
        attribute = "get_{}".format(attribute)
        return getattr(self.sumpot, attribute)(self.at), \
               getattr(self.pot1, attribute)(self.at) + getattr(self.pot2, attribute)(self.at)

    def test_energy(self):
        # self.assertAlmostEqual(*self.calcboth("energy"))
        self.assertAlmostEqual(self.sumpot.get_potential_energy(self.at),
                               self.pot1.get_potential_energy(self.at) + self.pot2.get_potential_energy(self.at))

    def test_forces(self):
        self.assertArrayAlmostEqual(*self.calcboth("forces"), tol=1E-06)

    def test_virial(self):
        self.assertArrayAlmostEqual(*self.calcboth("virial"), tol=1E-06)

    def test_stress(self):
        self.assertArrayAlmostEqual(*self.calcboth("stress"), tol=1E-06)

    def test_local_virial(self):
        self.assertArrayAlmostEqual(*self.calcboth("local_virial"), tol=1E-06)

    def test_local_energy(self):
        self.assertArrayAlmostEqual(*self.calcboth("local_energy"), tol=1E-06)

    def test_energies(self):
        self.assertArrayAlmostEqual(*self.calcboth("energies"), tol=1E-06)


if __name__ == '__main__':
    unittest.main()
