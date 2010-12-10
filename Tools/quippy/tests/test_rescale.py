# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
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

from quippy import *
import unittest, quippy
from quippytest import *

class TestRescale_Energy(QuippyTestCase):

    def setUp(self):
        self.p1 = Potential('IP Tersoff', param_str=quip_xml_parameters('Tersoff'))
        self.E_scale = 1.02
        self.bulk = diamond(5.44, 14)
        self.bulk.set_cutoff(self.p1.cutoff()+1.0)
        self.bulk.calc_connect()

        self.d1 = self.bulk.copy()
        self.d1.pos[:,:] += numpy.random.uniform(-0.5, 0.5, size=3*self.d1.n).reshape(3,self.d1.n)
        self.d1.calc_connect()
        self.d2 = self.d1.copy()

        self.p2 = Potential('IP Tersoff do_rescale_E E_scale=%f' % self.E_scale,
                            param_str=quip_xml_parameters('Tersoff'))

        self.p1.calc(self.d1, args_str="force energy virial local_virial local_energy")
        self.p2.calc(self.d2, args_str="force energy virial local_virial local_energy")

    def test_force(self):
        self.assertArrayAlmostEqual(self.d1.force*self.E_scale, self.d2.force)

    def test_energy(self):
        self.assertAlmostEqual(self.d1.energy*self.E_scale, self.d2.energy)

    def test_virial(self):
        self.assertArrayAlmostEqual(self.d1.virial*self.E_scale, self.d2.virial)
        
    def test_local_virial(self):
        self.assertArrayAlmostEqual(self.d1.local_virial*self.E_scale, self.d2.local_virial)

    def test_local_energy(self):
        self.assertArrayAlmostEqual(self.d1.local_energy*self.E_scale, self.d2.local_energy)

    def test_minim(self):
        bulk1 = self.bulk.copy()
        bulk2 = self.bulk.copy()
        verbosity_push_decrement()
        self.p1.minim(bulk1, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        self.p2.minim(bulk2, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        verbosity_pop()
        self.assertArrayAlmostEqual(bulk1.lattice, bulk2.lattice, tol=1e-4)
        self.assertArrayAlmostEqual(bulk1.pos, bulk2.pos, tol=1e-4)

    def test_bulk_modulus(self):
        b, v = self.p1.bulk_modulus(self.bulk, minimise_bulk=True)
        b2, v2 = self.p2.bulk_modulus(self.bulk, minimise_bulk=True)
        self.assertArrayAlmostEqual([b*self.E_scale, v], [b2, v2], tol=1e-4)
    

class TestRescale_Space(QuippyTestCase):

    def setUp(self):
        self.p1 = Potential('IP Tersoff', param_str=quip_xml_parameters('Tersoff'))
        self.r_scale = 0.99
        self.bulk = diamond(5.44, 14)
        self.bulk.set_cutoff(self.p1.cutoff()+1.0)
        self.bulk.calc_connect()

        self.d1 = self.bulk.copy()
        self.d1.pos[:,:] += numpy.random.uniform(-0.5, 0.5, size=3*self.d1.n).reshape(3,self.d1.n)
        self.d1.calc_connect()
        self.d2 = self.d1.copy()
        self.d2.set_lattice(self.d2.lattice/self.r_scale, True)
        self.p2 = Potential('IP Tersoff do_rescale_r r_scale=%f' % self.r_scale,
                            param_str=quip_xml_parameters('Tersoff'))

        self.p1.calc(self.d1, args_str="force energy virial local_virial local_energy")
        self.p2.calc(self.d2, args_str="force energy virial local_virial local_energy")

    def test_force(self):
        self.assertArrayAlmostEqual(self.d2.force, self.d1.force*self.r_scale)

    def test_energy(self):
        self.assertAlmostEqual(self.d2.energy, self.d1.energy)

    def test_virial(self):
        self.assertArrayAlmostEqual(self.d2.virial, self.d1.virial)
        
    def test_local_virial(self):
        self.assertArrayAlmostEqual(self.d2.local_virial, self.d1.local_virial)

    def test_local_energy(self):
        self.assertArrayAlmostEqual(self.d2.local_energy, self.d1.local_energy)

    def test_minim(self):
        bulk1 = self.bulk.copy()
        bulk2 = self.bulk.copy()
        verbosity_push_decrement()
        self.p1.minim(bulk1, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        self.p2.minim(bulk2, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        verbosity_pop()
        self.assertArrayAlmostEqual(bulk1.lattice/self.r_scale, bulk2.lattice, tol=1e-4)
        self.assertArrayAlmostEqual(bulk1.pos/self.r_scale, bulk2.pos, tol=1e-4)

    def test_bulk_modulus(self):
        b1, v1 = self.p1.bulk_modulus(self.bulk, minimise_bulk=True)
        b2, v2 = self.p2.bulk_modulus(self.bulk, minimise_bulk=True)
        self.assertArrayAlmostEqual([b2, v2], [b1*self.r_scale**3, v1/self.r_scale**3], tol=1e-3)


class TestRescale_Space_and_Energy(QuippyTestCase):

    def setUp(self):
        self.p1 = Potential('IP Tersoff', param_str=quip_xml_parameters('Tersoff'))
        self.r_scale = 0.99
        self.E_scale = 1.01
        self.bulk = diamond(5.44, 14)
        self.bulk.set_cutoff(self.p1.cutoff()+1.0)
        self.bulk.calc_connect()

        self.d1 = self.bulk.copy()
        self.d1.pos[:,:] += numpy.random.uniform(-0.5, 0.5, size=3*self.d1.n).reshape(3,self.d1.n)
        self.d1.calc_connect()
        self.d2 = self.d1.copy()
        self.d2.set_lattice(self.d2.lattice/self.r_scale, True)

        self.p2 = Potential('IP Tersoff do_rescale_r r_scale=%f do_rescale_E E_scale=%f' % (self.r_scale, self.E_scale),
                            param_str=quip_xml_parameters('Tersoff'))

        self.p1.calc(self.d1, args_str="force energy virial local_virial local_energy")
        self.p2.calc(self.d2, args_str="force energy virial local_virial local_energy")

    def test_force(self):
        self.assertArrayAlmostEqual(self.d2.force, self.d1.force*self.r_scale*self.E_scale)

    def test_energy(self):
        self.assertAlmostEqual(self.d2.energy, self.d1.energy*self.E_scale)

    def test_virial(self):
        self.assertArrayAlmostEqual(self.d2.virial, self.d1.virial*self.E_scale)
        
    def test_local_virial(self):
        self.assertArrayAlmostEqual(self.d2.local_virial, self.d1.local_virial*self.E_scale)

    def test_local_energy(self):
        self.assertArrayAlmostEqual(self.d2.local_energy, self.d1.local_energy*self.E_scale)

    def test_minim(self):
        bulk1 = self.bulk.copy()
        bulk2 = self.bulk.copy()
        verbosity_push_decrement()
        self.p1.minim(bulk1, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        self.p2.minim(bulk2, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        verbosity_pop()
        self.assertArrayAlmostEqual(bulk1.lattice/self.r_scale, bulk2.lattice, tol=1e-4)
        self.assertArrayAlmostEqual(bulk1.pos/self.r_scale, bulk2.pos, tol=1e-4)

    def test_bulk_modulus(self):
        b1, v1 = self.p1.bulk_modulus(self.bulk, minimise_bulk=True)
        b2, v2 = self.p2.bulk_modulus(self.bulk, minimise_bulk=True)
        self.assertArrayAlmostEqual([b2, v2], [b1*self.E_scale*self.r_scale**3, v1/self.r_scale**3], tol=1e-3)

class TestRescale_Automatic_Factors(QuippyTestCase):
    
    def setUp(self):
        self.target_vol = 165.0
        self.target_B = 105.0
        self.bulk = diamond(5.44, 14)
        verbosity_push_decrement()
        self.p = Potential('IP Tersoff do_rescale_r do_rescale_E target_vol=%f target_B=%f minimise_bulk' % (self.target_vol, self.target_B),
                           bulk_scale=self.bulk, param_str=quip_xml_parameters('Tersoff'))
        verbosity_pop()
        self.bulk.set_cutoff(self.p.cutoff()+1.0)
        self.bulk.calc_connect()        

    def test_factors(self):
        self.assertArrayAlmostEqual([self.p.r_scale, self.p.e_scale], [0.990373161231, 1.10443023495])

    def test_volume(self):
        verbosity_push_decrement()
        self.p.minim(self.bulk, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        verbosity_pop()
        self.assertArrayAlmostEqual([self.bulk.cell_volume()], [self.target_vol], tol=1e-4)

    def test_bulk_modulus(self):
        b, v = self.p.bulk_modulus(self.bulk, minimise_bulk=True)
        self.assertArrayAlmostEqual([b], [self.target_B], tol=1e-3)


if __name__ == '__main__':
   unittest.main()
