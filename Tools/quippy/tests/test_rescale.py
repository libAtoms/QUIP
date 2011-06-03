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

class Test_Tersoff_Rescale_Energy(QuippyTestCase):

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
    

class Test_Tersoff_Rescale_Space(QuippyTestCase):

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


class Test_Tersoff_Rescale_Space_and_Energy(QuippyTestCase):

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

class Test_Tersoff_Rescale_Automatic_Factors(QuippyTestCase):
    
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


class Test_Tersoff_Rescale_Energy(QuippyTestCase):

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
    

class Test_SW_Rescale_Space(QuippyTestCase):

    def setUp(self):
        self.p1 = Potential('IP SW', param_str=quip_xml_parameters('SW'))
        self.r_scale = 0.99
        self.bulk = diamond(5.44, 14)
        self.bulk.set_cutoff(self.p1.cutoff()+1.0)
        self.bulk.calc_connect()

        self.d1 = self.bulk.copy()
        self.d1.pos[:,:] += numpy.random.uniform(-0.5, 0.5, size=3*self.d1.n).reshape(3,self.d1.n)
        self.d1.calc_connect()
        self.d2 = self.d1.copy()
        self.d2.set_lattice(self.d2.lattice/self.r_scale, True)
        self.p2 = Potential('IP SW do_rescale_r r_scale=%f' % self.r_scale,
                            param_str=quip_xml_parameters('SW'))

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
        self.assertArrayAlmostEqual([b2, v2], [b1*self.r_scale**3, v1/self.r_scale**3], tol=1e-2)


class Test_SW_Rescale_Space_and_Energy(QuippyTestCase):

    def setUp(self):
        self.p1 = Potential('IP SW', param_str=quip_xml_parameters('SW'))
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

        self.p2 = Potential('IP SW do_rescale_r r_scale=%f do_rescale_E E_scale=%f' % (self.r_scale, self.E_scale),
                            param_str=quip_xml_parameters('SW'))

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
        self.assertArrayAlmostEqual([b2, v2], [b1*self.E_scale*self.r_scale**3, v1/self.r_scale**3], tol=1e-2)

class Test_SW_Rescale_Automatic_Factors(QuippyTestCase):
    
    def setUp(self):
        self.target_vol = 165.0
        self.target_B = 105.0
        self.bulk = diamond(5.44, 14)
        verbosity_push_decrement()
        self.p = Potential('IP SW do_rescale_r do_rescale_E target_vol=%f target_B=%f minimise_bulk' % (self.target_vol, self.target_B),
                           bulk_scale=self.bulk, param_str=quip_xml_parameters('SW'))
        verbosity_pop()
        self.bulk.set_cutoff(self.p.cutoff()+1.0)
        self.bulk.calc_connect()        

    def test_factors(self):
        self.assertArrayAlmostEqual([self.p.r_scale, self.p.e_scale], [0.990180734136, 1.06610169294])

    def test_volume(self):
        verbosity_push_decrement()
        self.p.minim(self.bulk, 'cg', 1e-6, 100, do_lat=True, do_pos=True)
        verbosity_pop()
        self.assertArrayAlmostEqual([self.bulk.cell_volume()], [self.target_vol], tol=1e-4)

    def test_bulk_modulus(self):
        b, v = self.p.bulk_modulus(self.bulk, minimise_bulk=True)
        self.assertArrayAlmostEqual([b], [self.target_B], tol=1e-3)


class TestRescale_ForceMixing(QuippyTestCase):

      def setUp(self):
         dia = diamond(5.44, 14)
         xml = '<quip_params>%s</quip_params>' % (quip_xml_parameters('SW') + quip_xml_parameters('Tersoff'))

         verbosity_push(PRINT_SILENT)
         self.pot = Potential('ForceMixing init_args_pot1={IP SW} init_args_pot2={IP Tersoff}', param_str=xml)
         self.pot_r = Potential('ForceMixing init_args_pot1={IP SW} init_args_pot2={IP Tersoff} do_rescale_r minimise_bulk', param_str=xml, bulk_scale=dia)
         verbosity_pop()

         self.at = supercell(dia, 4, 4, 4)
         randomise(self.at.pos, 0.1)
         self.at.set_cutoff(self.pot.cutoff()+1.0)
         self.at.calc_connect()

         self.at.add_property('hybrid', 0)

         self.embedlist = Table(4,0,0,0)
         self.embedlist.append((1,0,0,0))
         self.at.bfs_grow_list(self.embedlist, 2, nneighb_only=True)
         self.at.hybrid[self.embedlist.int[1,:]] = 1

         self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
         self.at.hybrid_mark[self.at.hybrid != 0] = HYBRID_ACTIVE_MARK

         # Ratio of Tersoff to SW lattice constants
         self.r_scale = 5.4320052041/5.4309497787

      def tearDown(self):
         if os.path.exists('clusters.xyz'): os.unlink('clusters.xyz')
         if os.path.exists('clusters.xyz.idx'): os.unlink('clusters.xyz.idx')

      def test_cluster(self):
         f = fzeros((3,self.at.n))
         self.pot.calc(self.at, force=f, args_str="qm_args_str={single_cluster cluster_calc_connect print_clusters terminate=F randomise_buffer=F}")
         self.pot_r.calc(self.at, force=f, args_str="qm_args_str={single_cluster cluster_calc_connect print_clusters terminate=F randomise_buffer=F}")

         # Load both clusters
         cluster, cluster_r = AtomsList('clusters.xyz')

         # Check non-zero lattice components and positions
         mask = cluster.lattice != 0.0
         lattice_ratio = cluster_r.lattice[mask]/cluster.lattice[mask]

         mask = cluster.pos != 0.0
         pos_ratio = cluster_r.pos[mask]/cluster.pos[mask]

         self.assert_(abs(lattice_ratio - self.r_scale).max() < 1.0e-5 and 
                      abs(pos_ratio - self.r_scale).max() < 1.0e-5)

      def test_force_cluster_atom_mask_no_rescale(self):
          cluster_force = fzeros((3,self.at.n))
          atom_mask_force = fzeros((3,self.at.n))

          self.pot.calc(self.at, force=cluster_force, args_str="qm_args_str={single_cluster cluster_calc_connect terminate=F randomise_buffer=F}")
          self.pot.calc(self.at, force=atom_mask_force, args_str="qm_args_str={atom_mask=active_plus_cutoff}")

          self.assertArrayAlmostEqual(cluster_force, atom_mask_force)

      def test_force_cluster_atom_mask_rescale(self):
          cluster_force = fzeros((3,self.at.n))
          atom_mask_force = fzeros((3,self.at.n))

          self.pot_r.calc(self.at, force=cluster_force, args_str="qm_args_str={single_cluster cluster_calc_connect terminate=F randomise_buffer=F}")
          self.pot_r.calc(self.at, force=atom_mask_force, args_str="qm_args_str={atom_mask=active_plus_cutoff}")

          self.assertArrayAlmostEqual(cluster_force, atom_mask_force)


if __name__ == '__main__':
   unittest.main()
