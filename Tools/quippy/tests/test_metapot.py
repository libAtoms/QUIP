from quippy import *
import unittest, quippy
from quippytest import *

if hasattr(quippy, 'Potential') and hasattr(quippy, 'MetaPotential'):

   class TestMetaPotentialSimple(QuippyTestCase):

      def setUp(self):

         xml="""
         <SW_params n_types="2" label="PRB_31_plus_H">
         <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984), extended for other elements </comment>
         <per_type_data type="1" atomic_num="1" />
         <per_type_data type="2" atomic_num="14" />
         <per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
         p="0" q="0" a="1.0" sigma="1.0" eps="0.0" />
         <per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
         p="4" q="0" a="1.25" sigma="2.537884" eps="2.1672" />
         <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
         p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

         <!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
         <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />

         <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         </SW_params>
         """

         system_reseed_rng(1)
         self.pot = Potential('IP SW', xml)
         self.metapot = MetaPotential('Simple', self.pot)
         self.at = diamond(5.44, 14)
         matrix_randomise(self.at.pos, 0.1)
         self.at.calc_connect()

         self.nsteps_ref = 6

         self.lat_ref = farray([[  5.43118500e+00,   6.52383282e-04,  -7.71644552e-05],
                                [  6.52353394e-04,   5.43113895e+00,  -1.05226108e-04],
                                [ -7.71047017e-05,  -1.05266395e-04,   5.43111965e+00]])

         self.pos_ref = farray([[ -1.46951886e-02,  -4.31724853e-04,  -2.98041952e-03],
                                [  1.34330077e+00,   1.35748083e+00,   1.35458260e+00],
                                [  2.70117370e+00,   2.71545008e+00,  -3.04193919e-03],
                                [ -1.37266267e+00,  -1.35840607e+00,   1.35462470e+00],
                                [  2.70085164e+00,  -1.50958027e-04,   2.71252244e+00],
                                [ -1.37227962e+00,   1.35724737e+00,  -1.36097640e+00],
                                [ -1.43610167e-02,   2.71511882e+00,   2.71250787e+00],
                                [  1.34296689e+00,  -1.35806961e+00,  -1.36098177e+00]])

         self.c = fzeros((6,6))
         self.c0 = fzeros((6,6))

         self.c_ref = farray([[  9.44769536e-01,   4.76813501e-01,   4.76813905e-01,    4.16718607e-09,   1.89783484e-08,  -4.69145249e-08],
                              [  4.76813501e-01,   9.44766827e-01,   4.76814875e-01,    1.31285521e-08,  -2.47515231e-07,   1.19759970e-07],
                              [  4.76813905e-01,   4.76814875e-01,   9.44765698e-01,   -5.57472759e-08,  -2.56955266e-07,   1.49840128e-07],
                              [  4.16718607e-09,   1.31285521e-08,  -5.57472759e-08,    3.52179803e-01,  -1.24244409e-07,   2.25322435e-07],
                              [  1.89783484e-08,  -2.47515231e-07,  -2.56955266e-07,   -1.24244409e-07,   3.52179506e-01,   6.48601834e-08],
                              [ -4.69145249e-08,   1.19759970e-07,   1.49840128e-07,    2.25322435e-07,   6.48601834e-08,   3.52179719e-01]])

         self.c0_ref = farray([[  9.44769536e-01,   4.76813501e-01,   4.76813905e-01,  -6.70046431e-08,  -7.14338269e-08,   1.84677402e-09],
                               [  4.76813501e-01,   9.44766827e-01,   4.76814875e-01,  -1.31381578e-08,  -2.47515237e-07,   1.19759950e-07],
                               [  4.76813905e-01,   4.76814875e-01,   9.44765698e-01,  -5.57472778e-08,  -2.56955271e-07,   1.49840135e-07],
                               [ -6.70046431e-08,  -1.31381578e-08,  -5.57472778e-08,   6.84792188e-01,   1.75895039e-07,  -9.05991095e-08],
                               [ -7.14338269e-08,  -2.47515237e-07,  -2.56955271e-07,   1.75895039e-07,   6.84792670e-01,  -2.66885472e-08],
                               [  1.84677402e-09,   1.19759950e-07,   1.49840135e-07,  -9.05991095e-08,  -2.66885472e-08,   6.84793194e-01]])


      def testcalc(self):
         v1 = fzeros((3,3))
         v2 = fzeros((3,3))

         self.at.calc_connect()
         self.pot.calc(self.at, calc_energy=True, calc_force=True, calc_df=True, virial=v1)

         cp = self.at.copy()
         self.metapot.calc(self.at, calc_energy=True, calc_force=True, calc_df=True, virial=v2)

         self.assertAlmostEqual(self.at.energy, cp.energy)
         self.assertArrayAlmostEqual(self.at.force, cp.force)
         self.assertArrayAlmostEqual(self.at.df, cp.df)
         self.assertArrayAlmostEqual(v1, v2)

      def testminim(self):
         verbosity_push(SILENT)
         nsteps = self.metapot.minim(self.at, 'cg', 1e-3, 100, do_pos=True, do_lat=True)
         verbosity_pop()

         self.assertEqual(nsteps, self.nsteps_ref)
         self.assertArrayAlmostEqual(self.at.lattice, self.lat_ref)
         self.assertArrayAlmostEqual(self.at.pos, self.pos_ref)

      if hasattr(MetaPotential, 'calc_elastic_constants'):
         def testelastic_constants(self):
            self.metapot.calc_elastic_constants(self.at, c=self.c, c0=self.c0, relax_initial=True)
            self.assertArrayAlmostEqual(self.c, self.c_ref)
            self.assertArrayAlmostEqual(self.c0, self.c0_ref)


   class TestMetaPotentialForceMixing(QuippyTestCase):

      def setUp(self):

         xml_sw ="""
         <SW_params n_types="2" label="PRB_31_plus_H">
         <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984), extended for other elements </comment>
         <per_type_data type="1" atomic_num="1" />
         <per_type_data type="2" atomic_num="14" />
         <per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
         p="0" q="0" a="1.0" sigma="1.0" eps="0.0" />
         <per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
         p="4" q="0" a="1.25" sigma="2.537884" eps="2.1672" />
         <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
         p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

         <!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
         <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />

         <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         </SW_params>
         """

         xml_sw_eps_23="""
         <SW_params n_types="2" label="PRB_31_plus_H">
         <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984), extended for other elements </comment>
         <per_type_data type="1" atomic_num="1" />
         <per_type_data type="2" atomic_num="14" />
         <per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
         p="0" q="0" a="1.0" sigma="1.0" eps="0.0" />
         <per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827" p="4" q="0" a="1.25" sigma="2.537884" eps="2.3" />
         <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584" p="4" q="0" a="1.80" sigma="2.0951" eps="2.3" />

         <!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
         <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />

         <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
         </SW_params>
         """

         system_reseed_rng(1)
         self.pot1 = Potential('IP SW', xml_sw)
         self.pot2 = Potential('IP SW', xml_sw_eps_23)
         self.metapot = MetaPotential('ForceMixing', pot=self.pot1, pot2=self.pot2)

         dia = diamond(5.44, 14)
         self.at = supercell(dia, 3, 3, 3)
         matrix_randomise(self.at.pos, 0.1)
         self.at.calc_connect()

      def test_force_mixing_abrupt(self):
         pass
         
         

if __name__ == '__main__':
   unittest.main()


