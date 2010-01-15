from quippy import *
import unittest, quippy
from quippytest import *

if hasattr(quippy, 'Potential'):

   class TestPotential(QuippyTestCase):

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

         self.at = diamond(5.44, 14)
         matrix_randomise(self.at.pos, 0.1)
         self.at.calc_connect()

         self.e, = fvar('e')
         self.f = fzeros((3,self.at.n))
         self.df = fzeros((3,self.at.n))
         self.le = fzeros(self.at.n)
         self.v = fzeros((3,3))

         self.e_ref = -34.5038375509

         self.le_ref = farray([-4.3144614,
                               -4.31461612,
                               -4.32476405,
                               -4.31437179,
                               -4.30909712,
                               -4.31681375,
                               -4.30561542,
                               -4.3040979 ])

         self.f_ref = farray([[ 0.89920374, -0.38025157, -0.38727027],
                              [ 0.36623356, -0.52403757,  0.7200206 ],
                              [-0.36952654,  0.12899529,  0.00458111],
                              [-0.19912365, -0.1632057 ,  1.08509495],
                              [-0.67565314, -0.59410498, -0.47921521],
                              [ 0.17097454,  0.5847822 , -0.31088749],
                              [ 0.43613712,  0.90811269,  0.1600328 ],
                              [-0.62824563,  0.03970963, -0.79235649]])

         self.v_ref = farray([[-0.34103601,  0.60925144, -0.02138795],
                              [ 0.60925144, -0.36145702, -0.19375487],
                              [-0.02138795, -0.19375487, -0.34640615]])

      def tearDown(self):
         del self.at
         del self.pot

      def testcalc(self):
         self.pot.calc(self.at)

      def testcalc2(self):
         self.pot.calc(self.at, f=self.f)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc3(self):
         self.pot.calc(self.at, local_e=self.le)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc4(self):
         self.pot.calc(self.at, local_e=self.le, f=self.f)
         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc5(self):
         self.pot.calc(self.at, e=self.e)
         self.assertAlmostEqual(self.e, self.e_ref)

      def testcalc6(self):
         self.pot.calc(self.at, e=self.e, f=self.f)
         self.assertAlmostEqual(e, self.e_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc7(self):
         self.pot.calc(self.at, e=self.e, local_e=self.le)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc8(self):
         self.pot.calc(self.at, e=self.e, local_e=self.le, f=self.f)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc9(self):
         self.pot.calc(self.at, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)

      def testcalc10(self):
         self.pot.calc(self.at, f=self.f, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc11(self):
         self.pot.calc(self.at, local_e=self.le, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc12(self):
         self.pot.calc(self.at, local_e=self.le, f=self.f, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc13(self):
         self.pot.calc(self.at, e=self.e, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)

      def testcalc14(self):
         self.pot.calc(self.at, e=self.e, f=self.f, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc15(self):
         self.pot.calc(self.at, e=self.e, local_e=self.le, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc16(self):
         self.pot.calc(self.at, e=self.e, local_e=self.le, f=self.f, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc17(self):
         self.pot.calc(self.at, calc_force=True)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc18(self):
         self.pot.calc(self.at, calc_local_e=True)
         self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)

      def testcalc19(self):
         self.pot.calc(self.at, calc_local_e=True, calc_force=True)
         self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc20(self):
         self.pot.calc(self.at, calc_energy=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)

      def testcalc21(self):
         self.pot.calc(self.at, calc_energy=True, calc_force=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc22(self):
         self.pot.calc(self.at, calc_energy=True, calc_local_e=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)

      def testcalc23(self):
         self.pot.calc(self.at, calc_energy=True, calc_local_e=True, calc_force=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc24(self):
         self.pot.calc(self.at, calc_energy=True, calc_local_e=True, calc_force=True, calc_virial=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)
         self.assertArrayAlmostEqual(self.at.virial, self.v_ref)

      def testcalc_df(self):
         self.pot.calc(self.at, df=self.df)
         self.assertArrayAlmostEqual(self.df, self.f_ref)

      def testcalc_df2(self):
         self.pot.calc(self.at, calc_df=True)
         self.assertArrayAlmostEqual(self.at.df, self.f_ref)

      def testcalc_force_both(self):
         self.pot.calc(self.at, f=self.f, calc_force=True)
         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc_all(self):
         self.pot.calc(self.at, f=self.f, df=self.df, e=self.e, local_e=self.le, virial=self.v,
                       calc_force=True, calc_energy=True, calc_local_e=True, calc_df=True, calc_virial=True)

         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

         self.assertArrayAlmostEqual(self.df, self.f_ref)
         self.assertArrayAlmostEqual(self.at.df, self.f_ref)

         self.assertArrayAlmostEqual(self.le, self.le_ref)
         self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)

         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertAlmostEqual(self.at.energy, self.e_ref)

         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.at.virial, self.v_ref)

if __name__ == '__main__':
   unittest.main()
