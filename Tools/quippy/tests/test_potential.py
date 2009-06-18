from quippy import *
import unittest

class TestPotential(unittest.TestCase):

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

   def assertArrayAlmostEqual(self, a, b, tol=1e-8):
      self.assert_(all((a - b) < tol))
      
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
      self.pot.calc(self.at, args_str="calc_force=T")
      self.assertArrayAlmostEqual(self.at.force, self.f_ref)

   def testcalc18(self):
      self.pot.calc(self.at, args_str="calc_local_e=T")
      self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)
     
   def testcalc19(self):
      self.pot.calc(self.at, args_str="calc_local_e=T calc_force=T")
      self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)
      self.assertArrayAlmostEqual(self.at.force, self.f_ref)

   def testcalc20(self):
      self.pot.calc(self.at, args_str="calc_energy=T")
      self.assertAlmostEqual(self.at.energy, self.e_ref)
    
   def testcalc21(self):
      self.pot.calc(self.at, args_str="calc_energy=T calc_force=T")
      self.assertAlmostEqual(self.at.energy, self.e_ref)
      self.assertArrayAlmostEqual(self.at.force, self.f_ref)

   def testcalc22(self):
      self.pot.calc(self.at, args_str="calc_energy=T calc_local_e=T")
      self.assertAlmostEqual(self.at.energy, self.e_ref)
      self.assertArrayAlmostEqual(self.at.local_e, self.le_ref)

   def testcalc23(self):
      self.pot.calc(self.at, args_str="calc_energy=T calc_local_e=T calc_force=T")
      self.assertAlmostEqual(self.at.energy, self.e_ref)
      self.assertArrayAlmostEqual(self.at.force, self.f_ref)

   def testcalc_df(self):
      df = fzeros((3,self.at.n))
      self.pot.calc(self.at, df=df)
      self.assertArrayAlmostEqual(df, self.f_ref)

   def testcalc_df2(self):
      self.pot.calc(self.at, args_str="calc_df=T")
      self.assertArrayAlmostEqual(self.at.df, self.f_ref)
      

class TestMetaPotential(unittest.TestCase):

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
      

   def assertArrayAlmostEqual(self, a, b, tol=1e-8):
      self.assert_(all((a - b) < tol))

   def testcalc(self):
      v1 = fzeros((3,3))
      v2 = fzeros((3,3))
      
      self.at.calc_connect()
      self.pot.calc(self.at, args_str="calc_energy=T calc_force=T calc_df=T", virial=v1)

      cp = self.at.copy()
      self.metapot.calc(self.at, args_str="calc_energy=T calc_force=T calc_df=T", virial=v2)

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

   def testelastic_constants(self):
      self.metapot.calc_elastic_constants(self.at, c=self.c, c0=self.c0, relax_initial=True)
      self.assertArrayAlmostEqual(self.c, self.c_ref)
      self.assertArrayAlmostEqual(self.c0, self.c0_ref)



      
def getTestSuite():
   return unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(TestPotential),
                              unittest.TestLoader().loadTestsFromTestCase(TestMetaPotential)])

if __name__ == '__main__':
   suite = getTestSuite()
   unittest.TextTestRunner(verbosity=2).run(suite)

