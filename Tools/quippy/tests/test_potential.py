from quippy import *
import unittest, quippy
from quippytest import *

if hasattr(quippy, 'Potential'):

   class TestPotential_SW(QuippyTestCase):

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

      def testcalc_force_non_contigous(self):
         self.at.add_property('force', 0.0, n_cols=3)
         self.assertRaises(ValueError, self.pot.calc, self.at, f=self.at.force)

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


   class TestPotential_NRL_TB(QuippyTestCase):

      def setUp(self):
         xml = """<eval_test_params>

<self_consistency tolerance="1e-8">
<U Z="14" U="2.0"/>
</self_consistency>

<KPoints N="4">
<point weight="2"> -0.25 -0.25 -0.25 </point>
<point weight="2"> -0.25 -0.25 0.25 </point>
<point weight="2"> -0.25 0.25 -0.25 </point>
<point weight="2"> -0.25 0.25 0.25 </point>
</KPoints>

<NRL_TB_params header_str="NN00001" label="Silicon">
<!-- from http://cst-www.nrl.navy.mil/bind/si.html, Phys. Rev. B v. 62, p. 4477 (2000), sp basis -->
  <defaults fermi_T="0.01"/>
  <header is_orthogonal="F" is_magnetic="F" has_pair_repulsion="F"
  overlap_zero_limit="T" force_harrison_signs="F"/>
  <n_types v="1"/>
  <cutoff v="12.5000000000000000"/>
  <per_type_data type="1" atomic_num="14" atomic_mass="28.0859999999999985"
    n_orbs="4" n_elecs="4" n_orb_sets="2" lambda_sq="1.2178584715159753">
    <orb_set_type>1 2   </orb_set_type>
  </per_type_data>
  <per_pair_data type1="1" type2="1" r_cut="12.5000000000000000"
    screen_l="0.5000000000000000">
    <abcd>  -0.0532334619024000  -0.9076427431860000  -8.8308491367399995  56.5661321469000029
         0.3578597152650000   0.3036476931010000   7.0922290356000000 -77.4785508398999951
        10.0000000000000000   0.0000000000000000   0.0000000000000000   0.0000000000000000
      </abcd>
    <H_coeff>  219.5608136509999895  -16.2132459617999984  -15.5048968096999999    1.5987058429226642
        10.1276876205999997  -4.4036811239600002   0.2266767833680000   0.8513235098690760
       -22.9590281075000000   1.7207707740500000   1.4191307713200001   1.0637223547809533
        10.2654492628000007   4.6718241428000002  -2.2161562720900001   1.2350950098834053
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
      </H_coeff>
    <S_coeff>  5.1575871864099998  0.6600093077760000 -0.0815441307353000  1.2279842062847821
         8.8736466648800008 -16.2407704748000015   5.1822969049500003   1.5392183365105399
        11.2504890092000007  -1.1701322928900000  -1.0591485021400000   1.2941988550186143
       -692.1842311450000125  396.1532489560000272  -13.8172106269999997    2.4727109467970014
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  0.0000000000000000  0.0000000000000000  1.0000000000000000
      </S_coeff>
  </per_pair_data>
</NRL_TB_params>
</eval_test_params>"""

         xyz = """8
Lattice="5.42883523318981 0 0 0 5.42883523318981 0 0 0 5.42883523318981" Properties=Z:I:1:pos:R:3
14 0.0210809907025043 -0.082432438809103 -0.0525403165481939
14 -0.0194362141345624 2.79575394428175 2.65958482185915
14 2.6322497362192 0.0428914772326525 2.75195107653885
14 2.7348145379646 2.61668595724592 -0.00569985609748024
14 1.33742072200028 1.26492764509989 1.32415402710619
14 1.27192747372937 4.12374105659993 4.07496695423226
14 3.97244199589577 1.36902339889138 4.0668447417454
14 4.09570476049115 4.02286216155155 1.27329051246382"""

         self.pot = Potential('TB NRL-TB', xml)
         self.at = Atoms(xyz, format='string')

         verbosity_push(SILENT)

      def tearDown(self):
         verbosity_pop()

      def test_energy(self):
         self.pot.calc(self.at, SCF_LOCAL_U=True, calc_energy=True)
         self.assertAlmostEqual(self.at.energy, 9.3810611711946219)

      def test_virial(self):
         self.pot.calc(self.at, SCF_LOCAL_U=True, calc_virial=True)
         self.assertArrayAlmostEqual(self.at.virial,
                                     farray([[ 0.5865603956735974E+00, -0.1377850228522192E+00, -0.3055872087626610E+00],
                                             [-0.1377850228522220E+00,  0.2216201633362491E+00,  0.1910908299371670E+01],
                                             [-0.3055872087626552E+00,  0.1910908299371667E+01,  0.1366230430253259E+01]]))

      def test_force(self):
         self.pot.calc(self.at, SCF_LOCAL_U=True, calc_force=True)
         self.assertArrayAlmostEqual(self.at.force,
                                     farray([[-1.77614476,    1.89628018,    0.55837291],
                                             [0.12109842,   -1.82740015,    0.52697574],
                                             [1.05222007,   -1.42686527,   -0.92338994],
                                             [-1.67217435,    1.53461226,    0.24978758],
                                             [-0.60877230,    1.25525465,    0.18449930], 
                                             [2.00297551,   -2.02323979,  -0.56112465], 
                                             [2.33071967,   -0.46343812,   -0.82128532], 
                                             [ -1.44992225,    1.05479624,    0.78616438]]))



if __name__ == '__main__':
   unittest.main()
