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
         self.pot = Potential('IP SW', param_str=xml)

         self.at = diamond(5.44, 14)
         matrix_randomise(self.at.pos, 0.1)
         self.at.set_cutoff(self.pot.cutoff())
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
         self.pot.calc(self.at, force=self.f)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc3(self):
         self.pot.calc(self.at, local_energy=self.le)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc4(self):
         self.pot.calc(self.at, local_energy=self.le, force=self.f)
         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc5(self):
         self.pot.calc(self.at, energy=self.e)
         self.assertAlmostEqual(self.e, self.e_ref)

      def testcalc6(self):
         self.pot.calc(self.at, energy=self.e, force=self.f)
         self.assertAlmostEqual(e, self.e_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc7(self):
         self.pot.calc(self.at, energy=self.e, local_energy=self.le)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc8(self):
         self.pot.calc(self.at, energy=self.e, local_energy=self.le, force=self.f)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc9(self):
         self.pot.calc(self.at, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)

      def testcalc10(self):
         self.pot.calc(self.at, force=self.f, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc11(self):
         self.pot.calc(self.at, local_energy=self.le, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc12(self):
         self.pot.calc(self.at, local_energy=self.le, force=self.f, virial=self.v)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc13(self):
         self.pot.calc(self.at, energy=self.e, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)

      def testcalc14(self):
         self.pot.calc(self.at, energy=self.e, force=self.f, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc15(self):
         self.pot.calc(self.at, energy=self.e, local_energy=self.le, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)

      def testcalc16(self):
         self.pot.calc(self.at, energy=self.e, local_energy=self.le, force=self.f, virial=self.v)
         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.le, self.le_ref)
         self.assertArrayAlmostEqual(self.f, self.f_ref)

      def testcalc17(self):
         self.pot.calc(self.at, args_str="force")
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc18(self):
         self.pot.calc(self.at, args_str="local_energy")
         self.assertArrayAlmostEqual(self.at.local_energy, self.le_ref)

      def testcalc19(self):
         self.pot.calc(self.at, args_str="local_energy force")
         self.assertArrayAlmostEqual(self.at.local_energy, self.le_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc20(self):
         self.pot.calc(self.at, args_str="energy")
         self.assertAlmostEqual(self.at.energy, self.e_ref)

      def testcalc21(self):
         self.pot.calc(self.at, args_str={'energy':True, 'force':True})
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc22(self):
         self.pot.calc(self.at, args_str={'energy':True, 'local_energy':True})
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.local_energy, self.le_ref)

      def testcalc23(self):
         self.pot.calc(self.at, energy=True, local_energy=True, force=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)
         self.assertArrayAlmostEqual(self.at.local_energy, self.le_ref)

      def testcalc24(self):
         self.pot.calc(self.at, args_str={'energy':True, 'local_energy':True, 'force':True, 'virial':True})
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)
         self.assertArrayAlmostEqual(self.at.virial, self.v_ref)

      def testcalc25(self):
         self.pot.calc(self.at, args_str="energy local_energy force virial")
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)
         self.assertArrayAlmostEqual(self.at.virial, self.v_ref)

      def testcalc26(self):
         self.pot.calc(self.at, energy=True, local_energy=True, force=True, virial=True)
         self.assertAlmostEqual(self.at.energy, self.e_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)
         self.assertArrayAlmostEqual(self.at.virial, self.v_ref)

      def testcalc_df(self):
         self.pot.calc(self.at, args_str="force fd_force")
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      def testcalc_force_both(self):
         self.pot.calc(self.at, force=self.f, args_str="force")
         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

      # The test below is no longer relevant -- after Atoms data structure change from
      # Table to dictionary, all properties *are* contigous in memory.
      
      #def testcalc_force_non_contigous(self):
      #   self.at.add_property('force', 0.0, n_cols=3)
      #   self.assertRaises(ValueError, self.pot.calc, self.at, force=self.at.force)

      def testcalc_all(self):
         self.pot.calc(self.at, force=self.f, energy=self.e, local_energy=self.le, virial=self.v,
                       args_str="force energy local_energy virial")

         self.assertArrayAlmostEqual(self.f, self.f_ref)
         self.assertArrayAlmostEqual(self.at.force, self.f_ref)

         self.assertArrayAlmostEqual(self.le, self.le_ref)
         self.assertArrayAlmostEqual(self.at.local_energy, self.le_ref)

         self.assertAlmostEqual(self.e, self.e_ref)
         self.assertAlmostEqual(self.at.energy, self.e_ref)

         self.assertArrayAlmostEqual(self.v, self.v_ref)
         self.assertArrayAlmostEqual(self.at.virial, self.v_ref)

   got_tight_binding = True
   tight_binding_xml = """<eval_test_params>

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
       <abcd>  -0.0532334619024000
       -0.9076427431860000 -8.8308491367399995  56.5661321469000029
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

   try:
      p = Potential('TB NRL-TB', param_str=tight_binding_xml)
   except RuntimeError:
      got_tight_binding = False

   if got_tight_binding:
      class TestPotential_NRL_TB(QuippyTestCase):

         def setUp(self):
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

            self.pot = Potential('TB NRL-TB', param_str=tight_binding_xml)
            self.at = Atoms(xyz, format='string')

            verbosity_push(PRINT_SILENT)

         def tearDown(self):
            verbosity_pop()

         def test_energy(self):
            self.pot.calc(self.at, args_str="energy SCF_LOCAL_U=T")
            self.assertAlmostEqual(self.at.energy, 9.3810611711946219)

         def test_virial(self):
            self.pot.calc(self.at, SCF_LOCAL_U=True, args_str="virial")
            self.assertArrayAlmostEqual(self.at.virial,
                                        farray([[ 0.5865603956735974E+00, -0.1377850228522192E+00, -0.3055872087626610E+00],
                                                [-0.1377850228522220E+00,  0.2216201633362491E+00,  0.1910908299371670E+01],
                                                [-0.3055872087626552E+00,  0.1910908299371667E+01,  0.1366230430253259E+01]]))

         def test_force(self):
            self.pot.calc(self.at, SCF_LOCAL_U=True, args_str="force")
            self.assertArrayAlmostEqual(self.at.force,
                                        farray([[-1.77614476,    1.89628018,    0.55837291],
                                                [0.12109842,   -1.82740015,    0.52697574],
                                                [1.05222007,   -1.42686527,   -0.92338994],
                                                [-1.67217435,    1.53461226,    0.24978758],
                                                [-0.60877230,    1.25525465,    0.18449930], 
                                                [2.00297551,   -2.02323979,  -0.56112465], 
                                                [2.33071967,   -0.46343812,   -0.82128532], 
                                                [ -1.44992225,    1.05479624,    0.78616438]]))


   class TestPotential_Callback(QuippyTestCase):

      @staticmethod
      def callback_1(at):
         if at.calc_energy:
            at.params['energy'] = 1.0
         if at.calc_force:
            at.add_property('force', 0.0, n_cols=3)
            at.force[:] = 1.0
         if at.calc_virial:
            virial = fzeros((3,3))
            virial[:] = 1.0
            at.params['virial'] = virial

      @staticmethod
      def callback_2(at):
         if at.calc_energy:
            at.params['energy'] = 2.0
         if at.calc_force:
            at.add_property('force', 0.0, n_cols=3)
            at.force[:] = 2.0
         if at.calc_virial:
            virial = fzeros((3,3))
            virial[:] = 2.0
            at.params['virial'] = virial

      def setUp(self):
         self.p = Potential('CallbackPot')
         self.p.set_callback(TestPotential_Callback.callback_1)

         self.a = diamond(5.44, 14)
         self.a.set_cutoff(self.p.cutoff())
         self.a.calc_connect()

      def test_energy(self):
         self.p.calc(self.a, args_str="energy")
         self.assertAlmostEqual(self.a.energy, 1.0)

      def test_force(self):
         self.p.calc(self.a, args_str="force")
         self.assertArrayAlmostEqual(self.a.force, 1.0*numpy.ones((3,self.a.n)))

      def test_virial(self):
         self.p.calc(self.a, args_str="virial")
         self.assertArrayAlmostEqual(self.a.virial, 1.0*numpy.ones((3,3)))

      def test_all(self):
         e = farray(0.0)
         f = fzeros((3,self.a.n))
         v = fzeros((3,3))
         self.p.calc(self.a, energy=e, force=f, virial=v)
         self.assertAlmostEqual(e, 1.0)
         self.assertArrayAlmostEqual(v, 1.0*numpy.ones((3,3)))
         self.assertArrayAlmostEqual(f, 1.0*numpy.ones((3,self.a.n)))

      def test_new_pot(self):
         p2 = Potential('CallbackPot')
         p2.set_callback(TestPotential_Callback.callback_2)
         p2.calc(self.a, energy=True)
         self.assertEqual(self.a.energy, 2.0)

         self.p.calc(self.a, energy=True)
         self.assertEqual(self.a.energy, 1.0)

      def test_bad_label(self):
         p3 = Potential('CallbackPot label=BAD')
         verbosity_push(PRINT_SILENT) # suppress warning about energy not being computed
         self.assertRaises(ValueError, p3.calc, self.a, args_str="energy")
         verbosity_pop()
         

   class TestPotential_ElasticMinim(QuippyTestCase):
      # Old TestMetaPotential_Simple class migrated from test_metapot.py

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
         # self.pot = Potential('IP SW', xml)
         # self.pot = Potential'Simple', self.pot)
         self.pot = Potential('IP SW', param_str=xml)
         self.at = diamond(5.44, 14)
         matrix_randomise(self.at.pos, 0.1)
         self.at.set_cutoff(self.pot.cutoff())
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

         self.pot.calc(self.at, virial=v1, args_str="energy force")

         cp = self.at.copy()
         self.pot.calc(self.at, virial=v2, args_str="energy force")

         self.assertAlmostEqual(self.at.energy, cp.energy)
         self.assertArrayAlmostEqual(self.at.force, cp.force)
         self.assertArrayAlmostEqual(v1, v2)

      def testminim(self):
         verbosity_push(PRINT_SILENT)
         nsteps = self.pot.minim(self.at, 'cg', 1e-3, 100, do_pos=True, do_lat=True)
         verbosity_pop()

         self.assertEqual(nsteps, self.nsteps_ref)
         self.assertArrayAlmostEqual(self.at.lattice, self.lat_ref)
	 self.at.map_into_cell() # self.pos_ref is already mapped
         self.assertArrayAlmostEqual(self.at.pos, self.pos_ref)

      def testelastic_constants(self):
         self.pot.calc_elastic_constants(self.at, fd=1e-3, c=self.c, c0=self.c0, relax_initial=True)
         self.assertArrayAlmostEqual(self.c, self.c_ref)
         self.assertArrayAlmostEqual(self.c0, self.c0_ref)

         


if __name__ == '__main__':
   unittest.main()
