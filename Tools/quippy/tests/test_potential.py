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
import numpy as np
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
         randomise(self.at.pos, 0.1)
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
                              [-0.62824563,  0.03970963, -0.79235649]]).T

         self.v_ref = farray([[-0.34103601,  0.60925144, -0.02138795],
                              [ 0.60925144, -0.36145702, -0.19375487],
                              [-0.02138795, -0.19375487, -0.34640615]]).T

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
         self.pot.calc(self.at, args_str="force force_using_fd=T force_fd_delta=1.0e-4")
         self.assertArrayAlmostEqual(self.at.force, self.f_ref, tol=1e-4)

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

   NRL_TB_tight_binding_xml = """<eval_test_params>

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

   DFTB_tight_binding_xml="""<eval_test_params>

   <self_consistency tolerance="1e-8">
   <U Z="6" U="1.0"/>
   <U Z="14" U="2.0"/>
   </self_consistency>

   <KPoints N="4">
   <point weight="2"> -0.25 -0.25 -0.25 </point>
   <point weight="2"> -0.25 -0.25 0.25 </point>
   <point weight="2"> -0.25 0.25 -0.25 </point>
   <point weight="2"> -0.25 0.25 0.25 </point>
   </KPoints>

<DFTB_params label="pbc-0-1">
<!-- from www.dftb.org pbc-0-1 parameter set -->
<defaults fermi_T="0.01" fermi_e="-3.4040794422705207" band_width="11.0"/>
  <n_types v="2"/>
  <cutoff v="10.4000000000000004"/>
  <max_n_orb_sets v="2"/>
  <per_type_data type="1" atomic_num="6" n_orbs="4" n_elecs="4" n_orb_sets="2">
    <orb_set_type>1 2   </orb_set_type>
    <E> -0.5048917200000000E+00 -0.1943551100000000E+00</E>
  </per_type_data>
  <per_type_data type="2" atomic_num="14" n_orbs="4" n_elecs="4" n_orb_sets="2">
    <orb_set_type>1 2   </orb_set_type>
    <E> -0.3957250600000000E+00 -0.1503138000000000E+00</E>
  </per_type_data>
  <per_pair_data type1="1" type2="1" SK_cutoff="10.4000000000000004"
    Vrep_cutoff="3.8384700000000000" SK_npts="521" Vrep_npts="61">
    <H_spline>
      <point r="0.000">  0.38384700E+01  0.00000000E+00  0.30153300E+00 -0.17241200E+00  0.32959000E-01 -0.12889900E+00 -0.11686400E+00  0.12010000E+02  0.18033900E-01  0.94725100E-01</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400"> -0.11499539E+01 -0.32050322E+00 -0.98576278E+00 -0.13221464E+01</point>
      <point r="0.420"> -0.11022958E+01 -0.30780693E+00 -0.94454972E+00 -0.13049945E+01</point>
      <point r="0.440"> -0.10581907E+01 -0.29355043E+00 -0.90322160E+00 -0.12875240E+01</point>
      <point r="0.460"> -0.10175005E+01 -0.27790765E+00 -0.86190127E+00 -0.12697757E+01</point>
      <point r="0.480"> -0.98007289E+00 -0.26104612E+00 -0.82070548E+00 -0.12517858E+01</point>
      <point r="0.500"> -0.94574855E+00 -0.24312361E+00 -0.77973906E+00 -0.12335925E+01</point>
      <point r="0.520"> -0.91436107E+00 -0.22429030E+00 -0.73909527E+00 -0.12152310E+01</point>
      <point r="0.540"> -0.88573928E+00 -0.20468939E+00 -0.69886054E+00 -0.11967324E+01</point>
      <point r="0.560"> -0.85971411E+00 -0.18445468E+00 -0.65911682E+00 -0.11781289E+01</point>
      <point r="0.580"> -0.83611945E+00 -0.16371109E+00 -0.61993993E+00 -0.11594510E+01</point>
      <point r="0.600"> -0.81478967E+00 -0.14257577E+00 -0.58139509E+00 -0.11407255E+01</point>
      <point r="0.620"> -0.79556076E+00 -0.12115741E+00 -0.54353521E+00 -0.11219771E+01</point>
      <point r="0.640"> -0.77827385E+00 -0.99556105E-01 -0.50640522E+00 -0.11032312E+01</point>
      <point r="0.660"> -0.76277694E+00 -0.77864656E-01 -0.47004809E+00 -0.10845115E+01</point>
      <point r="0.680"> -0.74892517E+00 -0.56169258E-01 -0.43450635E+00 -0.10658392E+01</point>
      <point r="0.700"> -0.73658068E+00 -0.34548739E-01 -0.39981855E+00 -0.10472335E+01</point>
      <point r="0.720"> -0.72561219E+00 -0.13074001E-01 -0.36601491E+00 -0.10287125E+01</point>
      <point r="0.740"> -0.71589411E+00  0.81913795E-02 -0.33311570E+00 -0.10102935E+01</point>
      <point r="0.760"> -0.70730615E+00  0.29190042E-01 -0.30113285E+00 -0.99199189E+00</point>
      <point r="0.780"> -0.69973402E+00  0.49870117E-01 -0.27007323E+00 -0.97382201E+00</point>
      <point r="0.800"> -0.69307074E+00  0.70184983E-01 -0.23994209E+00 -0.95579711E+00</point>
      <point r="0.820"> -0.68721373E+00  0.90091761E-01 -0.21073908E+00 -0.93793084E+00</point>
      <point r="0.840"> -0.68208568E+00  0.10955496E+00 -0.18250287E+00 -0.92023248E+00</point>
      <point r="0.860"> -0.67759358E+00  0.12854364E+00 -0.15521223E+00 -0.90271248E+00</point>
      <point r="0.880"> -0.67365754E+00  0.14703078E+00 -0.12886262E+00 -0.88537987E+00</point>
      <point r="0.900"> -0.67020374E+00  0.16499328E+00 -0.10344770E+00 -0.86824283E+00</point>
      <point r="0.920"> -0.66716436E+00  0.18241175E+00 -0.78959978E-01 -0.85130871E+00</point>
      <point r="0.940"> -0.66447733E+00  0.19927027E+00 -0.55390940E-01 -0.83458409E+00</point>
      <point r="0.960"> -0.66208588E+00  0.21555617E+00 -0.32730685E-01 -0.81807473E+00</point>
      <point r="0.980"> -0.65993806E+00  0.23125983E+00 -0.10967413E-01 -0.80178564E+00</point>
      <point r="1.000"> -0.65798620E+00  0.24637437E+00  0.99129084E-02 -0.78572113E+00</point>
      <point r="1.020"> -0.65618662E+00  0.26089547E+00  0.29926536E-01 -0.76988488E+00</point>
      <point r="1.040"> -0.65449950E+00  0.27482107E+00  0.49091518E-01 -0.75428011E+00</point>
      <point r="1.060"> -0.65288872E+00  0.28815114E+00  0.67426962E-01 -0.73890955E+00</point>
      <point r="1.080"> -0.65132185E+00  0.30088736E+00  0.84952210E-01 -0.72377545E+00</point>
      <point r="1.100"> -0.64976996E+00  0.31303292E+00  0.10168618E+00 -0.70887959E+00</point>
      <point r="1.120"> -0.64820747E+00  0.32459231E+00  0.11764712E+00 -0.69422323E+00</point>
      <point r="1.140"> -0.64661192E+00  0.33557129E+00  0.13285270E+00 -0.67980722E+00</point>
      <point r="1.160"> -0.64496372E+00  0.34597684E+00  0.14732054E+00 -0.66563201E+00</point>
      <point r="1.180"> -0.64324582E+00  0.35581714E+00  0.16106866E+00 -0.65169774E+00</point>
      <point r="1.200"> -0.64144336E+00  0.36510144E+00  0.17411585E+00 -0.63800424E+00</point>
      <point r="1.220"> -0.63954336E+00  0.37383984E+00  0.18648179E+00 -0.62455099E+00</point>
      <point r="1.240"> -0.63753445E+00  0.38204313E+00  0.19818701E+00 -0.61133716E+00</point>
      <point r="1.260"> -0.63540667E+00  0.38972259E+00  0.20925254E+00 -0.59836156E+00</point>
      <point r="1.280"> -0.63315151E+00  0.39688992E+00  0.21969973E+00 -0.58562274E+00</point>
      <point r="1.300"> -0.63076190E+00  0.40355722E+00  0.22954984E+00 -0.57311907E+00</point>
      <point r="1.320"> -0.62823226E+00  0.40973698E+00  0.23882379E+00 -0.56084879E+00</point>
      <point r="1.340"> -0.62555854E+00  0.41544205E+00  0.24754186E+00 -0.54881003E+00</point>
      <point r="1.360"> -0.62273810E+00  0.42068556E+00  0.25572352E+00 -0.53700086E+00</point>
      <point r="1.380"> -0.61976957E+00  0.42548081E+00  0.26338733E+00 -0.52541918E+00</point>
      <point r="1.400"> -0.61665265E+00  0.42984117E+00  0.27055100E+00 -0.51406269E+00</point>
      <point r="1.420"> -0.61338789E+00  0.43378004E+00  0.27723156E+00 -0.50292889E+00</point>
      <point r="1.440"> -0.60997657E+00  0.43731076E+00  0.28344553E+00 -0.49201505E+00</point>
      <point r="1.460"> -0.60642051E+00  0.44044666E+00  0.28920924E+00 -0.48131830E+00</point>
      <point r="1.480"> -0.60272206E+00  0.44320104E+00  0.29453895E+00 -0.47083561E+00</point>
      <point r="1.500"> -0.59888402E+00  0.44558724E+00  0.29945109E+00 -0.46056393E+00</point>
      <point r="1.520"> -0.59490958E+00  0.44761858E+00  0.30396222E+00 -0.45050019E+00</point>
      <point r="1.540"> -0.59080235E+00  0.44930840E+00  0.30808908E+00 -0.44064136E+00</point>
      <point r="1.560"> -0.58656624E+00  0.45066995E+00  0.31184844E+00 -0.43098442E+00</point>
      <point r="1.580"> -0.58220545E+00  0.45171637E+00  0.31525697E+00 -0.42152634E+00</point>
      <point r="1.600"> -0.57772441E+00  0.45246057E+00  0.31833109E+00 -0.41226411E+00</point>
      <point r="1.620"> -0.57312769E+00  0.45291521E+00  0.32108681E+00 -0.40319460E+00</point>
      <point r="1.640"> -0.56842001E+00  0.45309263E+00  0.32353955E+00 -0.39431466E+00</point>
      <point r="1.660"> -0.56359888E+00  0.45300992E+00  0.32574396E+00 -0.38562125E+00</point>
      <point r="1.680"> -0.55868263E+00  0.45266800E+00  0.32763685E+00 -0.37711016E+00</point>
      <point r="1.700"> -0.55367109E+00  0.45208349E+00  0.32926498E+00 -0.36877863E+00</point>
      <point r="1.720"> -0.54856947E+00  0.45126745E+00  0.33064086E+00 -0.36062340E+00</point>
      <point r="1.740"> -0.54338295E+00  0.45023063E+00  0.33177653E+00 -0.35264124E+00</point>
      <point r="1.760"> -0.53811673E+00  0.44898343E+00  0.33268354E+00 -0.34482890E+00</point>
      <point r="1.780"> -0.53277598E+00  0.44753595E+00  0.33337296E+00 -0.33718321E+00</point>
      <point r="1.800"> -0.52736584E+00  0.44589802E+00  0.33385542E+00 -0.32970099E+00</point>
      <point r="1.820"> -0.52189144E+00  0.44407918E+00  0.33414116E+00 -0.32237908E+00</point>
      <point r="1.840"> -0.51635786E+00  0.44208870E+00  0.33424003E+00 -0.31521437E+00</point>
      <point r="1.860"> -0.51077015E+00  0.43993562E+00  0.33416157E+00 -0.30820376E+00</point>
      <point r="1.880"> -0.50513329E+00  0.43762873E+00  0.33391497E+00 -0.30134417E+00</point>
      <point r="1.900"> -0.49945219E+00  0.43517656E+00  0.33350916E+00 -0.29463256E+00</point>
      <point r="1.920"> -0.49373165E+00  0.43258740E+00  0.33295273E+00 -0.28806590E+00</point>
      <point r="1.940"> -0.48797639E+00  0.42986927E+00  0.33225398E+00 -0.28164123E+00</point>
      <point r="1.960"> -0.48219101E+00  0.42702990E+00  0.33142086E+00 -0.27535561E+00</point>
      <point r="1.980"> -0.47637994E+00  0.42407674E+00  0.33046100E+00 -0.26920615E+00</point>
      <point r="2.000"> -0.47054751E+00  0.42101695E+00  0.32938167E+00 -0.26319000E+00</point>
      <point r="2.020"> -0.46469786E+00  0.41785736E+00  0.32818975E+00 -0.25730437E+00</point>
      <point r="2.040"> -0.45883498E+00  0.41460451E+00  0.32689176E+00 -0.25154650E+00</point>
      <point r="2.060"> -0.45296267E+00  0.41126463E+00  0.32549388E+00 -0.24591371E+00</point>
      <point r="2.080"> -0.44708460E+00  0.40784364E+00  0.32400189E+00 -0.24040332E+00</point>
      <point r="2.100"> -0.44120423E+00  0.40434718E+00  0.32242127E+00 -0.23501272E+00</point>
      <point r="2.120"> -0.43532489E+00  0.40078065E+00  0.32075718E+00 -0.22973935E+00</point>
      <point r="2.140"> -0.42944978E+00  0.39714917E+00  0.31901448E+00 -0.22458066E+00</point>
      <point r="2.160"> -0.42358194E+00  0.39345766E+00  0.31719780E+00 -0.21953419E+00</point>
      <point r="2.180"> -0.41772431E+00  0.38971083E+00  0.31531153E+00 -0.21459749E+00</point>
      <point r="2.200"> -0.41187971E+00  0.38591322E+00  0.31335989E+00 -0.20976817E+00</point>
      <point r="2.220"> -0.40605087E+00  0.38206922E+00  0.31134691E+00 -0.20504390E+00</point>
      <point r="2.240"> -0.40024044E+00  0.37818307E+00  0.30927648E+00 -0.20042239E+00</point>
      <point r="2.260"> -0.39445095E+00  0.37425885E+00  0.30715234E+00 -0.19590140E+00</point>
      <point r="2.280"> -0.38868489E+00  0.37030054E+00  0.30497815E+00 -0.19147876E+00</point>
      <point r="2.300"> -0.38294463E+00  0.36631197E+00  0.30275739E+00 -0.18715231E+00</point>
      <point r="2.320"> -0.37723246E+00  0.36229685E+00  0.30049347E+00 -0.18291997E+00</point>
      <point r="2.340"> -0.37155056E+00  0.35825874E+00  0.29818964E+00 -0.17877970E+00</point>
      <point r="2.360"> -0.36590100E+00  0.35420104E+00  0.29584904E+00 -0.17472949E+00</point>
      <point r="2.380"> -0.36028576E+00  0.35012701E+00  0.29347463E+00 -0.17076738E+00</point>
      <point r="2.400"> -0.35470669E+00  0.34603975E+00  0.29106927E+00 -0.16689144E+00</point>
      <point r="2.420"> -0.34916551E+00  0.34194220E+00  0.28863564E+00 -0.16309979E+00</point>
      <point r="2.440"> -0.34366384E+00  0.33783714E+00  0.28617627E+00 -0.15939059E+00</point>
      <point r="2.460"> -0.33820318E+00  0.33372717E+00  0.28369351E+00 -0.15576201E+00</point>
      <point r="2.480"> -0.33278490E+00  0.32961476E+00  0.28118959E+00 -0.15221228E+00</point>
      <point r="2.500"> -0.32741029E+00  0.32550222E+00  0.27866656E+00 -0.14873966E+00</point>
      <point r="2.520"> -0.32208053E+00  0.32139173E+00  0.27612634E+00 -0.14534246E+00</point>
      <point r="2.540"> -0.31679670E+00  0.31728533E+00  0.27357069E+00 -0.14201901E+00</point>
      <point r="2.560"> -0.31155982E+00  0.31318494E+00  0.27100125E+00 -0.13876771E+00</point>
      <point r="2.580"> -0.30637080E+00  0.30909239E+00  0.26841955E+00 -0.13558696E+00</point>
      <point r="2.600"> -0.30123052E+00  0.30500940E+00  0.26582702E+00 -0.13247523E+00</point>
      <point r="2.620"> -0.29613976E+00  0.30093758E+00  0.26322495E+00 -0.12943102E+00</point>
      <point r="2.640"> -0.29109926E+00  0.29687849E+00  0.26061460E+00 -0.12645287E+00</point>
      <point r="2.660"> -0.28610969E+00  0.29283359E+00  0.25799713E+00 -0.12353935E+00</point>
      <point r="2.680"> -0.28117169E+00  0.28880429E+00  0.25537364E+00 -0.12068908E+00</point>
      <point r="2.700"> -0.27628583E+00  0.28479191E+00  0.25274518E+00 -0.11790069E+00</point>
      <point r="2.720"> -0.27145262E+00  0.28079773E+00  0.25011277E+00 -0.11517287E+00</point>
      <point r="2.740"> -0.26667256E+00  0.27682297E+00  0.24747738E+00 -0.11250431E+00</point>
      <point r="2.760"> -0.26194606E+00  0.27286878E+00  0.24483996E+00 -0.10989376E+00</point>
      <point r="2.780"> -0.25727351E+00  0.26893627E+00  0.24220143E+00 -0.10733997E+00</point>
      <point r="2.800"> -0.25265525E+00  0.26502649E+00  0.23956271E+00 -0.10484173E+00</point>
      <point r="2.820"> -0.24809156E+00  0.26114046E+00  0.23692467E+00 -0.10239786E+00</point>
      <point r="2.840"> -0.24358271E+00  0.25727913E+00  0.23428819E+00 -0.10000720E+00</point>
      <point r="2.860"> -0.23912890E+00  0.25344341E+00  0.23165413E+00 -0.97668598E-01</point>
      <point r="2.880"> -0.23473031E+00  0.24963416E+00  0.22902330E+00 -0.95380957E-01</point>
      <point r="2.900"> -0.23038705E+00  0.24585222E+00  0.22639655E+00 -0.93143185E-01</point>
      <point r="2.920"> -0.22609924E+00  0.24209834E+00  0.22377464E+00 -0.90954215E-01</point>
      <point r="2.940"> -0.22186694E+00  0.23837327E+00  0.22115836E+00 -0.88813009E-01</point>
      <point r="2.960"> -0.21769016E+00  0.23467770E+00  0.21854845E+00 -0.86718548E-01</point>
      <point r="2.980"> -0.21356892E+00  0.23101228E+00  0.21594563E+00 -0.84669839E-01</point>
      <point r="3.000"> -0.20950318E+00  0.22737761E+00  0.21335058E+00 -0.82665913E-01</point>
      <point r="3.020"> -0.20549289E+00  0.22377427E+00  0.21076395E+00 -0.80705821E-01</point>
      <point r="3.040"> -0.20153794E+00  0.22020279E+00  0.20818638E+00 -0.78788638E-01</point>
      <point r="3.060"> -0.19763824E+00  0.21666366E+00  0.20561845E+00 -0.76913462E-01</point>
      <point r="3.080"> -0.19379363E+00  0.21315734E+00  0.20306071E+00 -0.75079411E-01</point>
      <point r="3.100"> -0.19000397E+00  0.20968425E+00  0.20051370E+00 -0.73285623E-01</point>
      <point r="3.120"> -0.18626905E+00  0.20624477E+00  0.19797789E+00 -0.71531257E-01</point>
      <point r="3.140"> -0.18258868E+00  0.20283925E+00  0.19545375E+00 -0.69815491E-01</point>
      <point r="3.160"> -0.17896261E+00  0.19946799E+00  0.19294171E+00 -0.68137522E-01</point>
      <point r="3.180"> -0.17539058E+00  0.19613129E+00  0.19044215E+00 -0.66496564E-01</point>
      <point r="3.200"> -0.17187233E+00  0.19282939E+00  0.18795547E+00 -0.64891847E-01</point>
      <point r="3.220"> -0.16840754E+00  0.18956251E+00  0.18548198E+00 -0.63322619E-01</point>
      <point r="3.240"> -0.16499590E+00  0.18633083E+00  0.18302203E+00 -0.61788143E-01</point>
      <point r="3.260"> -0.16163708E+00  0.18313451E+00  0.18057592E+00 -0.60287698E-01</point>
      <point r="3.280"> -0.15833071E+00  0.17997369E+00  0.17814392E+00 -0.58820577E-01</point>
      <point r="3.300"> -0.15507642E+00  0.17684848E+00  0.17572631E+00 -0.57386087E-01</point>
      <point r="3.320"> -0.15187383E+00  0.17375895E+00  0.17332333E+00 -0.55983550E-01</point>
      <point r="3.340"> -0.14872442E+00  0.17070802E+00  0.17093957E+00 -0.54612339E-01</point>
      <point r="3.360"> -0.14562330E+00  0.16768869E+00  0.16856400E+00 -0.53271695E-01</point>
      <point r="3.380"> -0.14257261E+00  0.16470519E+00  0.16620381E+00 -0.51961047E-01</point>
      <point r="3.400"> -0.13957195E+00  0.16175755E+00  0.16385930E+00 -0.50679772E-01</point>
      <point r="3.420"> -0.13662089E+00  0.15884581E+00  0.16153072E+00 -0.49427261E-01</point>
      <point r="3.440"> -0.13371898E+00  0.15596998E+00  0.15921833E+00 -0.48202917E-01</point>
      <point r="3.460"> -0.13086578E+00  0.15313005E+00  0.15692239E+00 -0.47006156E-01</point>
      <point r="3.480"> -0.12806085E+00  0.15032599E+00  0.15464312E+00 -0.45836405E-01</point>
      <point r="3.500"> -0.12530372E+00  0.14755778E+00  0.15238077E+00 -0.44693105E-01</point>
      <point r="3.520"> -0.12259392E+00  0.14482535E+00  0.15013553E+00 -0.43575705E-01</point>
      <point r="3.540"> -0.11993097E+00  0.14212865E+00  0.14790762E+00 -0.42483667E-01</point>
      <point r="3.560"> -0.11731440E+00  0.13946758E+00  0.14569724E+00 -0.41416466E-01</point>
      <point r="3.580"> -0.11474372E+00  0.13684207E+00  0.14350456E+00 -0.40373584E-01</point>
      <point r="3.600"> -0.11221844E+00  0.13425199E+00  0.14132977E+00 -0.39354516E-01</point>
      <point r="3.620"> -0.10973806E+00  0.13169724E+00  0.13917303E+00 -0.38358766E-01</point>
      <point r="3.640"> -0.10730207E+00  0.12917768E+00  0.13703450E+00 -0.37385848E-01</point>
      <point r="3.660"> -0.10490999E+00  0.12669318E+00  0.13491434E+00 -0.36435287E-01</point>
      <point r="3.680"> -0.10256130E+00  0.12424358E+00  0.13281268E+00 -0.35506616E-01</point>
      <point r="3.700"> -0.10025548E+00  0.12182872E+00  0.13072966E+00 -0.34599378E-01</point>
      <point r="3.720"> -0.97992035E-01  0.11944844E+00  0.12866542E+00 -0.33713125E-01</point>
      <point r="3.740"> -0.95770441E-01  0.11710256E+00  0.12662006E+00 -0.32847419E-01</point>
      <point r="3.760"> -0.93590186E-01  0.11479088E+00  0.12459371E+00 -0.32001829E-01</point>
      <point r="3.780"> -0.91450752E-01  0.11251322E+00  0.12258648E+00 -0.31175933E-01</point>
      <point r="3.800"> -0.89351622E-01  0.11026938E+00  0.12059846E+00 -0.30369318E-01</point>
      <point r="3.820"> -0.87292277E-01  0.10805914E+00  0.11862977E+00 -0.29581579E-01</point>
      <point r="3.840"> -0.85272201E-01  0.10588229E+00  0.11668047E+00 -0.28812318E-01</point>
      <point r="3.860"> -0.83290873E-01  0.10373861E+00  0.11475067E+00 -0.28061146E-01</point>
      <point r="3.880"> -0.81347777E-01  0.10162787E+00  0.11284044E+00 -0.27327680E-01</point>
      <point r="3.900"> -0.79442395E-01  0.99549835E-01  0.11094986E+00 -0.26611547E-01</point>
      <point r="3.920"> -0.77574211E-01  0.97504268E-01  0.10907898E+00 -0.25912380E-01</point>
      <point r="3.940"> -0.75742708E-01  0.95490922E-01  0.10722788E+00 -0.25229817E-01</point>
      <point r="3.960"> -0.73947373E-01  0.93509547E-01  0.10539661E+00 -0.24563507E-01</point>
      <point r="3.980"> -0.72187693E-01  0.91559887E-01  0.10358521E+00 -0.23913103E-01</point>
      <point r="4.000"> -0.70463156E-01  0.89641682E-01  0.10179374E+00 -0.23278266E-01</point>
      <point r="4.020"> -0.68773252E-01  0.87754669E-01  0.10002223E+00 -0.22658662E-01</point>
      <point r="4.040"> -0.67117474E-01  0.85898576E-01  0.98270713E-01 -0.22053966E-01</point>
      <point r="4.060"> -0.65495316E-01  0.84073132E-01  0.96539210E-01 -0.21463856E-01</point>
      <point r="4.080"> -0.63906275E-01  0.82278059E-01  0.94827743E-01 -0.20888020E-01</point>
      <point r="4.100"> -0.62349850E-01  0.80513076E-01  0.93136324E-01 -0.20326149E-01</point>
      <point r="4.120"> -0.60825542E-01  0.78777898E-01  0.91464960E-01 -0.19777942E-01</point>
      <point r="4.140"> -0.59332857E-01  0.77072238E-01  0.89813652E-01 -0.19243102E-01</point>
      <point r="4.160"> -0.57871301E-01  0.75395803E-01  0.88182393E-01 -0.18721340E-01</point>
      <point r="4.180"> -0.56440385E-01  0.73748301E-01  0.86571171E-01 -0.18212370E-01</point>
      <point r="4.200"> -0.55039621E-01  0.72129434E-01  0.84979969E-01 -0.17715913E-01</point>
      <point r="4.220"> -0.53668528E-01  0.70538903E-01  0.83408761E-01 -0.17231695E-01</point>
      <point r="4.240"> -0.52326624E-01  0.68976406E-01  0.81857518E-01 -0.16759449E-01</point>
      <point r="4.260"> -0.51013433E-01  0.67441639E-01  0.80326204E-01 -0.16298911E-01</point>
      <point r="4.280"> -0.49728483E-01  0.65934296E-01  0.78814777E-01 -0.15849823E-01</point>
      <point r="4.300"> -0.48471303E-01  0.64454068E-01  0.77323189E-01 -0.15411932E-01</point>
      <point r="4.320"> -0.47241430E-01  0.63000647E-01  0.75851388E-01 -0.14984990E-01</point>
      <point r="4.340"> -0.46038400E-01  0.61573722E-01  0.74399317E-01 -0.14568754E-01</point>
      <point r="4.360"> -0.44861757E-01  0.60172981E-01  0.72966910E-01 -0.14162984E-01</point>
      <point r="4.380"> -0.43711047E-01  0.58798110E-01  0.71554101E-01 -0.13767448E-01</point>
      <point r="4.400"> -0.42585821E-01  0.57448795E-01  0.70160816E-01 -0.13381916E-01</point>
      <point r="4.420"> -0.41485634E-01  0.56124722E-01  0.68786977E-01 -0.13006163E-01</point>
      <point r="4.440"> -0.40410045E-01  0.54825575E-01  0.67432502E-01 -0.12639969E-01</point>
      <point r="4.460"> -0.39358618E-01  0.53551039E-01  0.66097302E-01 -0.12283118E-01</point>
      <point r="4.480"> -0.38330920E-01  0.52300798E-01  0.64781288E-01 -0.11935399E-01</point>
      <point r="4.500"> -0.37326525E-01  0.51074535E-01  0.63484364E-01 -0.11596603E-01</point>
      <point r="4.520"> -0.36345010E-01  0.49871936E-01  0.62206429E-01 -0.11266528E-01</point>
      <point r="4.540"> -0.35385956E-01  0.48692684E-01  0.60947382E-01 -0.10944973E-01</point>
      <point r="4.560"> -0.34448950E-01  0.47536465E-01  0.59707113E-01 -0.10631744E-01</point>
      <point r="4.580"> -0.33533583E-01  0.46402963E-01  0.58485513E-01 -0.10326649E-01</point>
      <point r="4.600"> -0.32639451E-01  0.45291864E-01  0.57282467E-01 -0.10029500E-01</point>
      <point r="4.620"> -0.31766155E-01  0.44202856E-01  0.56097858E-01 -0.97401134E-02</point>
      <point r="4.640"> -0.30913300E-01  0.43135626E-01  0.54931563E-01 -0.94583090E-02</point>
      <point r="4.660"> -0.30080497E-01  0.42089863E-01  0.53783459E-01 -0.91839101E-02</point>
      <point r="4.680"> -0.29267361E-01  0.41065255E-01  0.52653418E-01 -0.89167437E-02</point>
      <point r="4.700"> -0.28473512E-01  0.40061495E-01  0.51541310E-01 -0.86566402E-02</point>
      <point r="4.720"> -0.27698575E-01  0.39078274E-01  0.50447001E-01 -0.84034336E-02</point>
      <point r="4.740"> -0.26942181E-01  0.38115287E-01  0.49370355E-01 -0.81569612E-02</point>
      <point r="4.760"> -0.26203963E-01  0.37172227E-01  0.48311233E-01 -0.79170636E-02</point>
      <point r="4.780"> -0.25483563E-01  0.36248792E-01  0.47269494E-01 -0.76835847E-02</point>
      <point r="4.800"> -0.24780623E-01  0.35344680E-01  0.46244994E-01 -0.74563715E-02</point>
      <point r="4.820"> -0.24094795E-01  0.34459592E-01  0.45237586E-01 -0.72352745E-02</point>
      <point r="4.840"> -0.23425733E-01  0.33593229E-01  0.44247122E-01 -0.70201469E-02</point>
      <point r="4.860"> -0.22773095E-01  0.32745294E-01  0.43273450E-01 -0.68108451E-02</point>
      <point r="4.880"> -0.22136547E-01  0.31915494E-01  0.42316418E-01 -0.66072285E-02</point>
      <point r="4.900"> -0.21515757E-01  0.31103535E-01  0.41375871E-01 -0.64091596E-02</point>
      <point r="4.920"> -0.20910400E-01  0.30309129E-01  0.40451650E-01 -0.62165034E-02</point>
      <point r="4.940"> -0.20320154E-01  0.29531987E-01  0.39543598E-01 -0.60291280E-02</point>
      <point r="4.960"> -0.19744703E-01  0.28771822E-01  0.38651553E-01 -0.58469044E-02</point>
      <point r="4.980"> -0.19183736E-01  0.28028352E-01  0.37775353E-01 -0.56697060E-02</point>
      <point r="5.000"> -0.18636947E-01  0.27301294E-01  0.36914835E-01 -0.54974092E-02</point>
      <point r="5.020"> -0.18104032E-01  0.26590370E-01  0.36069833E-01 -0.53298929E-02</point>
      <point r="5.040"> -0.17584695E-01  0.25895303E-01  0.35240179E-01 -0.51670385E-02</point>
      <point r="5.060"> -0.17078643E-01  0.25215819E-01  0.34425707E-01 -0.50087300E-02</point>
      <point r="5.080"> -0.16585590E-01  0.24551646E-01  0.33626246E-01 -0.48548541E-02</point>
      <point r="5.100"> -0.16105250E-01  0.23902514E-01  0.32841626E-01 -0.47052998E-02</point>
      <point r="5.120"> -0.15637347E-01  0.23268158E-01  0.32071677E-01 -0.45599584E-02</point>
      <point r="5.140"> -0.15181606E-01  0.22648313E-01  0.31316225E-01 -0.44187237E-02</point>
      <point r="5.160"> -0.14737757E-01  0.22042717E-01  0.30575098E-01 -0.42814919E-02</point>
      <point r="5.180"> -0.14305537E-01  0.21451113E-01  0.29848122E-01 -0.41481613E-02</point>
      <point r="5.200"> -0.13884684E-01  0.20873243E-01  0.29135123E-01 -0.40186325E-02</point>
      <point r="5.220"> -0.13474943E-01  0.20308855E-01  0.28435926E-01 -0.38928084E-02</point>
      <point r="5.240"> -0.13076062E-01  0.19757697E-01  0.27750354E-01 -0.37705939E-02</point>
      <point r="5.260"> -0.12687794E-01  0.19219523E-01  0.27078234E-01 -0.36518962E-02</point>
      <point r="5.280"> -0.12309897E-01  0.18694088E-01  0.26419388E-01 -0.35366244E-02</point>
      <point r="5.300"> -0.11942132E-01  0.18181149E-01  0.25773641E-01 -0.34246898E-02</point>
      <point r="5.320"> -0.11584265E-01  0.17680466E-01  0.25140816E-01 -0.33160055E-02</point>
      <point r="5.340"> -0.11236067E-01  0.17191804E-01  0.24520738E-01 -0.32104869E-02</point>
      <point r="5.360"> -0.10897311E-01  0.16714930E-01  0.23913229E-01 -0.31080510E-02</point>
      <point r="5.380"> -0.10567775E-01  0.16249612E-01  0.23318115E-01 -0.30086168E-02</point>
      <point r="5.400"> -0.10247244E-01  0.15795623E-01  0.22735218E-01 -0.29121052E-02</point>
      <point r="5.420"> -0.99355021E-02  0.15352738E-01  0.22164365E-01 -0.28184389E-02</point>
      <point r="5.440"> -0.96323414E-02  0.14920736E-01  0.21605380E-01 -0.27275424E-02</point>
      <point r="5.460"> -0.93375561E-02  0.14499397E-01  0.21058087E-01 -0.26393420E-02</point>
      <point r="5.480"> -0.90509449E-02  0.14088506E-01  0.20522314E-01 -0.25537655E-02</point>
      <point r="5.500"> -0.87723101E-02  0.13687851E-01  0.19997886E-01 -0.24707427E-02</point>
      <point r="5.520"> -0.85014581E-02  0.13297220E-01  0.19484630E-01 -0.23902047E-02</point>
      <point r="5.540"> -0.82381990E-02  0.12916407E-01  0.18982375E-01 -0.23120846E-02</point>
      <point r="5.560"> -0.79823466E-02  0.12545208E-01  0.18490949E-01 -0.22363169E-02</point>
      <point r="5.580"> -0.77337184E-02  0.12183422E-01  0.18010181E-01 -0.21628375E-02</point>
      <point r="5.600"> -0.74921357E-02  0.11830851E-01  0.17539902E-01 -0.20915841E-02</point>
      <point r="5.620"> -0.72574232E-02  0.11487300E-01  0.17079943E-01 -0.20224957E-02</point>
      <point r="5.640"> -0.70294094E-02  0.11152577E-01  0.16630135E-01 -0.19555130E-02</point>
      <point r="5.660"> -0.68079260E-02  0.10826491E-01  0.16190314E-01 -0.18905778E-02</point>
      <point r="5.680"> -0.65928084E-02  0.10508857E-01  0.15760311E-01 -0.18276336E-02</point>
      <point r="5.700"> -0.63838955E-02  0.10199492E-01  0.15339965E-01 -0.17666251E-02</point>
      <point r="5.720"> -0.61810293E-02  0.98982149E-02  0.14929109E-01 -0.17074984E-02</point>
      <point r="5.740"> -0.59840552E-02  0.96048478E-02  0.14527584E-01 -0.16502010E-02</point>
      <point r="5.760"> -0.57928220E-02  0.93192161E-02  0.14135227E-01 -0.15946816E-02</point>
      <point r="5.780"> -0.56071816E-02  0.90411479E-02  0.13751880E-01 -0.15408902E-02</point>
      <point r="5.800"> -0.54269892E-02  0.87704742E-02  0.13377384E-01 -0.14887781E-02</point>
      <point r="5.820"> -0.52521030E-02  0.85070287E-02  0.13011581E-01 -0.14382976E-02</point>
      <point r="5.840"> -0.50823843E-02  0.82506478E-02  0.12654318E-01 -0.13894025E-02</point>
      <point r="5.860"> -0.49176974E-02  0.80011711E-02  0.12305439E-01 -0.13420475E-02</point>
      <point r="5.880"> -0.47579098E-02  0.77584404E-02  0.11964791E-01 -0.12961887E-02</point>
      <point r="5.900"> -0.46028916E-02  0.75223008E-02  0.11632225E-01 -0.12517829E-02</point>
      <point r="5.920"> -0.44525161E-02  0.72925997E-02  0.11307590E-01 -0.12087885E-02</point>
      <point r="5.940"> -0.43066593E-02  0.70691874E-02  0.10990738E-01 -0.11671646E-02</point>
      <point r="5.960"> -0.41651999E-02  0.68519170E-02  0.10681523E-01 -0.11268715E-02</point>
      <point r="5.980"> -0.40280196E-02  0.66406439E-02  0.10379800E-01 -0.10878705E-02</point>
      <point r="6.000"> -0.38950025E-02  0.64352265E-02  0.10085425E-01 -0.10501239E-02</point>
      <point r="6.020"> -0.37660357E-02  0.62355255E-02  0.97982567E-02 -0.10135949E-02</point>
      <point r="6.040"> -0.36410086E-02  0.60414045E-02  0.95181553E-02 -0.97824775E-03</point>
      <point r="6.060"> -0.35198132E-02  0.58527293E-02  0.92449821E-02 -0.94404763E-03</point>
      <point r="6.080"> -0.34023443E-02  0.56693685E-02  0.89786005E-02 -0.91096059E-03</point>
      <point r="6.100"> -0.32884988E-02  0.54911931E-02  0.87188754E-02 -0.87895360E-03</point>
      <point r="6.120"> -0.31781764E-02  0.53180767E-02  0.84656734E-02 -0.84799449E-03</point>
      <point r="6.140"> -0.30712789E-02  0.51498951E-02  0.82188627E-02 -0.81805195E-03</point>
      <point r="6.160"> -0.29677104E-02  0.49865266E-02  0.79783134E-02 -0.78909550E-03</point>
      <point r="6.180"> -0.28673777E-02  0.48278522E-02  0.77438973E-02 -0.76109548E-03</point>
      <point r="6.200"> -0.27701894E-02  0.46737548E-02  0.75154878E-02 -0.73402301E-03</point>
      <point r="6.220"> -0.26760566E-02  0.45241199E-02  0.72929601E-02 -0.70785002E-03</point>
      <point r="6.240"> -0.25848924E-02  0.43788352E-02  0.70761911E-02 -0.68254917E-03</point>
      <point r="6.260"> -0.24966121E-02  0.42377909E-02  0.68650593E-02 -0.65809386E-03</point>
      <point r="6.280"> -0.24111332E-02  0.41008791E-02  0.66594452E-02 -0.63445825E-03</point>
      <point r="6.300"> -0.23283750E-02  0.39679944E-02  0.64592309E-02 -0.61161716E-03</point>
      <point r="6.320"> -0.22482590E-02  0.38390334E-02  0.62643002E-02 -0.58954613E-03</point>
      <point r="6.340"> -0.21707087E-02  0.37138949E-02  0.60745387E-02 -0.56822136E-03</point>
      <point r="6.360"> -0.20956493E-02  0.35924800E-02  0.58898336E-02 -0.54761973E-03</point>
      <point r="6.380"> -0.20230082E-02  0.34746916E-02  0.57100739E-02 -0.52771874E-03</point>
      <point r="6.400"> -0.19527142E-02  0.33604349E-02  0.55351505E-02 -0.50849651E-03</point>
      <point r="6.420"> -0.18846984E-02  0.32496171E-02  0.53649557E-02 -0.48993181E-03</point>
      <point r="6.440"> -0.18188935E-02  0.31421472E-02  0.51993836E-02 -0.47200396E-03</point>
      <point r="6.460"> -0.17552337E-02  0.30379365E-02  0.50383302E-02 -0.45469289E-03</point>
      <point r="6.480"> -0.16936552E-02  0.29368980E-02  0.48816929E-02 -0.43797910E-03</point>
      <point r="6.500"> -0.16340958E-02  0.28389467E-02  0.47293711E-02 -0.42184362E-03</point>
      <point r="6.520"> -0.15764948E-02  0.27439995E-02  0.45812655E-02 -0.40626805E-03</point>
      <point r="6.540"> -0.15207933E-02  0.26519751E-02  0.44372789E-02 -0.39123448E-03</point>
      <point r="6.560"> -0.14669338E-02  0.25627941E-02  0.42973154E-02 -0.37672556E-03</point>
      <point r="6.580"> -0.14148603E-02  0.24763788E-02  0.41612808E-02 -0.36272440E-03</point>
      <point r="6.600"> -0.13645184E-02  0.23926533E-02  0.40290828E-02 -0.34921461E-03</point>
      <point r="6.620"> -0.13158552E-02  0.23115435E-02  0.39006304E-02 -0.33618028E-03</point>
      <point r="6.640"> -0.12688191E-02  0.22329769E-02  0.37758343E-02 -0.32360598E-03</point>
      <point r="6.660"> -0.12233599E-02  0.21568828E-02  0.36546070E-02 -0.31147669E-03</point>
      <point r="6.680"> -0.11794289E-02  0.20831921E-02  0.35368623E-02 -0.29977787E-03</point>
      <point r="6.700"> -0.11370372E-02  0.20119554E-02  0.34227488E-02 -0.28851249E-03</point>
      <point r="6.720"> -0.10960393E-02  0.19429068E-02  0.33117913E-02 -0.27763414E-03</point>
      <point r="6.740"> -0.10564305E-02  0.18760627E-02  0.32040656E-02 -0.26714501E-03</point>
      <point r="6.760"> -0.10181670E-02  0.18113601E-02  0.30994916E-02 -0.25703220E-03</point>
      <point r="6.780"> -0.98120652E-03  0.17487377E-02  0.29979906E-02 -0.24728317E-03</point>
      <point r="6.800"> -0.94550762E-03  0.16881355E-02  0.28994855E-02 -0.23788579E-03</point>
      <point r="6.820"> -0.91103022E-03  0.16294950E-02  0.28039009E-02 -0.22882827E-03</point>
      <point r="6.840"> -0.87773533E-03  0.15727593E-02  0.27111625E-02 -0.22009920E-03</point>
      <point r="6.860"> -0.84558509E-03  0.15178726E-02  0.26211978E-02 -0.21168750E-03</point>
      <point r="6.880"> -0.81454269E-03  0.14647808E-02  0.25339355E-02 -0.20358242E-03</point>
      <point r="6.900"> -0.78457239E-03  0.14134310E-02  0.24493060E-02 -0.19577357E-03</point>
      <point r="6.920"> -0.75563949E-03  0.13637715E-02  0.23672408E-02 -0.18825085E-03</point>
      <point r="6.940"> -0.72771029E-03  0.13157522E-02  0.22876730E-02 -0.18100449E-03</point>
      <point r="6.960"> -0.70075207E-03  0.12693241E-02  0.22105371E-02 -0.17402501E-03</point>
      <point r="6.980"> -0.67473306E-03  0.12244394E-02  0.21357690E-02 -0.16730324E-03</point>
      <point r="7.000"> -0.64962244E-03  0.11810516E-02  0.20633058E-02 -0.16083027E-03</point>
      <point r="7.020"> -0.62539029E-03  0.11391154E-02  0.19930860E-02 -0.15459749E-03</point>
      <point r="7.040"> -0.60200757E-03  0.10985868E-02  0.19250495E-02 -0.14859657E-03</point>
      <point r="7.060"> -0.57944614E-03  0.10594227E-02  0.18591373E-02 -0.14281941E-03</point>
      <point r="7.080"> -0.55767867E-03  0.10215813E-02  0.17952920E-02 -0.13725820E-03</point>
      <point r="7.100"> -0.53667865E-03  0.98502193E-03  0.17334572E-02 -0.13190535E-03</point>
      <point r="7.120"> -0.51642041E-03  0.94970493E-03  0.16735779E-02 -0.12675354E-03</point>
      <point r="7.140"> -0.49687902E-03  0.91559171E-03  0.16156001E-02 -0.12179567E-03</point>
      <point r="7.160"> -0.47803033E-03  0.88264475E-03  0.15594714E-02 -0.11702485E-03</point>
      <point r="7.180"> -0.45985092E-03  0.85082753E-03  0.15051402E-02 -0.11243446E-03</point>
      <point r="7.200"> -0.44231812E-03  0.82010452E-03  0.14525562E-02 -0.10801804E-03</point>
      <point r="7.220"> -0.42540992E-03  0.79044116E-03  0.14016705E-02 -0.10376937E-03</point>
      <point r="7.240"> -0.40910502E-03  0.76180387E-03  0.13524350E-02 -0.99682445E-04</point>
      <point r="7.260"> -0.39338279E-03  0.73415996E-03  0.13048029E-02 -0.95751426E-04</point>
      <point r="7.280"> -0.37822322E-03  0.70747768E-03  0.12587284E-02 -0.91970683E-04</point>
      <point r="7.300"> -0.36360696E-03  0.68172616E-03  0.12141669E-02 -0.88334771E-04</point>
      <point r="7.320"> -0.34951525E-03  0.65687539E-03  0.11710748E-02 -0.84838422E-04</point>
      <point r="7.340"> -0.33592994E-03  0.63289621E-03  0.11294095E-02 -0.81476545E-04</point>
      <point r="7.360"> -0.32283345E-03  0.60976030E-03  0.10891295E-02 -0.78244219E-04</point>
      <point r="7.380"> -0.31020877E-03  0.58744014E-03  0.10501942E-02 -0.75136688E-04</point>
      <point r="7.400"> -0.29803943E-03  0.56590899E-03  0.10125643E-02 -0.72149353E-04</point>
      <point r="7.420"> -0.28630951E-03  0.54514090E-03  0.97620096E-03 -0.69277775E-04</point>
      <point r="7.440"> -0.27500357E-03  0.52511066E-03  0.94106674E-03 -0.66517662E-04</point>
      <point r="7.460"> -0.26410672E-03  0.50579379E-03  0.90712493E-03 -0.63864868E-04</point>
      <point r="7.480"> -0.25360452E-03  0.48716653E-03  0.87433978E-03 -0.61315391E-04</point>
      <point r="7.500"> -0.24348301E-03  0.46920584E-03  0.84267643E-03 -0.58865363E-04</point>
      <point r="7.520"> -0.23372871E-03  0.45188932E-03  0.81210092E-03 -0.56511051E-04</point>
      <point r="7.540"> -0.22432856E-03  0.43519527E-03  0.78258013E-03 -0.54248850E-04</point>
      <point r="7.560"> -0.21526996E-03  0.41910261E-03  0.75408181E-03 -0.52075281E-04</point>
      <point r="7.580"> -0.20654072E-03  0.40359091E-03  0.72657453E-03 -0.49986985E-04</point>
      <point r="7.600"> -0.19812904E-03  0.38864035E-03  0.70002768E-03 -0.47980720E-04</point>
      <point r="7.620"> -0.19002355E-03  0.37423171E-03  0.67441142E-03 -0.46053358E-04</point>
      <point r="7.640"> -0.18221323E-03  0.36034635E-03  0.64969674E-03 -0.44201882E-04</point>
      <point r="7.660"> -0.17468746E-03  0.34696620E-03  0.62585534E-03 -0.42423379E-04</point>
      <point r="7.680"> -0.16743598E-03  0.33407374E-03  0.60285971E-03 -0.40715042E-04</point>
      <point r="7.700"> -0.16044885E-03  0.32165200E-03  0.58068304E-03 -0.39074163E-04</point>
      <point r="7.720"> -0.15371652E-03  0.30968453E-03  0.55929927E-03 -0.37498131E-04</point>
      <point r="7.740"> -0.14722972E-03  0.29815539E-03  0.53868300E-03 -0.35984426E-04</point>
      <point r="7.760"> -0.14097953E-03  0.28704915E-03  0.51880956E-03 -0.34530623E-04</point>
      <point r="7.780"> -0.13495733E-03  0.27635083E-03  0.49965491E-03 -0.33134381E-04</point>
      <point r="7.800"> -0.12915480E-03  0.26604597E-03  0.48119570E-03 -0.31793447E-04</point>
      <point r="7.820"> -0.12356393E-03  0.25612053E-03  0.46340919E-03 -0.30505647E-04</point>
      <point r="7.840"> -0.11817695E-03  0.24656093E-03  0.44627329E-03 -0.29268889E-04</point>
      <point r="7.860"> -0.11298640E-03  0.23735404E-03  0.42976653E-03 -0.28081157E-04</point>
      <point r="7.880"> -0.10798508E-03  0.22848713E-03  0.41386801E-03 -0.26940508E-04</point>
      <point r="7.900"> -0.10316602E-03  0.21994789E-03  0.39855743E-03 -0.25845072E-04</point>
      <point r="7.920"> -0.98522523E-04  0.21172441E-03  0.38381509E-03 -0.24793050E-04</point>
      <point r="7.940"> -0.94048129E-04  0.20380516E-03  0.36962180E-03 -0.23782708E-04</point>
      <point r="7.960"> -0.89736604E-04  0.19617901E-03  0.35595896E-03 -0.22812377E-04</point>
      <point r="7.980"> -0.85581940E-04  0.18883518E-03  0.34280847E-03 -0.21880452E-04</point>
      <point r="8.000"> -0.81578349E-04  0.18176325E-03  0.33015278E-03 -0.20985388E-04</point>
      <point r="8.020"> -0.77720251E-04  0.17495313E-03  0.31797484E-03 -0.20125698E-04</point>
      <point r="8.040"> -0.74002271E-04  0.16839511E-03  0.30625808E-03 -0.19299953E-04</point>
      <point r="8.060"> -0.70419231E-04  0.16207977E-03  0.29498646E-03 -0.18506777E-04</point>
      <point r="8.080"> -0.66966142E-04  0.15599802E-03  0.28414436E-03 -0.17744848E-04</point>
      <point r="8.100"> -0.63638203E-04  0.15014108E-03  0.27371666E-03 -0.17012895E-04</point>
      <point r="8.120"> -0.60430786E-04  0.14450048E-03  0.26368869E-03 -0.16309695E-04</point>
      <point r="8.140"> -0.57339439E-04  0.13906801E-03  0.25404620E-03 -0.15634075E-04</point>
      <point r="8.160"> -0.54359876E-04  0.13383578E-03  0.24477540E-03 -0.14984904E-04</point>
      <point r="8.180"> -0.51487969E-04  0.12879616E-03  0.23586289E-03 -0.14361099E-04</point>
      <point r="8.200"> -0.48719748E-04  0.12394177E-03  0.22729570E-03 -0.13761618E-04</point>
      <point r="8.220"> -0.46051392E-04  0.11926552E-03  0.21906125E-03 -0.13185460E-04</point>
      <point r="8.240"> -0.43479226E-04  0.11476055E-03  0.21114736E-03 -0.12631664E-04</point>
      <point r="8.260"> -0.40999713E-04  0.11042024E-03  0.20354222E-03 -0.12099308E-04</point>
      <point r="8.280"> -0.38609452E-04  0.10623823E-03  0.19623440E-03 -0.11587506E-04</point>
      <point r="8.300"> -0.36305174E-04  0.10220837E-03  0.18921281E-03 -0.11095407E-04</point>
      <point r="8.320"> -0.34083733E-04  0.98324724E-04  0.18246674E-03 -0.10622197E-04</point>
      <point r="8.340"> -0.31942106E-04  0.94581595E-04  0.17598581E-03 -0.10167093E-04</point>
      <point r="8.360"> -0.29877387E-04  0.90973479E-04  0.16975997E-03 -0.97293426E-05</point>
      <point r="8.380"> -0.27886783E-04  0.87495079E-04  0.16377951E-03 -0.93082275E-05</point>
      <point r="8.400"> -0.25967608E-04  0.84141292E-04  0.15803501E-03 -0.89030565E-05</point>
      <point r="8.420"> -0.24117285E-04  0.80907206E-04  0.15251739E-03 -0.85131678E-05</point>
      <point r="8.440"> -0.22333334E-04  0.77788092E-04  0.14721785E-03 -0.81379266E-05</point>
      <point r="8.460"> -0.20613375E-04  0.74779399E-04  0.14212788E-03 -0.77767249E-05</point>
      <point r="8.480"> -0.18955122E-04  0.71876750E-04  0.13723928E-03 -0.74289799E-05</point>
      <point r="8.500"> -0.17356377E-04  0.69075932E-04  0.13254410E-03 -0.70941332E-05</point>
      <point r="8.520"> -0.15815032E-04  0.66372895E-04  0.12803467E-03 -0.67716499E-05</point>
      <point r="8.540"> -0.14329061E-04  0.63763744E-04  0.12370358E-03 -0.64610179E-05</point>
      <point r="8.560"> -0.12896520E-04  0.61244737E-04  0.11954368E-03 -0.61617468E-05</point>
      <point r="8.580"> -0.11515542E-04  0.58812274E-04  0.11554805E-03 -0.58733669E-05</point>
      <point r="8.600"> -0.10184334E-04  0.56462901E-04  0.11171005E-03 -0.55954289E-05</point>
      <point r="8.620"> -0.89011758E-05  0.54193296E-04  0.10802323E-03 -0.53275026E-05</point>
      <point r="8.640"> -0.76644169E-05  0.52000272E-04  0.10448138E-03 -0.50691766E-05</point>
      <point r="8.660"> -0.64724719E-05  0.49880768E-04  0.10107854E-03 -0.48200570E-05</point>
      <point r="8.680"> -0.53238200E-05  0.47831844E-04  0.97808914E-04 -0.45797674E-05</point>
      <point r="8.700"> -0.42170014E-05  0.45850683E-04  0.94666955E-04 -0.43479475E-05</point>
      <point r="8.720"> -0.31506156E-05  0.43934580E-04  0.91647302E-04 -0.41242530E-05</point>
      <point r="8.740"> -0.21233183E-05  0.42080940E-04  0.88744791E-04 -0.39083548E-05</point>
      <point r="8.760"> -0.11338198E-05  0.40287275E-04  0.85954448E-04 -0.36999382E-05</point>
      <point r="8.780"> -0.18088249E-06  0.38551203E-04  0.83271486E-04 -0.34987026E-05</point>
      <point r="8.800">  0.73668127E-06  0.36870437E-04  0.80691296E-04 -0.33043604E-05</point>
      <point r="8.820">  0.16200110E-05  0.35242788E-04  0.78209440E-04 -0.31166374E-05</point>
      <point r="8.840">  0.24702003E-05  0.33666159E-04  0.75821653E-04 -0.29352713E-05</point>
      <point r="8.860">  0.32882985E-05  0.32138542E-04  0.73523830E-04 -0.27600118E-05</point>
      <point r="8.880">  0.40753127E-05  0.30658014E-04  0.71312024E-04 -0.25906196E-05</point>
      <point r="8.900">  0.48322094E-05  0.29222736E-04  0.69182441E-04 -0.24268665E-05</point>
      <point r="8.920">  0.55599161E-05  0.27830946E-04  0.67131437E-04 -0.22685346E-05</point>
      <point r="8.940">  0.62593233E-05  0.26480960E-04  0.65155511E-04 -0.21154159E-05</point>
      <point r="8.960">  0.69312858E-05  0.25171168E-04  0.63251298E-04 -0.19673118E-05</point>
      <point r="8.980">  0.75766242E-05  0.23900029E-04  0.61415572E-04 -0.18240330E-05</point>
      <point r="9.000">  0.81961265E-05  0.22666072E-04  0.59645233E-04 -0.16853987E-05</point>
      <point r="9.020">  0.87905496E-05  0.21467891E-04  0.57937310E-04 -0.15512365E-05</point>
      <point r="9.040">  0.93606201E-05  0.20304141E-04  0.56288953E-04 -0.14213818E-05</point>
      <point r="9.060">  0.99070366E-05  0.19173540E-04  0.54697427E-04 -0.12956780E-05</point>
      <point r="9.080">  0.10430470E-04  0.18074862E-04  0.53160114E-04 -0.11739752E-05</point>
      <point r="9.100">  0.10931565E-04  0.17006937E-04  0.51674504E-04 -0.10561309E-05</point>
      <point r="9.120">  0.11410942E-04  0.15968650E-04  0.50238195E-04 -0.94200905E-06</point>
      <point r="9.140">  0.11869197E-04  0.14958936E-04  0.48848885E-04 -0.83147983E-06</point>
      <point r="9.160">  0.12306903E-04  0.13976780E-04  0.47504373E-04 -0.72441959E-06</point>
      <point r="9.180">  0.12724612E-04  0.13021212E-04  0.46202552E-04 -0.62071039E-06</point>
      <point r="9.200">  0.13122855E-04  0.12091310E-04  0.44941409E-04 -0.52023979E-06</point>
      <point r="9.220">  0.13502143E-04  0.11186194E-04  0.43719018E-04 -0.42290059E-06</point>
      <point r="9.240">  0.13862969E-04  0.10305025E-04  0.42533541E-04 -0.32859059E-06</point>
      <point r="9.260">  0.14205807E-04  0.94470059E-05  0.41383221E-04 -0.23721234E-06</point>
      <point r="9.280">  0.14531114E-04  0.86113751E-05  0.40266381E-04 -0.14867291E-06</point>
      <point r="9.300">  0.14839331E-04  0.77974086E-05  0.39181421E-04 -0.62883708E-07</point>
      <point r="9.320">  0.15130884E-04  0.70044172E-05  0.38126816E-04  0.20239759E-07</point>
      <point r="9.340">  0.15406184E-04  0.62317447E-05  0.37101111E-04  0.10077807E-06</point>
      <point r="9.360">  0.15665626E-04  0.54787669E-05  0.36102919E-04  0.17880810E-06</point>
      <point r="9.380">  0.15909594E-04  0.47448898E-05  0.35130922E-04  0.25440317E-06</point>
      <point r="9.400">  0.16138458E-04  0.40295486E-05  0.34183863E-04  0.32763325E-06</point>
      <point r="9.420">  0.16352577E-04  0.33322059E-05  0.33260546E-04  0.39856511E-06</point>
      <point r="9.440">  0.16552297E-04  0.26523509E-05  0.32359836E-04  0.46726248E-06</point>
      <point r="9.460">  0.16737954E-04  0.19894979E-05  0.31480652E-04  0.53378623E-06</point>
      <point r="9.480">  0.16909873E-04  0.13431849E-05  0.30621969E-04  0.59819449E-06</point>
      <point r="9.500">  0.17068370E-04  0.71297325E-06  0.29782814E-04  0.66054278E-06</point>
      <point r="9.520">  0.17213751E-04  0.98445550E-07  0.28962263E-04  0.72088419E-06</point>
      <point r="9.540">  0.17346312E-04 -0.50079470E-06  0.28159442E-04  0.77926946E-06</point>
      <point r="9.560">  0.17466343E-04 -0.10851244E-05  0.27373520E-04  0.83574713E-06</point>
      <point r="9.580">  0.17574124E-04 -0.16549017E-05  0.26603713E-04  0.89036365E-06</point>
      <point r="9.600">  0.17669928E-04 -0.22104672E-05  0.25849279E-04  0.94316348E-06</point>
      <point r="9.620">  0.17754020E-04 -0.27521445E-05  0.25109516E-04  0.99418924E-06</point>
      <point r="9.640">  0.17826660E-04 -0.32802413E-05  0.24383760E-04  0.10434817E-05</point>
      <point r="9.660">  0.17888100E-04 -0.37950502E-05  0.23671386E-04  0.10910801E-05</point>
      <point r="9.680">  0.17938586E-04 -0.42968493E-05  0.22971805E-04  0.11370220E-05</point>
      <point r="9.700">  0.17978358E-04 -0.47859034E-05  0.22284459E-04  0.11813435E-05</point>
      <point r="9.720">  0.18007651E-04 -0.52624644E-05  0.21608827E-04  0.12240791E-05</point>
      <point r="9.740">  0.18026695E-04 -0.57267719E-05  0.20944415E-04  0.12652624E-05</point>
      <point r="9.760">  0.18035714E-04 -0.61790543E-05  0.20290761E-04  0.13049253E-05</point>
      <point r="9.780">  0.18034929E-04 -0.66195291E-05  0.19647433E-04  0.13430989E-05</point>
      <point r="9.800">  0.18024554E-04 -0.70484036E-05  0.19014021E-04  0.13798128E-05</point>
      <point r="9.820">  0.18004800E-04 -0.74658756E-05  0.18390147E-04  0.14150959E-05</point>
      <point r="9.840">  0.17975876E-04 -0.78721338E-05  0.17775453E-04  0.14489760E-05</point>
      <point r="9.860">  0.17937983E-04 -0.82673586E-05  0.17169607E-04  0.14814798E-05</point>
      <point r="9.880">  0.17891323E-04 -0.86517222E-05  0.16572298E-04  0.15126334E-05</point>
      <point r="9.900">  0.17836090E-04 -0.90253895E-05  0.15983237E-04  0.15424619E-05</point>
      <point r="9.920">  0.17772479E-04 -0.93885186E-05  0.15402155E-04  0.15709898E-05</point>
      <point r="9.940">  0.17700679E-04 -0.97412608E-05  0.14828802E-04  0.15982407E-05</point>
      <point r="9.960">  0.17620876E-04 -0.10083761E-04  0.14262947E-04  0.16242376E-05</point>
      <point r="9.980">  0.17533256E-04 -0.10416160E-04  0.13704375E-04  0.16490029E-05</point>
      <point r="10.000">  0.17438001E-04 -0.10738591E-04  0.13152888E-04  0.16725584E-05</point>
      <point r="10.020">  0.17335287E-04 -0.11051184E-04  0.12608304E-04  0.16949253E-05</point>
      <point r="10.040">  0.17225294E-04 -0.11354063E-04  0.12070456E-04  0.17161243E-05</point>
      <point r="10.060">  0.17108194E-04 -0.11647350E-04  0.11539189E-04  0.17361757E-05</point>
      <point r="10.080">  0.16984161E-04 -0.11931161E-04  0.11014362E-04  0.17550992E-05</point>
      <point r="10.100">  0.16853363E-04 -0.12205609E-04  0.10495847E-04  0.17729143E-05</point>
      <point r="10.120">  0.16715970E-04 -0.12470804E-04  0.99835279E-05  0.17896397E-05</point>
      <point r="10.140">  0.16572146E-04 -0.12726853E-04  0.94772981E-05  0.18052942E-05</point>
      <point r="10.160">  0.16422057E-04 -0.12973860E-04  0.89770622E-05  0.18198960E-05</point>
      <point r="10.180">  0.16265865E-04 -0.13211925E-04  0.84827343E-05  0.18334630E-05</point>
      <point r="10.200">  0.16103730E-04 -0.13441149E-04  0.79942378E-05  0.18460128E-05</point>
      <point r="10.220">  0.15935811E-04 -0.13661628E-04  0.75115042E-05  0.18575628E-05</point>
      <point r="10.240">  0.15762266E-04 -0.13873458E-04  0.70344732E-05  0.18681300E-05</point>
      <point r="10.260">  0.15583251E-04 -0.14076731E-04  0.65630919E-05  0.18777314E-05</point>
      <point r="10.280">  0.15398921E-04 -0.14271540E-04  0.60973141E-05  0.18863835E-05</point>
      <point r="10.300">  0.15209427E-04 -0.14457976E-04  0.56371002E-05  0.18941028E-05</point>
      <point r="10.320">  0.15014921E-04 -0.14636127E-04  0.51824165E-05  0.19009053E-05</point>
      <point r="10.340">  0.14815553E-04 -0.14806082E-04  0.47332347E-05  0.19068073E-05</point>
      <point r="10.360">  0.14611473E-04 -0.14967928E-04  0.42895317E-05  0.19118244E-05</point>
      <point r="10.380">  0.14402826E-04 -0.15121752E-04  0.38512889E-05  0.19159724E-05</point>
      <point r="10.400">  0.00000000E+00</point>
    </H_spline>
    <S_spline>
      <point r="0.000">  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.35000000E+01</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400">  0.92141382E+00 -0.12355975E+00  0.83908379E+00  0.94453015E+00</point>
      <point r="0.420">  0.91559199E+00 -0.13124306E+00  0.82402247E+00  0.93915227E+00</point>
      <point r="0.440">  0.90975515E+00 -0.13904607E+00  0.80847704E+00  0.93356757E+00</point>
      <point r="0.460">  0.90391471E+00 -0.14695910E+00  0.79247519E+00  0.92778259E+00</point>
      <point r="0.480">  0.89808044E+00 -0.15497188E+00  0.77604463E+00  0.92180392E+00</point>
      <point r="0.500">  0.89226062E+00 -0.16307370E+00  0.75921297E+00  0.91563821E+00</point>
      <point r="0.520">  0.88646211E+00 -0.17125345E+00  0.74200768E+00  0.90929213E+00</point>
      <point r="0.540">  0.88069048E+00 -0.17949976E+00  0.72445603E+00  0.90277241E+00</point>
      <point r="0.560">  0.87495012E+00 -0.18780102E+00  0.70658503E+00  0.89608575E+00</point>
      <point r="0.580">  0.86924438E+00 -0.19614551E+00  0.68842139E+00  0.88923888E+00</point>
      <point r="0.600">  0.86357563E+00 -0.20452142E+00  0.66999147E+00  0.88223850E+00</point>
      <point r="0.620">  0.85794538E+00 -0.21291691E+00  0.65132124E+00  0.87509130E+00</point>
      <point r="0.640">  0.85235438E+00 -0.22132020E+00  0.63243624E+00  0.86780394E+00</point>
      <point r="0.660">  0.84680269E+00 -0.22971959E+00  0.61336156E+00  0.86038302E+00</point>
      <point r="0.680">  0.84128977E+00 -0.23810351E+00  0.59412179E+00  0.85283512E+00</point>
      <point r="0.700">  0.83581456E+00 -0.24646055E+00  0.57474101E+00  0.84516675E+00</point>
      <point r="0.720">  0.83037554E+00 -0.25477953E+00  0.55524277E+00  0.83738436E+00</point>
      <point r="0.740">  0.82497083E+00 -0.26304947E+00  0.53565005E+00  0.82949433E+00</point>
      <point r="0.760">  0.81959819E+00 -0.27125969E+00  0.51598526E+00  0.82150298E+00</point>
      <point r="0.780">  0.81425514E+00 -0.27939979E+00  0.49627021E+00  0.81341652E+00</point>
      <point r="0.800">  0.80893896E+00 -0.28745966E+00  0.47652611E+00  0.80524110E+00</point>
      <point r="0.820">  0.80364679E+00 -0.29542955E+00  0.45677356E+00  0.79698277E+00</point>
      <point r="0.840">  0.79837560E+00 -0.30330003E+00  0.43703250E+00  0.78864749E+00</point>
      <point r="0.860">  0.79312231E+00 -0.31106205E+00  0.41732227E+00  0.78024113E+00</point>
      <point r="0.880">  0.78788377E+00 -0.31870691E+00  0.39766155E+00  0.77176945E+00</point>
      <point r="0.900">  0.78265679E+00 -0.32622630E+00  0.37806838E+00  0.76323810E+00</point>
      <point r="0.920">  0.77743822E+00 -0.33361229E+00  0.35856016E+00  0.75465263E+00</point>
      <point r="0.940">  0.77222490E+00 -0.34085734E+00  0.33915363E+00  0.74601848E+00</point>
      <point r="0.960">  0.76701376E+00 -0.34795432E+00  0.31986487E+00  0.73734099E+00</point>
      <point r="0.980">  0.76180178E+00 -0.35489649E+00  0.30070935E+00  0.72862536E+00</point>
      <point r="1.000">  0.75658604E+00 -0.36167749E+00  0.28170185E+00  0.71987671E+00</point>
      <point r="1.020">  0.75136370E+00 -0.36829138E+00  0.26285653E+00  0.71110001E+00</point>
      <point r="1.040">  0.74613207E+00 -0.37473260E+00  0.24418692E+00  0.70230014E+00</point>
      <point r="1.060">  0.74088856E+00 -0.38099600E+00  0.22570590E+00  0.69348184E+00</point>
      <point r="1.080">  0.73563075E+00 -0.38707681E+00  0.20742574E+00  0.68464974E+00</point>
      <point r="1.100">  0.73035632E+00 -0.39297065E+00  0.18935807E+00  0.67580836E+00</point>
      <point r="1.120">  0.72506313E+00 -0.39867350E+00  0.17151393E+00  0.66696208E+00</point>
      <point r="1.140">  0.71974920E+00 -0.40418175E+00  0.15390376E+00  0.65811518E+00</point>
      <point r="1.160">  0.71441268E+00 -0.40949212E+00  0.13653738E+00  0.64927182E+00</point>
      <point r="1.180">  0.70905191E+00 -0.41460172E+00  0.11942404E+00  0.64043601E+00</point>
      <point r="1.200">  0.70366537E+00 -0.41950799E+00  0.10257244E+00  0.63161168E+00</point>
      <point r="1.220">  0.69825171E+00 -0.42420873E+00  0.85990673E-01  0.62280261E+00</point>
      <point r="1.240">  0.69280976E+00 -0.42870207E+00  0.69686313E-01  0.61401249E+00</point>
      <point r="1.260">  0.68733848E+00 -0.43298647E+00  0.53666379E-01  0.60524487E+00</point>
      <point r="1.280">  0.68183700E+00 -0.43706070E+00  0.37937364E-01  0.59650319E+00</point>
      <point r="1.300">  0.67630463E+00 -0.44092385E+00  0.22505242E-01  0.58779078E+00</point>
      <point r="1.320">  0.67074079E+00 -0.44457529E+00  0.73754832E-02  0.57911084E+00</point>
      <point r="1.340">  0.66514509E+00 -0.44801470E+00 -0.74469366E-02  0.57046648E+00</point>
      <point r="1.360">  0.65951726E+00 -0.45124202E+00 -0.21957522E-01  0.56186067E+00</point>
      <point r="1.380">  0.65385719E+00 -0.45425747E+00 -0.36152245E-01  0.55329629E+00</point>
      <point r="1.400">  0.64816490E+00 -0.45706153E+00 -0.50027533E-01  0.54477611E+00</point>
      <point r="1.420">  0.64244054E+00 -0.45965493E+00 -0.63580257E-01  0.53630277E+00</point>
      <point r="1.440">  0.63668440E+00 -0.46203862E+00 -0.76807718E-01  0.52787883E+00</point>
      <point r="1.460">  0.63089689E+00 -0.46421380E+00 -0.89707630E-01  0.51950672E+00</point>
      <point r="1.480">  0.62507851E+00 -0.46618187E+00 -0.10227812E+00  0.51118878E+00</point>
      <point r="1.500">  0.61922993E+00 -0.46794446E+00 -0.11451769E+00  0.50292726E+00</point>
      <point r="1.520">  0.61335186E+00 -0.46950338E+00 -0.12642524E+00  0.49472428E+00</point>
      <point r="1.540">  0.60744517E+00 -0.47086063E+00 -0.13800001E+00  0.48658188E+00</point>
      <point r="1.560">  0.60151080E+00 -0.47201839E+00 -0.14924163E+00  0.47850200E+00</point>
      <point r="1.580">  0.59554978E+00 -0.47297902E+00 -0.16015004E+00  0.47048649E+00</point>
      <point r="1.600">  0.58956322E+00 -0.47374503E+00 -0.17072551E+00  0.46253709E+00</point>
      <point r="1.620">  0.58355233E+00 -0.47431907E+00 -0.18096864E+00  0.45465546E+00</point>
      <point r="1.640">  0.57751839E+00 -0.47470395E+00 -0.19088032E+00  0.44684318E+00</point>
      <point r="1.660">  0.57146272E+00 -0.47490261E+00 -0.20046174E+00  0.43910171E+00</point>
      <point r="1.680">  0.56538675E+00 -0.47491809E+00 -0.20971436E+00  0.43143246E+00</point>
      <point r="1.700">  0.55929193E+00 -0.47475357E+00 -0.21863992E+00  0.42383673E+00</point>
      <point r="1.720">  0.55317978E+00 -0.47441233E+00 -0.22724040E+00  0.41631576E+00</point>
      <point r="1.740">  0.54705188E+00 -0.47389774E+00 -0.23551804E+00  0.40887068E+00</point>
      <point r="1.760">  0.54090984E+00 -0.47321328E+00 -0.24347529E+00  0.40150257E+00</point>
      <point r="1.780">  0.53475530E+00 -0.47236248E+00 -0.25111485E+00  0.39421241E+00</point>
      <point r="1.800">  0.52858995E+00 -0.47134897E+00 -0.25843961E+00  0.38700112E+00</point>
      <point r="1.820">  0.52241550E+00 -0.47017644E+00 -0.26545267E+00  0.37986955E+00</point>
      <point r="1.840">  0.51623371E+00 -0.46884864E+00 -0.27215732E+00  0.37281846E+00</point>
      <point r="1.860">  0.51004633E+00 -0.46736939E+00 -0.27855701E+00  0.36584856E+00</point>
      <point r="1.880">  0.50385513E+00 -0.46574251E+00 -0.28465538E+00  0.35896048E+00</point>
      <point r="1.900">  0.49766192E+00 -0.46397192E+00 -0.29045624E+00  0.35215478E+00</point>
      <point r="1.920">  0.49146850E+00 -0.46206153E+00 -0.29596351E+00  0.34543198E+00</point>
      <point r="1.940">  0.48527667E+00 -0.46001530E+00 -0.30118128E+00  0.33879252E+00</point>
      <point r="1.960">  0.47908824E+00 -0.45783720E+00 -0.30611376E+00  0.33223676E+00</point>
      <point r="1.980">  0.47290504E+00 -0.45553123E+00 -0.31076529E+00  0.32576505E+00</point>
      <point r="2.000">  0.46672885E+00 -0.45310139E+00 -0.31514032E+00  0.31937764E+00</point>
      <point r="2.020">  0.46056149E+00 -0.45055168E+00 -0.31924340E+00  0.31307475E+00</point>
      <point r="2.040">  0.45440474E+00 -0.44788613E+00 -0.32307916E+00  0.30685654E+00</point>
      <point r="2.060">  0.44826038E+00 -0.44510874E+00 -0.32665236E+00  0.30072310E+00</point>
      <point r="2.080">  0.44213017E+00 -0.44222351E+00 -0.32996779E+00  0.29467450E+00</point>
      <point r="2.100">  0.43601585E+00 -0.43923443E+00 -0.33303036E+00  0.28871074E+00</point>
      <point r="2.120">  0.42991914E+00 -0.43614546E+00 -0.33584500E+00  0.28283179E+00</point>
      <point r="2.140">  0.42384175E+00 -0.43296057E+00 -0.33841673E+00  0.27703756E+00</point>
      <point r="2.160">  0.41778535E+00 -0.42968367E+00 -0.34075060E+00  0.27132792E+00</point>
      <point r="2.180">  0.41175160E+00 -0.42631866E+00 -0.34285172E+00  0.26570271E+00</point>
      <point r="2.200">  0.40574210E+00 -0.42286942E+00 -0.34472522E+00  0.26016171E+00</point>
      <point r="2.220">  0.39975846E+00 -0.41933978E+00 -0.34637627E+00  0.25470468E+00</point>
      <point r="2.240">  0.39380223E+00 -0.41573353E+00 -0.34781006E+00  0.24933132E+00</point>
      <point r="2.260">  0.38787495E+00 -0.41205443E+00 -0.34903181E+00  0.24404132E+00</point>
      <point r="2.280">  0.38197809E+00 -0.40830619E+00 -0.35004673E+00  0.23883433E+00</point>
      <point r="2.300">  0.37611313E+00 -0.40449249E+00 -0.35086005E+00  0.23370994E+00</point>
      <point r="2.320">  0.37028147E+00 -0.40061694E+00 -0.35147701E+00  0.22866775E+00</point>
      <point r="2.340">  0.36448451E+00 -0.39668312E+00 -0.35190283E+00  0.22370729E+00</point>
      <point r="2.360">  0.35872359E+00 -0.39269454E+00 -0.35214273E+00  0.21882809E+00</point>
      <point r="2.380">  0.35300001E+00 -0.38865466E+00 -0.35220190E+00  0.21402965E+00</point>
      <point r="2.400">  0.34731503E+00 -0.38456691E+00 -0.35208555E+00  0.20931141E+00</point>
      <point r="2.420">  0.34166989E+00 -0.38043463E+00 -0.35179884E+00  0.20467284E+00</point>
      <point r="2.440">  0.33606577E+00 -0.37626112E+00 -0.35134689E+00  0.20011334E+00</point>
      <point r="2.460">  0.33050379E+00 -0.37204961E+00 -0.35073482E+00  0.19563230E+00</point>
      <point r="2.480">  0.32498508E+00 -0.36780327E+00 -0.34996770E+00  0.19122911E+00</point>
      <point r="2.500">  0.31951068E+00 -0.36352521E+00 -0.34905057E+00  0.18690310E+00</point>
      <point r="2.520">  0.31408161E+00 -0.35921849E+00 -0.34798841E+00  0.18265362E+00</point>
      <point r="2.540">  0.30869884E+00 -0.35488608E+00 -0.34678618E+00  0.17847998E+00</point>
      <point r="2.560">  0.30336330E+00 -0.35053090E+00 -0.34544878E+00  0.17438148E+00</point>
      <point r="2.580">  0.29807589E+00 -0.34615581E+00 -0.34398105E+00  0.17035739E+00</point>
      <point r="2.600">  0.29283743E+00 -0.34176360E+00 -0.34238780E+00  0.16640699E+00</point>
      <point r="2.620">  0.28764874E+00 -0.33735697E+00 -0.34067376E+00  0.16252953E+00</point>
      <point r="2.640">  0.28251058E+00 -0.33293859E+00 -0.33884362E+00  0.15872424E+00</point>
      <point r="2.660">  0.27742367E+00 -0.32851104E+00 -0.33690200E+00  0.15499036E+00</point>
      <point r="2.680">  0.27238867E+00 -0.32407684E+00 -0.33485345E+00  0.15132710E+00</point>
      <point r="2.700">  0.26740624E+00 -0.31963843E+00 -0.33270247E+00  0.14773367E+00</point>
      <point r="2.720">  0.26247695E+00 -0.31519821E+00 -0.33045349E+00  0.14420926E+00</point>
      <point r="2.740">  0.25760137E+00 -0.31075849E+00 -0.32811087E+00  0.14075306E+00</point>
      <point r="2.760">  0.25278001E+00 -0.30632151E+00 -0.32567888E+00  0.13736426E+00</point>
      <point r="2.780">  0.24801335E+00 -0.30188945E+00 -0.32316176E+00  0.13404202E+00</point>
      <point r="2.800">  0.24330181E+00 -0.29746444E+00 -0.32056364E+00  0.13078551E+00</point>
      <point r="2.820">  0.23864581E+00 -0.29304851E+00 -0.31788860E+00  0.12759389E+00</point>
      <point r="2.840">  0.23404569E+00 -0.28864365E+00 -0.31514062E+00  0.12446632E+00</point>
      <point r="2.860">  0.22950177E+00 -0.28425178E+00 -0.31232364E+00  0.12140196E+00</point>
      <point r="2.880">  0.22501435E+00 -0.27987473E+00 -0.30944148E+00  0.11839995E+00</point>
      <point r="2.900">  0.22058368E+00 -0.27551431E+00 -0.30649792E+00  0.11545943E+00</point>
      <point r="2.920">  0.21620996E+00 -0.27117222E+00 -0.30349664E+00  0.11257956E+00</point>
      <point r="2.940">  0.21189338E+00 -0.26685013E+00 -0.30044125E+00  0.10975947E+00</point>
      <point r="2.960">  0.20763408E+00 -0.26254963E+00 -0.29733528E+00  0.10699831E+00</point>
      <point r="2.980">  0.20343219E+00 -0.25827225E+00 -0.29418217E+00  0.10429522E+00</point>
      <point r="3.000">  0.19928777E+00 -0.25401946E+00 -0.29098531E+00  0.10164934E+00</point>
      <point r="3.020">  0.19520088E+00 -0.24979268E+00 -0.28774797E+00  0.99059805E-01</point>
      <point r="3.040">  0.19117154E+00 -0.24559325E+00 -0.28447337E+00  0.96525769E-01</point>
      <point r="3.060">  0.18719973E+00 -0.24142247E+00 -0.28116464E+00  0.94046374E-01</point>
      <point r="3.080">  0.18328543E+00 -0.23728157E+00 -0.27782484E+00  0.91620766E-01</point>
      <point r="3.100">  0.17942855E+00 -0.23317172E+00 -0.27445692E+00  0.89248096E-01</point>
      <point r="3.120">  0.17562901E+00 -0.22909405E+00 -0.27106379E+00  0.86927516E-01</point>
      <point r="3.140">  0.17188668E+00 -0.22504962E+00 -0.26764825E+00  0.84658183E-01</point>
      <point r="3.160">  0.16820142E+00 -0.22103945E+00 -0.26421305E+00  0.82439257E-01</point>
      <point r="3.180">  0.16457304E+00 -0.21706449E+00 -0.26076084E+00  0.80269902E-01</point>
      <point r="3.200">  0.16100137E+00 -0.21312564E+00 -0.25729420E+00  0.78149286E-01</point>
      <point r="3.220">  0.15748617E+00 -0.20922376E+00 -0.25381563E+00  0.76076583E-01</point>
      <point r="3.240">  0.15402720E+00 -0.20535966E+00 -0.25032755E+00  0.74050973E-01</point>
      <point r="3.260">  0.15062420E+00 -0.20153409E+00 -0.24683233E+00  0.72071639E-01</point>
      <point r="3.280">  0.14727689E+00 -0.19774776E+00 -0.24333224E+00  0.70137772E-01</point>
      <point r="3.300">  0.14398495E+00 -0.19400132E+00 -0.23982947E+00  0.68248567E-01</point>
      <point r="3.320">  0.14074807E+00 -0.19029540E+00 -0.23632615E+00  0.66403228E-01</point>
      <point r="3.340">  0.13756589E+00 -0.18663055E+00 -0.23282435E+00  0.64600963E-01</point>
      <point r="3.360">  0.13443807E+00 -0.18300732E+00 -0.22932604E+00  0.62840988E-01</point>
      <point r="3.380">  0.13136421E+00 -0.17942617E+00 -0.22583315E+00  0.61122526E-01</point>
      <point r="3.400">  0.12834393E+00 -0.17588756E+00 -0.22234750E+00  0.59444808E-01</point>
      <point r="3.420">  0.12537682E+00 -0.17239188E+00 -0.21887089E+00  0.57807070E-01</point>
      <point r="3.440">  0.12246245E+00 -0.16893949E+00 -0.21540502E+00  0.56208558E-01</point>
      <point r="3.460">  0.11960038E+00 -0.16553072E+00 -0.21195153E+00  0.54648525E-01</point>
      <point r="3.480">  0.11679017E+00 -0.16216585E+00 -0.20851199E+00  0.53126231E-01</point>
      <point r="3.500">  0.11403134E+00 -0.15884513E+00 -0.20508792E+00  0.51640946E-01</point>
      <point r="3.520">  0.11132343E+00 -0.15556878E+00 -0.20168077E+00  0.50191947E-01</point>
      <point r="3.540">  0.10866593E+00 -0.15233697E+00 -0.19829192E+00  0.48778520E-01</point>
      <point r="3.560">  0.10605836E+00 -0.14914986E+00 -0.19492270E+00  0.47399957E-01</point>
      <point r="3.580">  0.10350021E+00 -0.14600756E+00 -0.19157436E+00  0.46055562E-01</point>
      <point r="3.600">  0.10099095E+00 -0.14291015E+00 -0.18824813E+00  0.44744646E-01</point>
      <point r="3.620">  0.98530071E-01 -0.13985769E+00 -0.18494515E+00  0.43466528E-01</point>
      <point r="3.640">  0.96117027E-01 -0.13685019E+00 -0.18166650E+00  0.42220538E-01</point>
      <point r="3.660">  0.93751280E-01 -0.13388767E+00 -0.17841323E+00  0.41006012E-01</point>
      <point r="3.680">  0.91432280E-01 -0.13097009E+00 -0.17518633E+00  0.39822297E-01</point>
      <point r="3.700">  0.89159474E-01 -0.12809738E+00 -0.17198671E+00  0.38668748E-01</point>
      <point r="3.720">  0.86932300E-01 -0.12526948E+00 -0.16881527E+00  0.37544730E-01</point>
      <point r="3.740">  0.84750193E-01 -0.12248628E+00 -0.16567283E+00  0.36449616E-01</point>
      <point r="3.760">  0.82612582E-01 -0.11974764E+00 -0.16256016E+00  0.35382788E-01</point>
      <point r="3.780">  0.80518893E-01 -0.11705341E+00 -0.15947801E+00  0.34343638E-01</point>
      <point r="3.800">  0.78468547E-01 -0.11440343E+00 -0.15642705E+00  0.33331567E-01</point>
      <point r="3.820">  0.76460962E-01 -0.11179750E+00 -0.15340793E+00  0.32345984E-01</point>
      <point r="3.840">  0.74495556E-01 -0.10923541E+00 -0.15042124E+00  0.31386308E-01</point>
      <point r="3.860">  0.72571741E-01 -0.10671693E+00 -0.14746754E+00  0.30451966E-01</point>
      <point r="3.880">  0.70688928E-01 -0.10424180E+00 -0.14454732E+00  0.29542396E-01</point>
      <point r="3.900">  0.68846529E-01 -0.10180977E+00 -0.14166107E+00  0.28657044E-01</point>
      <point r="3.920">  0.67043954E-01 -0.99420537E-01 -0.13880921E+00  0.27795363E-01</point>
      <point r="3.940">  0.65280610E-01 -0.97073817E-01 -0.13599213E+00  0.26956820E-01</point>
      <point r="3.960">  0.63555906E-01 -0.94769292E-01 -0.13321019E+00  0.26140885E-01</point>
      <point r="3.980">  0.61869253E-01 -0.92506636E-01 -0.13046370E+00  0.25347042E-01</point>
      <point r="4.000">  0.60220058E-01 -0.90285507E-01 -0.12775295E+00  0.24574780E-01</point>
      <point r="4.020">  0.58607733E-01 -0.88105552E-01 -0.12507818E+00  0.23823599E-01</point>
      <point r="4.040">  0.57031689E-01 -0.85966407E-01 -0.12243962E+00  0.23093007E-01</point>
      <point r="4.060">  0.55491339E-01 -0.83867694E-01 -0.11983743E+00  0.22382522E-01</point>
      <point r="4.080">  0.53986100E-01 -0.81809026E-01 -0.11727178E+00  0.21691669E-01</point>
      <point r="4.100">  0.52515387E-01 -0.79790008E-01 -0.11474279E+00  0.21019982E-01</point>
      <point r="4.120">  0.51078621E-01 -0.77810232E-01 -0.11225055E+00  0.20367003E-01</point>
      <point r="4.140">  0.49675225E-01 -0.75869285E-01 -0.10979513E+00  0.19732284E-01</point>
      <point r="4.160">  0.48304623E-01 -0.73966744E-01 -0.10737658E+00  0.19115385E-01</point>
      <point r="4.180">  0.46966245E-01 -0.72102178E-01 -0.10499489E+00  0.18515873E-01</point>
      <point r="4.200">  0.45659523E-01 -0.70275151E-01 -0.10265007E+00  0.17933324E-01</point>
      <point r="4.220">  0.44383893E-01 -0.68485219E-01 -0.10034208E+00  0.17367322E-01</point>
      <point r="4.240">  0.43138794E-01 -0.66731935E-01 -0.98070865E-01  0.16817460E-01</point>
      <point r="4.260">  0.41923672E-01 -0.65014843E-01 -0.95836341E-01  0.16283338E-01</point>
      <point r="4.280">  0.40737973E-01 -0.63333484E-01 -0.93638411E-01  0.15764564E-01</point>
      <point r="4.300">  0.39581152E-01 -0.61687397E-01 -0.91476954E-01  0.15260754E-01</point>
      <point r="4.320">  0.38452665E-01 -0.60076113E-01 -0.89351832E-01  0.14771531E-01</point>
      <point r="4.340">  0.37351975E-01 -0.58499163E-01 -0.87262887E-01  0.14296527E-01</point>
      <point r="4.360">  0.36278549E-01 -0.56956073E-01 -0.85209943E-01  0.13835380E-01</point>
      <point r="4.380">  0.35231860E-01 -0.55446369E-01 -0.83192808E-01  0.13387737E-01</point>
      <point r="4.400">  0.34211386E-01 -0.53969573E-01 -0.81211275E-01  0.12953252E-01</point>
      <point r="4.420">  0.33216610E-01 -0.52525205E-01 -0.79265120E-01  0.12531584E-01</point>
      <point r="4.440">  0.32247020E-01 -0.51112786E-01 -0.77354105E-01  0.12122403E-01</point>
      <point r="4.460">  0.31302111E-01 -0.49731835E-01 -0.75477980E-01  0.11725383E-01</point>
      <point r="4.480">  0.30381383E-01 -0.48381871E-01 -0.73636478E-01  0.11340205E-01</point>
      <point r="4.500">  0.29484342E-01 -0.47062410E-01 -0.71829325E-01  0.10966560E-01</point>
      <point r="4.520">  0.28610499E-01 -0.45772973E-01 -0.70056231E-01  0.10604141E-01</point>
      <point r="4.540">  0.27759374E-01 -0.44513077E-01 -0.68316897E-01  0.10252653E-01</point>
      <point r="4.560">  0.26930489E-01 -0.43282244E-01 -0.66611014E-01  0.99118027E-02</point>
      <point r="4.580">  0.26123374E-01 -0.42079994E-01 -0.64938263E-01  0.95813062E-02</point>
      <point r="4.600">  0.25337566E-01 -0.40905849E-01 -0.63298315E-01  0.92608850E-02</point>
      <point r="4.620">  0.24572608E-01 -0.39759334E-01 -0.61690834E-01  0.89502670E-02</point>
      <point r="4.640">  0.23828047E-01 -0.38639974E-01 -0.60115475E-01  0.86491861E-02</point>
      <point r="4.660">  0.23103439E-01 -0.37547297E-01 -0.58571886E-01  0.83573826E-02</point>
      <point r="4.680">  0.22398346E-01 -0.36480835E-01 -0.57059709E-01  0.80746023E-02</point>
      <point r="4.700">  0.21712335E-01 -0.35440120E-01 -0.55578578E-01  0.78005973E-02</point>
      <point r="4.720">  0.21044981E-01 -0.34424689E-01 -0.54128124E-01  0.75351252E-02</point>
      <point r="4.740">  0.20395863E-01 -0.33434080E-01 -0.52707971E-01  0.72779495E-02</point>
      <point r="4.760">  0.19764569E-01 -0.32467836E-01 -0.51317737E-01  0.70288388E-02</point>
      <point r="4.780">  0.19150693E-01 -0.31525502E-01 -0.49957038E-01  0.67875678E-02</point>
      <point r="4.800">  0.18553834E-01 -0.30606628E-01 -0.48625486E-01  0.65539161E-02</point>
      <point r="4.820">  0.17973599E-01 -0.29710767E-01 -0.47322689E-01  0.63276686E-02</point>
      <point r="4.840">  0.17409600E-01 -0.28837475E-01 -0.46048251E-01  0.61086158E-02</point>
      <point r="4.860">  0.16861457E-01 -0.27986315E-01 -0.44801775E-01  0.58965528E-02</point>
      <point r="4.880">  0.16328795E-01 -0.27156850E-01 -0.43582862E-01  0.56912798E-02</point>
      <point r="4.900">  0.15811246E-01 -0.26348650E-01 -0.42391109E-01  0.54926022E-02</point>
      <point r="4.920">  0.15308448E-01 -0.25561290E-01 -0.41226114E-01  0.53003299E-02</point>
      <point r="4.940">  0.14820046E-01 -0.24794347E-01 -0.40087473E-01  0.51142775E-02</point>
      <point r="4.960">  0.14345691E-01 -0.24047405E-01 -0.38974780E-01  0.49342643E-02</point>
      <point r="4.980">  0.13885040E-01 -0.23320051E-01 -0.37887632E-01  0.47601143E-02</point>
      <point r="5.000">  0.13437756E-01 -0.22611878E-01 -0.36825621E-01  0.45916556E-02</point>
      <point r="5.020">  0.13003509E-01 -0.21922483E-01 -0.35788344E-01  0.44287210E-02</point>
      <point r="5.040">  0.12581974E-01 -0.21251469E-01 -0.34775395E-01  0.42711473E-02</point>
      <point r="5.060">  0.12172833E-01 -0.20598444E-01 -0.33786371E-01  0.41187755E-02</point>
      <point r="5.080">  0.11775775E-01 -0.19963019E-01 -0.32820868E-01  0.39714509E-02</point>
      <point r="5.100">  0.11390492E-01 -0.19344812E-01 -0.31878484E-01  0.38290226E-02</point>
      <point r="5.120">  0.11016684E-01 -0.18743447E-01 -0.30958820E-01  0.36913437E-02</point>
      <point r="5.140">  0.10654059E-01 -0.18158551E-01 -0.30061476E-01  0.35582711E-02</point>
      <point r="5.160">  0.10302325E-01 -0.17589758E-01 -0.29186055E-01  0.34296656E-02</point>
      <point r="5.180">  0.99612020E-02 -0.17036706E-01 -0.28332164E-01  0.33053915E-02</point>
      <point r="5.200">  0.96304117E-02 -0.16499039E-01 -0.27499410E-01  0.31853167E-02</point>
      <point r="5.220">  0.93096830E-02 -0.15976407E-01 -0.26687403E-01  0.30693129E-02</point>
      <point r="5.240">  0.89987502E-02 -0.15468464E-01 -0.25895755E-01  0.29572548E-02</point>
      <point r="5.260">  0.86973531E-02 -0.14974871E-01 -0.25124083E-01  0.28490209E-02</point>
      <point r="5.280">  0.84052368E-02 -0.14495293E-01 -0.24372006E-01  0.27444927E-02</point>
      <point r="5.300">  0.81221520E-02 -0.14029400E-01 -0.23639144E-01  0.26435550E-02</point>
      <point r="5.320">  0.78478548E-02 -0.13576869E-01 -0.22925124E-01  0.25460958E-02</point>
      <point r="5.340">  0.75821063E-02 -0.13137381E-01 -0.22229573E-01  0.24520063E-02</point>
      <point r="5.360">  0.73246730E-02 -0.12710623E-01 -0.21552124E-01  0.23611803E-02</point>
      <point r="5.380">  0.70753265E-02 -0.12296288E-01 -0.20892412E-01  0.22735149E-02</point>
      <point r="5.400">  0.68338434E-02 -0.11894073E-01 -0.20250077E-01  0.21889100E-02</point>
      <point r="5.420">  0.66000054E-02 -0.11503681E-01 -0.19624761E-01  0.21072682E-02</point>
      <point r="5.440">  0.63735989E-02 -0.11124820E-01 -0.19016113E-01  0.20284949E-02</point>
      <point r="5.460">  0.61544155E-02 -0.10757204E-01 -0.18423783E-01  0.19524982E-02</point>
      <point r="5.480">  0.59422512E-02 -0.10400551E-01 -0.17847427E-01  0.18791888E-02</point>
      <point r="5.500">  0.57369070E-02 -0.10054585E-01 -0.17286704E-01  0.18084798E-02</point>
      <point r="5.520">  0.55381883E-02 -0.97190359E-02 -0.16741278E-01  0.17402870E-02</point>
      <point r="5.540">  0.53459052E-02 -0.93936368E-02 -0.16210816E-01  0.16745286E-02</point>
      <point r="5.560">  0.51598723E-02 -0.90781275E-02 -0.15694992E-01  0.16111249E-02</point>
      <point r="5.580">  0.49799084E-02 -0.87722523E-02 -0.15193482E-01  0.15499988E-02</point>
      <point r="5.600">  0.48058369E-02 -0.84757605E-02 -0.14705968E-01  0.14910753E-02</point>
      <point r="5.620">  0.46374854E-02 -0.81884066E-02 -0.14232134E-01  0.14342817E-02</point>
      <point r="5.640">  0.44746856E-02 -0.79099497E-02 -0.13771671E-01  0.13795474E-02</point>
      <point r="5.660">  0.43172733E-02 -0.76401540E-02 -0.13324275E-01  0.13268037E-02</point>
      <point r="5.680">  0.41650886E-02 -0.73787884E-02 -0.12889643E-01  0.12759843E-02</point>
      <point r="5.700">  0.40179752E-02 -0.71256268E-02 -0.12467481E-01  0.12270246E-02</point>
      <point r="5.720">  0.38757811E-02 -0.68804475E-02 -0.12057496E-01  0.11798620E-02</point>
      <point r="5.740">  0.37383578E-02 -0.66430336E-02 -0.11659402E-01  0.11344358E-02</point>
      <point r="5.760">  0.36055607E-02 -0.64131730E-02 -0.11272916E-01  0.10906870E-02</point>
      <point r="5.780">  0.34772491E-02 -0.61906579E-02 -0.10897761E-01  0.10485587E-02</point>
      <point r="5.800">  0.33532854E-02 -0.59752851E-02 -0.10533663E-01  0.10079953E-02</point>
      <point r="5.820">  0.32335362E-02 -0.57668560E-02 -0.10180355E-01  0.96894340E-03</point>
      <point r="5.840">  0.31178710E-02 -0.55651760E-02 -0.98375709E-02  0.93135079E-03</point>
      <point r="5.860">  0.30061631E-02 -0.53700554E-02 -0.95050530E-02  0.89516711E-03</point>
      <point r="5.880">  0.28982891E-02 -0.51813082E-02 -0.91825461E-02  0.86034350E-03</point>
      <point r="5.900">  0.27941286E-02 -0.49987529E-02 -0.88698000E-02  0.82683261E-03</point>
      <point r="5.920">  0.26935648E-02 -0.48222122E-02 -0.85665689E-02  0.79458857E-03</point>
      <point r="5.940">  0.25964837E-02 -0.46515128E-02 -0.82726117E-02  0.76356697E-03</point>
      <point r="5.960">  0.25027747E-02 -0.44864853E-02 -0.79876918E-02  0.73372477E-03</point>
      <point r="5.980">  0.24123301E-02 -0.43269646E-02 -0.77115768E-02  0.70502031E-03</point>
      <point r="6.000">  0.23250449E-02 -0.41727891E-02 -0.74440389E-02  0.67741323E-03</point>
      <point r="6.020">  0.22408174E-02 -0.40238014E-02 -0.71848547E-02  0.65086448E-03</point>
      <point r="6.040">  0.21595485E-02 -0.38798478E-02 -0.69338051E-02  0.62533625E-03</point>
      <point r="6.060">  0.20811420E-02 -0.37407780E-02 -0.66906754E-02  0.60079194E-03</point>
      <point r="6.080">  0.20055041E-02 -0.36064459E-02 -0.64552551E-02  0.57719614E-03</point>
      <point r="6.100">  0.19325441E-02 -0.34767085E-02 -0.62273380E-02  0.55451458E-03</point>
      <point r="6.120">  0.18621735E-02 -0.33514266E-02 -0.60067220E-02  0.53271411E-03</point>
      <point r="6.140">  0.17943066E-02 -0.32304645E-02 -0.57932093E-02  0.51176266E-03</point>
      <point r="6.160">  0.17288600E-02 -0.31136898E-02 -0.55866062E-02  0.49162921E-03</point>
      <point r="6.180">  0.16657527E-02 -0.30009735E-02 -0.53867229E-02  0.47228377E-03</point>
      <point r="6.200">  0.16049062E-02 -0.28921897E-02 -0.51933739E-02  0.45369733E-03</point>
      <point r="6.220">  0.15462443E-02 -0.27872162E-02 -0.50063774E-02  0.43584186E-03</point>
      <point r="6.240">  0.14896928E-02 -0.26859334E-02 -0.48255559E-02  0.41869025E-03</point>
      <point r="6.260">  0.14351801E-02 -0.25882252E-02 -0.46507353E-02  0.40221631E-03</point>
      <point r="6.280">  0.13826363E-02 -0.24939785E-02 -0.44817459E-02  0.38639473E-03</point>
      <point r="6.300">  0.13319940E-02 -0.24030829E-02 -0.43184212E-02  0.37120104E-03</point>
      <point r="6.320">  0.12831875E-02 -0.23154314E-02 -0.41605989E-02  0.35661161E-03</point>
      <point r="6.340">  0.12361534E-02 -0.22309194E-02 -0.40081202E-02  0.34260362E-03</point>
      <point r="6.360">  0.11908299E-02 -0.21494455E-02 -0.38608299E-02  0.32915503E-03</point>
      <point r="6.380">  0.11471574E-02 -0.20709108E-02 -0.37185764E-02  0.31624454E-03</point>
      <point r="6.400">  0.11050780E-02 -0.19952193E-02 -0.35812117E-02  0.30385161E-03</point>
      <point r="6.420">  0.10645356E-02 -0.19222774E-02 -0.34485912E-02  0.29195640E-03</point>
      <point r="6.440">  0.10254759E-02 -0.18519943E-02 -0.33205738E-02  0.28053975E-03</point>
      <point r="6.460">  0.98784618E-03 -0.17842816E-02 -0.31970217E-02  0.26958318E-03</point>
      <point r="6.480">  0.95159558E-03 -0.17190535E-02 -0.30778004E-02  0.25906887E-03</point>
      <point r="6.500">  0.91667473E-03 -0.16562266E-02 -0.29627787E-02  0.24897959E-03</point>
      <point r="6.520">  0.88303588E-03 -0.15957197E-02 -0.28518287E-02  0.23929877E-03</point>
      <point r="6.540">  0.85063278E-03 -0.15374543E-02 -0.27448254E-02  0.23001038E-03</point>
      <point r="6.560">  0.81942070E-03 -0.14813536E-02 -0.26416472E-02  0.22109900E-03</point>
      <point r="6.580">  0.78935637E-03 -0.14273437E-02 -0.25421753E-02  0.21254974E-03</point>
      <point r="6.600">  0.76039792E-03 -0.13753522E-02 -0.24462940E-02  0.20434826E-03</point>
      <point r="6.620">  0.73250486E-03 -0.13253093E-02 -0.23538907E-02  0.19648072E-03</point>
      <point r="6.640">  0.70563804E-03 -0.12771471E-02 -0.22648553E-02  0.18893380E-03</point>
      <point r="6.660">  0.67975962E-03 -0.12307996E-02 -0.21790809E-02  0.18169465E-03</point>
      <point r="6.680">  0.65483301E-03 -0.11862030E-02 -0.20964631E-02  0.17475092E-03</point>
      <point r="6.700">  0.63082285E-03 -0.11432953E-02 -0.20169006E-02  0.16809067E-03</point>
      <point r="6.720">  0.60769500E-03 -0.11020164E-02 -0.19402943E-02  0.16170244E-03</point>
      <point r="6.740">  0.58541644E-03 -0.10623080E-02 -0.18665481E-02  0.15557517E-03</point>
      <point r="6.760">  0.56395532E-03 -0.10241136E-02 -0.17955684E-02  0.14969823E-03</point>
      <point r="6.780">  0.54328086E-03 -0.98737849E-03 -0.17272639E-02  0.14406137E-03</point>
      <point r="6.800">  0.52336336E-03 -0.95204963E-03 -0.16615460E-02  0.13865474E-03</point>
      <point r="6.820">  0.50417413E-03 -0.91807564E-03 -0.15983285E-02  0.13346885E-03</point>
      <point r="6.840">  0.48568551E-03 -0.88540677E-03 -0.15375275E-02  0.12849459E-03</point>
      <point r="6.860">  0.46787079E-03 -0.85399483E-03 -0.14790615E-02  0.12372317E-03</point>
      <point r="6.880">  0.45070421E-03 -0.82379319E-03 -0.14228512E-02  0.11914615E-03</point>
      <point r="6.900">  0.43416094E-03 -0.79475668E-03 -0.13688195E-02  0.11475541E-03</point>
      <point r="6.920">  0.41821703E-03 -0.76684163E-03 -0.13168915E-02  0.11054316E-03</point>
      <point r="6.940">  0.40284938E-03 -0.74000575E-03 -0.12669946E-02  0.10650189E-03</point>
      <point r="6.960">  0.38803573E-03 -0.71420817E-03 -0.12190581E-02  0.10262440E-03</point>
      <point r="6.980">  0.37375464E-03 -0.68940933E-03 -0.11730133E-02  0.98903756E-04</point>
      <point r="7.000">  0.35998545E-03 -0.66557101E-03 -0.11287936E-02  0.95333305E-04</point>
      <point r="7.020">  0.34670825E-03 -0.64265624E-03 -0.10863344E-02  0.91906657E-04</point>
      <point r="7.040">  0.33390387E-03 -0.62062931E-03 -0.10455730E-02  0.88617676E-04</point>
      <point r="7.060">  0.32155386E-03 -0.59945570E-03 -0.10064483E-02  0.85460468E-04</point>
      <point r="7.080">  0.30964045E-03 -0.57910206E-03 -0.96890126E-03  0.82429377E-04</point>
      <point r="7.100">  0.29814655E-03 -0.55953621E-03 -0.93287463E-03  0.79518974E-04</point>
      <point r="7.120">  0.28705571E-03 -0.54072703E-03 -0.89831280E-03  0.76724051E-04</point>
      <point r="7.140">  0.27635210E-03 -0.52264451E-03 -0.86516187E-03  0.74039608E-04</point>
      <point r="7.160">  0.26602049E-03 -0.50525969E-03 -0.83336961E-03  0.71460853E-04</point>
      <point r="7.180">  0.25604627E-03 -0.48854460E-03 -0.80288538E-03  0.68983189E-04</point>
      <point r="7.200">  0.24641534E-03 -0.47247229E-03 -0.77366013E-03  0.66602208E-04</point>
      <point r="7.220">  0.23711420E-03 -0.45701673E-03 -0.74564633E-03  0.64313686E-04</point>
      <point r="7.240">  0.22812985E-03 -0.44215285E-03 -0.71879797E-03  0.62113574E-04</point>
      <point r="7.260">  0.21944981E-03 -0.42785647E-03 -0.69307046E-03  0.59997992E-04</point>
      <point r="7.280">  0.21106208E-03 -0.41410430E-03 -0.66842067E-03  0.57963225E-04</point>
      <point r="7.300">  0.20295516E-03 -0.40087388E-03 -0.64480683E-03  0.56005711E-04</point>
      <point r="7.320">  0.19511799E-03 -0.38814359E-03 -0.62218853E-03  0.54122044E-04</point>
      <point r="7.340">  0.18753996E-03 -0.37589260E-03 -0.60052667E-03  0.52308961E-04</point>
      <point r="7.360">  0.18021089E-03 -0.36410086E-03 -0.57978344E-03  0.50563339E-04</point>
      <point r="7.380">  0.17312103E-03 -0.35274908E-03 -0.55992227E-03  0.48882189E-04</point>
      <point r="7.400">  0.16626099E-03 -0.34181868E-03 -0.54090781E-03  0.47262654E-04</point>
      <point r="7.420">  0.15962181E-03 -0.33129181E-03 -0.52270588E-03  0.45701998E-04</point>
      <point r="7.440">  0.15319486E-03 -0.32115127E-03 -0.50528347E-03  0.44197608E-04</point>
      <point r="7.460">  0.14697190E-03 -0.31138055E-03 -0.48860867E-03  0.42746985E-04</point>
      <point r="7.480">  0.14094501E-03 -0.30196378E-03 -0.47265066E-03  0.41347738E-04</point>
      <point r="7.500">  0.13510662E-03 -0.29288569E-03 -0.45737971E-03  0.39997587E-04</point>
      <point r="7.520">  0.12944948E-03 -0.28413164E-03 -0.44276707E-03  0.38694349E-04</point>
      <point r="7.540">  0.12396663E-03 -0.27568754E-03 -0.42878504E-03  0.37435943E-04</point>
      <point r="7.560">  0.11865143E-03 -0.26753988E-03 -0.41540686E-03  0.36220381E-04</point>
      <point r="7.580">  0.11349750E-03 -0.25967569E-03 -0.40260673E-03  0.35045763E-04</point>
      <point r="7.600">  0.10849875E-03 -0.25208254E-03 -0.39035977E-03  0.33910279E-04</point>
      <point r="7.620">  0.10364935E-03 -0.24474848E-03 -0.37864197E-03  0.32812200E-04</point>
      <point r="7.640">  0.98943732E-04 -0.23766207E-03 -0.36743023E-03  0.31749878E-04</point>
      <point r="7.660">  0.94376561E-04 -0.23081234E-03 -0.35670225E-03  0.30721740E-04</point>
      <point r="7.680">  0.89942740E-04 -0.22418878E-03 -0.34643657E-03  0.29726289E-04</point>
      <point r="7.700">  0.85637398E-04 -0.21778132E-03 -0.33661251E-03  0.28762094E-04</point>
      <point r="7.720">  0.81455879E-04 -0.21158031E-03 -0.32721017E-03  0.27827796E-04</point>
      <point r="7.740">  0.77393736E-04 -0.20557653E-03 -0.31821041E-03  0.26922097E-04</point>
      <point r="7.760">  0.73446719E-04 -0.19976114E-03 -0.30959477E-03  0.26043764E-04</point>
      <point r="7.780">  0.69610773E-04 -0.19412570E-03 -0.30134553E-03  0.25191620E-04</point>
      <point r="7.800">  0.65882021E-04 -0.18866213E-03 -0.29344564E-03  0.24364546E-04</point>
      <point r="7.820">  0.62256766E-04 -0.18336270E-03 -0.28587870E-03  0.23561478E-04</point>
      <point r="7.840">  0.58731479E-04 -0.17822004E-03 -0.27862898E-03  0.22781402E-04</point>
      <point r="7.860">  0.55302789E-04 -0.17322710E-03 -0.27168133E-03  0.22023357E-04</point>
      <point r="7.880">  0.51967484E-04 -0.16837716E-03 -0.26502121E-03  0.21286425E-04</point>
      <point r="7.900">  0.48722497E-04 -0.16366379E-03 -0.25863469E-03  0.20569736E-04</point>
      <point r="7.920">  0.45564903E-04 -0.15908088E-03 -0.25250837E-03  0.19872463E-04</point>
      <point r="7.940">  0.42491914E-04 -0.15462257E-03 -0.24662940E-03  0.19193820E-04</point>
      <point r="7.960">  0.39500869E-04 -0.15028331E-03 -0.24098546E-03  0.18533060E-04</point>
      <point r="7.980">  0.36589234E-04 -0.14605780E-03 -0.23556475E-03  0.17889476E-04</point>
      <point r="8.000">  0.33754590E-04 -0.14194098E-03 -0.23035596E-03  0.17262395E-04</point>
      <point r="8.020">  0.30994634E-04 -0.13792804E-03 -0.22534823E-03  0.16651177E-04</point>
      <point r="8.040">  0.28307168E-04 -0.13401443E-03 -0.22053119E-03  0.16055219E-04</point>
      <point r="8.060">  0.25690099E-04 -0.13019579E-03 -0.21589491E-03  0.15473944E-04</point>
      <point r="8.080">  0.23141431E-04 -0.12646798E-03 -0.21142989E-03  0.14906809E-04</point>
      <point r="8.100">  0.20659263E-04 -0.12282710E-03 -0.20712703E-03  0.14353298E-04</point>
      <point r="8.120">  0.18241783E-04 -0.11926940E-03 -0.20297765E-03  0.13812921E-04</point>
      <point r="8.140">  0.15887261E-04 -0.11579136E-03 -0.19897345E-03  0.13285214E-04</point>
      <point r="8.160">  0.13594052E-04 -0.11238962E-03 -0.19510651E-03  0.12769738E-04</point>
      <point r="8.180">  0.11360586E-04 -0.10906100E-03 -0.19136925E-03  0.12266078E-04</point>
      <point r="8.200">  0.91853652E-05 -0.10580250E-03 -0.18775447E-03  0.11773838E-04</point>
      <point r="8.220">  0.70669637E-05 -0.10261126E-03 -0.18425528E-03  0.11292648E-04</point>
      <point r="8.240">  0.50040201E-05 -0.99484584E-04 -0.18086513E-03  0.10822153E-04</point>
      <point r="8.260">  0.29952363E-05 -0.96419930E-04 -0.17757776E-03  0.10362019E-04</point>
      <point r="8.280">  0.10393735E-05 -0.93414884E-04 -0.17438723E-03  0.99119301E-05</point>
      <point r="8.300"> -0.86475052E-06 -0.90467173E-04 -0.17128790E-03  0.94715877E-05</point>
      <point r="8.320"> -0.27182644E-05 -0.87574649E-04 -0.16827437E-03  0.90407084E-05</point>
      <point r="8.340"> -0.45222461E-05 -0.84735288E-04 -0.16534155E-03  0.86190245E-05</point>
      <point r="8.360"> -0.62777258E-05 -0.81947184E-04 -0.16248457E-03  0.82062830E-05</point>
      <point r="8.380"> -0.79856881E-05 -0.79208539E-04 -0.15969885E-03  0.78022442E-05</point>
      <point r="8.400"> -0.96470747E-05 -0.76517664E-04 -0.15698001E-03  0.74066817E-05</point>
      <point r="8.420"> -0.11262787E-04 -0.73872970E-04 -0.15432391E-03  0.70193810E-05</point>
      <point r="8.440"> -0.12833688E-04 -0.71272966E-04 -0.15172666E-03  0.66401395E-05</point>
      <point r="8.460"> -0.14360604E-04 -0.68716250E-04 -0.14918453E-03  0.62687653E-05</point>
      <point r="8.480"> -0.15844329E-04 -0.66201509E-04 -0.14669405E-03  0.59050769E-05</point>
      <point r="8.500"> -0.17285624E-04 -0.63727512E-04 -0.14425189E-03  0.55489024E-05</point>
      <point r="8.520"> -0.18685219E-04 -0.61293108E-04 -0.14185495E-03  0.52000792E-05</point>
      <point r="8.540"> -0.20043817E-04 -0.58897218E-04 -0.13950029E-03  0.48584532E-05</point>
      <point r="8.560"> -0.21362094E-04 -0.56538836E-04 -0.13718514E-03  0.45238785E-05</point>
      <point r="8.580"> -0.22640703E-04 -0.54217023E-04 -0.13490691E-03  0.41962166E-05</point>
      <point r="8.600"> -0.23880269E-04 -0.51930903E-04 -0.13266315E-03  0.38753365E-05</point>
      <point r="8.620"> -0.25081400E-04 -0.49679659E-04 -0.13045158E-03  0.35611135E-05</point>
      <point r="8.640"> -0.26244680E-04 -0.47462533E-04 -0.12827005E-03  0.32534294E-05</point>
      <point r="8.660"> -0.27370678E-04 -0.45278818E-04 -0.12611655E-03  0.29521720E-05</point>
      <point r="8.680"> -0.28459940E-04 -0.43127859E-04 -0.12398920E-03  0.26572343E-05</point>
      <point r="8.700"> -0.29513000E-04 -0.41009050E-04 -0.12188627E-03  0.23685148E-05</point>
      <point r="8.720"> -0.30530374E-04 -0.38921826E-04 -0.11980612E-03  0.20859164E-05</point>
      <point r="8.740"> -0.31512565E-04 -0.36865668E-04 -0.11774724E-03  0.18093469E-05</point>
      <point r="8.760"> -0.32460061E-04 -0.34840095E-04 -0.11570823E-03  0.15387179E-05</point>
      <point r="8.780"> -0.33373340E-04 -0.32844662E-04 -0.11368780E-03  0.12739452E-05</point>
      <point r="8.800"> -0.34252867E-04 -0.30878961E-04 -0.11168474E-03  0.10149479E-05</point>
      <point r="8.820"> -0.35099098E-04 -0.28942615E-04 -0.10969796E-03  0.76164858E-06</point>
      <point r="8.840"> -0.35912476E-04 -0.27035278E-04 -0.10772643E-03  0.51397299E-06</point>
      <point r="8.860"> -0.36693440E-04 -0.25156630E-04 -0.10576925E-03  0.27184961E-06</point>
      <point r="8.880"> -0.37442416E-04 -0.23306382E-04 -0.10382556E-03  0.35209608E-07</point>
      <point r="8.900"> -0.38159827E-04 -0.21484264E-04 -0.10189459E-03 -0.19601344E-06</point>
      <point r="8.920"> -0.38846085E-04 -0.19690033E-04 -0.99975666E-04 -0.42188371E-06</point>
      <point r="8.940"> -0.39501599E-04 -0.17923465E-04 -0.98068150E-04 -0.64246335E-06</point>
      <point r="8.960"> -0.40126770E-04 -0.16184354E-04 -0.96171490E-04 -0.85781262E-06</point>
      <point r="8.980"> -0.40721996E-04 -0.14472513E-04 -0.94285190E-04 -0.10679901E-05</point>
      <point r="9.000"> -0.41287670E-04 -0.12787773E-04 -0.92408817E-04 -0.12730528E-05</point>
      <point r="9.020"> -0.41824178E-04 -0.11129977E-04 -0.90541991E-04 -0.14730562E-05</point>
      <point r="9.040"> -0.42331906E-04 -0.94989818E-05 -0.88684385E-04 -0.16680549E-05</point>
      <point r="9.060"> -0.42811236E-04 -0.78946580E-05 -0.86835722E-04 -0.18581018E-05</point>
      <point r="9.080"> -0.43262544E-04 -0.63168858E-05 -0.84995771E-04 -0.20432492E-05</point>
      <point r="9.100"> -0.43686207E-04 -0.47655558E-05 -0.83164345E-04 -0.22235483E-05</point>
      <point r="9.120"> -0.44082599E-04 -0.32405671E-05 -0.81341299E-04 -0.23990496E-05</point>
      <point r="9.140"> -0.44452089E-04 -0.17418266E-05 -0.79526525E-04 -0.25698027E-05</point>
      <point r="9.160"> -0.44795049E-04 -0.26924796E-06 -0.77719952E-04 -0.27358568E-05</point>
      <point r="9.180"> -0.45111846E-04  0.11772495E-05 -0.75921542E-04 -0.28972603E-05</point>
      <point r="9.200"> -0.45402847E-04  0.25977412E-05 -0.74131289E-04 -0.30540613E-05</point>
      <point r="9.220"> -0.45668418E-04  0.39922988E-05 -0.72349216E-04 -0.32063076E-05</point>
      <point r="9.240"> -0.45908923E-04  0.53609904E-05 -0.70575371E-04 -0.33540463E-05</point>
      <point r="9.260"> -0.46124727E-04  0.67038814E-05 -0.68809831E-04 -0.34973246E-05</point>
      <point r="9.280"> -0.46316193E-04  0.80210355E-05 -0.67052693E-04 -0.36361893E-05</point>
      <point r="9.300"> -0.46483684E-04  0.93125147E-05 -0.65304076E-04 -0.37706871E-05</point>
      <point r="9.320"> -0.46627563E-04  0.10578381E-04 -0.63564118E-04 -0.39008646E-05</point>
      <point r="9.340"> -0.46748192E-04  0.11818695E-04 -0.61832976E-04 -0.40267681E-05</point>
      <point r="9.360"> -0.46845934E-04  0.13033518E-04 -0.60110823E-04 -0.41484441E-05</point>
      <point r="9.380"> -0.46921149E-04  0.14222914E-04 -0.58397846E-04 -0.42659390E-05</point>
      <point r="9.400"> -0.46974202E-04  0.15386946E-04 -0.56694245E-04 -0.43792993E-05</point>
      <point r="9.420"> -0.47005452E-04  0.16525679E-04 -0.55000232E-04 -0.44885713E-05</point>
      <point r="9.440"> -0.47015263E-04  0.17639182E-04 -0.53316031E-04 -0.45938017E-05</point>
      <point r="9.460"> -0.47003997E-04  0.18727524E-04 -0.51641874E-04 -0.46950370E-05</point>
      <point r="9.480"> -0.46972015E-04  0.19790778E-04 -0.49978000E-04 -0.47923240E-05</point>
      <point r="9.500"> -0.46919679E-04  0.20829019E-04 -0.48324656E-04 -0.48857095E-05</point>
      <point r="9.520"> -0.46847351E-04  0.21842328E-04 -0.46682096E-04 -0.49752405E-05</point>
      <point r="9.540"> -0.46755394E-04  0.22830787E-04 -0.45050577E-04 -0.50609641E-05</point>
      <point r="9.560"> -0.46644169E-04  0.23794482E-04 -0.43430360E-04 -0.51429276E-05</point>
      <point r="9.580"> -0.46514039E-04  0.24733505E-04 -0.41821710E-04 -0.52211784E-05</point>
      <point r="9.600"> -0.46365364E-04  0.25647950E-04 -0.40224894E-04 -0.52957641E-05</point>
      <point r="9.620"> -0.46198508E-04  0.26537916E-04 -0.38640178E-04 -0.53667324E-05</point>
      <point r="9.640"> -0.46013831E-04  0.27403506E-04 -0.37067831E-04 -0.54341312E-05</point>
      <point r="9.660"> -0.45811694E-04  0.28244830E-04 -0.35508122E-04 -0.54980086E-05</point>
      <point r="9.680"> -0.45592459E-04  0.29062000E-04 -0.33961316E-04 -0.55584128E-05</point>
      <point r="9.700"> -0.45356486E-04  0.29855133E-04 -0.32427679E-04 -0.56153921E-05</point>
      <point r="9.720"> -0.45104135E-04  0.30624352E-04 -0.30907474E-04 -0.56689951E-05</point>
      <point r="9.740"> -0.44835766E-04  0.31369785E-04 -0.29400961E-04 -0.57192704E-05</point>
      <point r="9.760"> -0.44551738E-04  0.32091565E-04 -0.27908398E-04 -0.57662668E-05</point>
      <point r="9.780"> -0.44252409E-04  0.32789827E-04 -0.26430038E-04 -0.58100332E-05</point>
      <point r="9.800"> -0.43938136E-04  0.33464715E-04 -0.24966129E-04 -0.58506185E-05</point>
      <point r="9.820"> -0.43609277E-04  0.34116376E-04 -0.23516915E-04 -0.58880720E-05</point>
      <point r="9.840"> -0.43266187E-04  0.34744963E-04 -0.22082637E-04 -0.59224429E-05</point>
      <point r="9.860"> -0.42909222E-04  0.35350631E-04 -0.20663527E-04 -0.59537803E-05</point>
      <point r="9.880"> -0.42538734E-04  0.35933544E-04 -0.19259814E-04 -0.59821338E-05</point>
      <point r="9.900"> -0.42155077E-04  0.36493869E-04 -0.17871719E-04 -0.60075526E-05</point>
      <point r="9.920"> -0.41758602E-04  0.37031776E-04 -0.16499458E-04 -0.60300863E-05</point>
      <point r="9.940"> -0.41349659E-04  0.37547442E-04 -0.15143239E-04 -0.60497843E-05</point>
      <point r="9.960"> -0.40928596E-04  0.38041049E-04 -0.13803265E-04 -0.60666962E-05</point>
      <point r="9.980"> -0.40495761E-04  0.38512783E-04 -0.12479729E-04 -0.60808714E-05</point>
      <point r="10.000"> -0.40051499E-04  0.38962833E-04 -0.11172821E-04 -0.60923594E-05</point>
      <point r="10.020"> -0.39596153E-04  0.39391393E-04 -0.98827191E-05 -0.61012097E-05</point>
      <point r="10.040"> -0.39130066E-04  0.39798664E-04 -0.86095964E-05 -0.61074717E-05</point>
      <point r="10.060"> -0.38653578E-04  0.40184848E-04 -0.73536176E-05 -0.61111947E-05</point>
      <point r="10.080"> -0.38167026E-04  0.40550153E-04 -0.61149396E-05 -0.61124282E-05</point>
      <point r="10.100"> -0.37670748E-04  0.40894790E-04 -0.48937112E-05 -0.61112211E-05</point>
      <point r="10.120"> -0.37165076E-04  0.41218974E-04 -0.36900734E-05 -0.61076227E-05</point>
      <point r="10.140"> -0.36650342E-04  0.41522924E-04 -0.25041590E-05 -0.61016819E-05</point>
      <point r="10.160"> -0.36126877E-04  0.41806863E-04 -0.13360926E-05 -0.60934475E-05</point>
      <point r="10.180"> -0.35595006E-04  0.42071017E-04 -0.18599093E-06 -0.60829683E-05</point>
      <point r="10.200"> -0.35055055E-04  0.42315616E-04  0.94603780E-06 -0.60702928E-05</point>
      <point r="10.220"> -0.34507345E-04  0.42540892E-04  0.20598934E-05 -0.60554694E-05</point>
      <point r="10.240"> -0.33952196E-04  0.42747083E-04  0.31554839E-05 -0.60385462E-05</point>
      <point r="10.260"> -0.33389925E-04  0.42934426E-04  0.42327254E-05 -0.60195712E-05</point>
      <point r="10.280"> -0.32820846E-04  0.43103165E-04  0.52915423E-05 -0.59985921E-05</point>
      <point r="10.300"> -0.32245269E-04  0.43253544E-04  0.63318671E-05 -0.59756565E-05</point>
      <point r="10.320"> -0.31663503E-04  0.43385810E-04  0.73536401E-05 -0.59508116E-05</point>
      <point r="10.340"> -0.31075854E-04  0.43500213E-04  0.83568099E-05 -0.59241045E-05</point>
      <point r="10.360"> -0.30482623E-04  0.43597005E-04  0.93413328E-05 -0.58955818E-05</point>
      <point r="10.380"> -0.29884109E-04  0.43676442E-04  0.10307173E-04 -0.58652901E-05</point>
      <point r="10.400">  0.00000000E+00</point>
    </S_spline>
    <Vrep_spline>
      <point r="0.000">  0.27342623E+03</point>
      <point r="0.100">  0.20504684E+03</point>
      <point r="0.200">  0.15357627E+03</point>
      <point r="0.300">  0.11483334E+03</point>
      <point r="0.400">  0.85670737E+02</point>
      <point r="0.500">  0.63719456E+02</point>
      <point r="0.600">  0.47196280E+02</point>
      <point r="0.700">  0.34758949E+02</point>
      <point r="0.800">  0.25397118E+02</point>
      <point r="0.900">  0.18350277E+02</point>
      <point r="1.000">  0.13045976E+02</point>
      <point r="1.057">  0.10655245E+02</point>
      <point r="1.114">  0.86843958E+01</point>
      <point r="1.170">  0.70663867E+01</point>
      <point r="1.227">  0.57434070E+01</point>
      <point r="1.284">  0.46658413E+01</point>
      <point r="1.341">  0.37913209E+01</point>
      <point r="1.397">  0.30838599E+01</point>
      <point r="1.454">  0.25130679E+01</point>
      <point r="1.511">  0.20534375E+01</point>
      <point r="1.568">  0.16836989E+01</point>
      <point r="1.624">  0.13862412E+01</point>
      <point r="1.681">  0.11465913E+01</point>
      <point r="1.738">  0.95295090E+00</point>
      <point r="1.795">  0.79578377E+00</point>
      <point r="1.852">  0.66745221E+00</point>
      <point r="1.908">  0.56189727E+00</point>
      <point r="1.965">  0.47435986E+00</point>
      <point r="2.022">  0.40113923E+00</point>
      <point r="2.079">  0.33938539E+00</point>
      <point r="2.135">  0.28692244E+00</point>
      <point r="2.192">  0.24209978E+00</point>
      <point r="2.249">  0.20366827E+00</point>
      <point r="2.306">  0.17067870E+00</point>
      <point r="2.362">  0.14239994E+00</point>
      <point r="2.419">  0.11825438E+00</point>
      <point r="2.476">  0.97768319E-01</point>
      <point r="2.533">  0.80535170E-01</point>
      <point r="2.590">  0.66189522E-01</point>
      <point r="2.646">  0.54390131E-01</point>
      <point r="2.703">  0.44810172E-01</point>
      <point r="2.760">  0.37133141E-01</point>
      <point r="2.817">  0.31053017E-01</point>
      <point r="2.873">  0.26277351E-01</point>
      <point r="2.930">  0.22532178E-01</point>
      <point r="2.987">  0.19567720E-01</point>
      <point r="3.044">  0.17164050E-01</point>
      <point r="3.100">  0.15135976E-01</point>
      <point r="3.157">  0.13336593E-01</point>
      <point r="3.214">  0.11659073E-01</point>
      <point r="3.271">  0.10036399E-01</point>
      <point r="3.328">  0.84389139E-02</point>
      <point r="3.384">  0.68696777E-02</point>
      <point r="3.441">  0.53577890E-02</point>
      <point r="3.498">  0.39499487E-02</point>
      <point r="3.555">  0.27007077E-02</point>
      <point r="3.611">  0.16619699E-02</point>
      <point r="3.668">  0.87247140E-03</point>
      <point r="3.725">  0.34809868E-03</point>
      <point r="3.782">  0.74051485E-04</point>
      <point r="3.838">  0.00000000E+00</point>
    </Vrep_spline>
  </per_pair_data>
  <per_pair_data type1="1" type2="2" SK_cutoff="10.4000000000000004"
    Vrep_cutoff="3.9226500000000000" SK_npts="521" Vrep_npts="61">
    <H_spline>
      <point r="0.000">  0.39226500E+01  0.00000000E+00  0.11682400E+01 -0.42589700E+00  0.61469600E-01  0.10886800E+01 -0.15677300E+01  0.00000000E+00  0.39891300E-01 -0.29448900E+00</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400"> -0.87819926E+00  0.82689748E+00 -0.41515211E+00 -0.60566109E-01</point>
      <point r="0.420"> -0.86277382E+00  0.79795758E+00 -0.44333612E+00 -0.78141717E-01</point>
      <point r="0.440"> -0.84534046E+00  0.76799390E+00 -0.46938649E+00 -0.95351112E-01</point>
      <point r="0.460"> -0.82636805E+00  0.73744479E+00 -0.49327114E+00 -0.11214844E+00</point>
      <point r="0.480"> -0.80627264E+00  0.70668840E+00 -0.51497917E+00 -0.12849461E+00</point>
      <point r="0.500"> -0.78541999E+00  0.67604837E+00 -0.53451735E+00 -0.14435665E+00</point>
      <point r="0.520"> -0.76412878E+00  0.64579919E+00 -0.55190696E+00 -0.15970719E+00</point>
      <point r="0.540"> -0.74267422E+00  0.61617181E+00 -0.56718169E+00 -0.17452396E+00</point>
      <point r="0.560"> -0.72129134E+00  0.58735812E+00 -0.58038549E+00 -0.18878925E+00</point>
      <point r="0.580"> -0.70017849E+00  0.55951425E+00 -0.59157070E+00 -0.20248945E+00</point>
      <point r="0.600"> -0.67950210E+00  0.53276565E+00 -0.60079893E+00 -0.21561490E+00</point>
      <point r="0.620"> -0.65939994E+00  0.50721086E+00 -0.60813930E+00 -0.22815919E+00</point>
      <point r="0.640"> -0.63998286E+00  0.48292313E+00 -0.61366366E+00 -0.24011863E+00</point>
      <point r="0.660"> -0.62133803E+00  0.45995409E+00 -0.61744606E+00 -0.25149234E+00</point>
      <point r="0.680"> -0.60353313E+00  0.43833881E+00 -0.61956525E+00 -0.26228216E+00</point>
      <point r="0.700"> -0.58661858E+00  0.41809813E+00 -0.62010387E+00 -0.27249184E+00</point>
      <point r="0.720"> -0.57062861E+00  0.39923877E+00 -0.61914437E+00 -0.28212669E+00</point>
      <point r="0.740"> -0.55558297E+00  0.38175424E+00 -0.61676712E+00 -0.29119372E+00</point>
      <point r="0.760"> -0.54148947E+00  0.36562722E+00 -0.61305190E+00 -0.29970171E+00</point>
      <point r="0.780"> -0.52834646E+00  0.35083240E+00 -0.60808030E+00 -0.30766070E+00</point>
      <point r="0.800"> -0.51614466E+00  0.33733845E+00 -0.60193605E+00 -0.31508170E+00</point>
      <point r="0.820"> -0.50487097E+00  0.32510888E+00 -0.59471334E+00 -0.32197673E+00</point>
      <point r="0.840"> -0.49449886E+00  0.31410246E+00 -0.58647387E+00 -0.32835774E+00</point>
      <point r="0.860"> -0.48500150E+00  0.30427312E+00 -0.57730004E+00 -0.33423841E+00</point>
      <point r="0.880"> -0.47634752E+00  0.29557186E+00 -0.56726776E+00 -0.33963257E+00</point>
      <point r="0.900"> -0.46850230E+00  0.28794762E+00 -0.55645056E+00 -0.34455445E+00</point>
      <point r="0.920"> -0.46142861E+00  0.28134811E+00 -0.54491913E+00 -0.34901858E+00</point>
      <point r="0.940"> -0.45508739E+00  0.27572024E+00 -0.53274138E+00 -0.35303975E+00</point>
      <point r="0.960"> -0.44943844E+00  0.27101035E+00 -0.51998274E+00 -0.35663294E+00</point>
      <point r="0.980"> -0.44444103E+00  0.26716456E+00 -0.50670631E+00 -0.35981321E+00</point>
      <point r="1.000"> -0.44005424E+00  0.26412918E+00 -0.49297274E+00 -0.36259567E+00</point>
      <point r="1.020"> -0.43623719E+00  0.26185118E+00 -0.47883993E+00 -0.36499538E+00</point>
      <point r="1.040"> -0.43294933E+00  0.26027852E+00 -0.46436279E+00 -0.36702730E+00</point>
      <point r="1.060"> -0.43015064E+00  0.25936040E+00 -0.44959325E+00 -0.36870627E+00</point>
      <point r="1.080"> -0.42780194E+00  0.25904737E+00 -0.43458038E+00 -0.37004696E+00</point>
      <point r="1.100"> -0.42586515E+00  0.25929139E+00 -0.41937061E+00 -0.37106387E+00</point>
      <point r="1.120"> -0.42430341E+00  0.26004596E+00 -0.40400787E+00 -0.37177128E+00</point>
      <point r="1.140"> -0.42308120E+00  0.26126623E+00 -0.38853363E+00 -0.37218317E+00</point>
      <point r="1.160"> -0.42216440E+00  0.26290908E+00 -0.37298687E+00 -0.37231328E+00</point>
      <point r="1.180"> -0.42152034E+00  0.26493320E+00 -0.35740408E+00 -0.37217501E+00</point>
      <point r="1.200"> -0.42111783E+00  0.26729911E+00 -0.34181929E+00 -0.37178146E+00</point>
      <point r="1.220"> -0.42092723E+00  0.26996915E+00 -0.32626406E+00 -0.37114540E+00</point>
      <point r="1.240"> -0.42092043E+00  0.27290751E+00 -0.31076765E+00 -0.37027929E+00</point>
      <point r="1.260"> -0.42107087E+00  0.27608019E+00 -0.29535707E+00 -0.36919526E+00</point>
      <point r="1.280"> -0.42135355E+00  0.27945503E+00 -0.28005730E+00 -0.36790506E+00</point>
      <point r="1.300"> -0.42174495E+00  0.28300167E+00 -0.26489131E+00 -0.36642007E+00</point>
      <point r="1.320"> -0.42222308E+00  0.28669154E+00 -0.24988027E+00 -0.36475130E+00</point>
      <point r="1.340"> -0.42276743E+00  0.29049786E+00 -0.23504353E+00 -0.36290942E+00</point>
      <point r="1.360"> -0.42335894E+00  0.29439557E+00 -0.22039877E+00 -0.36090473E+00</point>
      <point r="1.380"> -0.42397994E+00  0.29836129E+00 -0.20596196E+00 -0.35874718E+00</point>
      <point r="1.400"> -0.42461409E+00  0.30237320E+00 -0.19174745E+00 -0.35644639E+00</point>
      <point r="1.420"> -0.42524635E+00  0.30641094E+00 -0.17776803E+00 -0.35401161E+00</point>
      <point r="1.440"> -0.42586282E+00  0.31045555E+00 -0.16403500E+00 -0.35145175E+00</point>
      <point r="1.460"> -0.42645079E+00  0.31448938E+00 -0.15055828E+00 -0.34877535E+00</point>
      <point r="1.480"> -0.42699864E+00  0.31849609E+00 -0.13734652E+00 -0.34599059E+00</point>
      <point r="1.500"> -0.42749584E+00  0.32246063E+00 -0.12440715E+00 -0.34310537E+00</point>
      <point r="1.520"> -0.42793291E+00  0.32636919E+00 -0.11174654E+00 -0.34012725E+00</point>
      <point r="1.540"> -0.42830143E+00  0.33020926E+00 -0.99370034E-01 -0.33706349E+00</point>
      <point r="1.560"> -0.42859393E+00  0.33396950E+00 -0.87282032E-01 -0.33392108E+00</point>
      <point r="1.580"> -0.42880385E+00  0.33763973E+00 -0.75486051E-01 -0.33070671E+00</point>
      <point r="1.600"> -0.42892550E+00  0.34121081E+00 -0.63984773E-01 -0.32742677E+00</point>
      <point r="1.620"> -0.42895394E+00  0.34467458E+00 -0.52780092E-01 -0.32408735E+00</point>
      <point r="1.640"> -0.42888492E+00  0.34802374E+00 -0.41873153E-01 -0.32069427E+00</point>
      <point r="1.660"> -0.42870293E+00  0.35123163E+00 -0.31308214E-01 -0.31725282E+00</point>
      <point r="1.680"> -0.42843228E+00  0.35433701E+00 -0.20988447E-01 -0.31376869E+00</point>
      <point r="1.700"> -0.42805468E+00  0.35731043E+00 -0.10965758E-01 -0.31024664E+00</point>
      <point r="1.720"> -0.42756813E+00  0.36014741E+00 -0.12391277E-02 -0.30669146E+00</point>
      <point r="1.740"> -0.42697121E+00  0.36284408E+00  0.81928854E-02 -0.30310770E+00</point>
      <point r="1.760"> -0.42626294E+00  0.36539713E+00  0.17332137E-01 -0.29949969E+00</point>
      <point r="1.780"> -0.42544284E+00  0.36780389E+00  0.26180898E-01 -0.29587153E+00</point>
      <point r="1.800"> -0.42451079E+00  0.37006220E+00  0.34741861E-01 -0.29222712E+00</point>
      <point r="1.820"> -0.42346704E+00  0.37217046E+00  0.43018125E-01 -0.28857014E+00</point>
      <point r="1.840"> -0.42231218E+00  0.37412759E+00  0.51013185E-01 -0.28490408E+00</point>
      <point r="1.860"> -0.42104709E+00  0.37593297E+00  0.58730899E-01 -0.28123223E+00</point>
      <point r="1.880"> -0.41967291E+00  0.37758642E+00  0.66175452E-01 -0.27755772E+00</point>
      <point r="1.900"> -0.41819106E+00  0.37908816E+00  0.73351305E-01 -0.27388347E+00</point>
      <point r="1.920"> -0.41660317E+00  0.38043877E+00  0.80263150E-01 -0.27021225E+00</point>
      <point r="1.940"> -0.41491110E+00  0.38163914E+00  0.86915855E-01 -0.26654668E+00</point>
      <point r="1.960"> -0.41311692E+00  0.38269046E+00  0.93314422E-01 -0.26288920E+00</point>
      <point r="1.980"> -0.41122285E+00  0.38359418E+00  0.99463941E-01 -0.25924213E+00</point>
      <point r="2.000"> -0.40923129E+00  0.38435197E+00  0.10536956E+00 -0.25560762E+00</point>
      <point r="2.020"> -0.40714477E+00  0.38496568E+00  0.11103647E+00 -0.25198771E+00</point>
      <point r="2.040"> -0.40496592E+00  0.38543738E+00  0.11646987E+00 -0.24838429E+00</point>
      <point r="2.060"> -0.40269748E+00  0.38576927E+00  0.12167497E+00 -0.24479914E+00</point>
      <point r="2.080"> -0.40034224E+00  0.38596370E+00  0.12665700E+00 -0.24123391E+00</point>
      <point r="2.100"> -0.39790307E+00  0.38602317E+00  0.13142118E+00 -0.23769015E+00</point>
      <point r="2.120"> -0.39538288E+00  0.38595028E+00  0.13597274E+00 -0.23416927E+00</point>
      <point r="2.140"> -0.39278461E+00  0.38574775E+00  0.14031693E+00 -0.23067260E+00</point>
      <point r="2.160"> -0.39011125E+00  0.38541839E+00  0.14445898E+00 -0.22720137E+00</point>
      <point r="2.180"> -0.38736579E+00  0.38496511E+00  0.14840415E+00 -0.22375672E+00</point>
      <point r="2.200"> -0.38455125E+00  0.38439090E+00  0.15215768E+00 -0.22033968E+00</point>
      <point r="2.220"> -0.38167065E+00  0.38369883E+00  0.15572479E+00 -0.21695123E+00</point>
      <point r="2.240"> -0.37872704E+00  0.38289202E+00  0.15911069E+00 -0.21359223E+00</point>
      <point r="2.260"> -0.37572344E+00  0.38197365E+00  0.16232059E+00 -0.21026349E+00</point>
      <point r="2.280"> -0.37266288E+00  0.38094693E+00  0.16535961E+00 -0.20696574E+00</point>
      <point r="2.300"> -0.36954838E+00  0.37981512E+00  0.16823285E+00 -0.20369966E+00</point>
      <point r="2.320"> -0.36638292E+00  0.37858148E+00  0.17094533E+00 -0.20046582E+00</point>
      <point r="2.340"> -0.36316948E+00  0.37724928E+00  0.17350202E+00 -0.19726476E+00</point>
      <point r="2.360"> -0.35991097E+00  0.37582179E+00  0.17590776E+00 -0.19409696E+00</point>
      <point r="2.380"> -0.35661029E+00  0.37430228E+00  0.17816733E+00 -0.19096283E+00</point>
      <point r="2.400"> -0.35327026E+00  0.37269397E+00  0.18028538E+00 -0.18786274E+00</point>
      <point r="2.420"> -0.34989368E+00  0.37100009E+00  0.18226645E+00 -0.18479699E+00</point>
      <point r="2.440"> -0.34648325E+00  0.36922380E+00  0.18411496E+00 -0.18176586E+00</point>
      <point r="2.460"> -0.34304166E+00  0.36736824E+00  0.18583521E+00 -0.17876956E+00</point>
      <point r="2.480"> -0.33957149E+00  0.36543649E+00  0.18743136E+00 -0.17580828E+00</point>
      <point r="2.500"> -0.33607527E+00  0.36343159E+00  0.18890747E+00 -0.17288215E+00</point>
      <point r="2.520"> -0.33255549E+00  0.36135653E+00  0.19026744E+00 -0.16999129E+00</point>
      <point r="2.540"> -0.32901452E+00  0.35921422E+00  0.19151508E+00 -0.16713576E+00</point>
      <point r="2.560"> -0.32545471E+00  0.35700753E+00  0.19265407E+00 -0.16431560E+00</point>
      <point r="2.580"> -0.32187831E+00  0.35473927E+00  0.19368799E+00 -0.16153083E+00</point>
      <point r="2.600"> -0.31828752E+00  0.35241219E+00  0.19462028E+00 -0.15878143E+00</point>
      <point r="2.620"> -0.31468445E+00  0.35002896E+00  0.19545432E+00 -0.15606736E+00</point>
      <point r="2.640"> -0.31107119E+00  0.34759223E+00  0.19619337E+00 -0.15338853E+00</point>
      <point r="2.660"> -0.30744971E+00  0.34510455E+00  0.19684059E+00 -0.15074487E+00</point>
      <point r="2.680"> -0.30382195E+00  0.34256844E+00  0.19739906E+00 -0.14813626E+00</point>
      <point r="2.700"> -0.30018978E+00  0.33998635E+00  0.19787180E+00 -0.14556256E+00</point>
      <point r="2.720"> -0.29655501E+00  0.33736071E+00  0.19826171E+00 -0.14302363E+00</point>
      <point r="2.740"> -0.29291941E+00  0.33469386E+00  0.19857164E+00 -0.14051930E+00</point>
      <point r="2.760"> -0.28928466E+00  0.33198812E+00  0.19880436E+00 -0.13804938E+00</point>
      <point r="2.780"> -0.28565242E+00  0.32924576E+00  0.19896258E+00 -0.13561368E+00</point>
      <point r="2.800"> -0.28202427E+00  0.32646899E+00  0.19904893E+00 -0.13321199E+00</point>
      <point r="2.820"> -0.27840176E+00  0.32365999E+00  0.19906597E+00 -0.13084408E+00</point>
      <point r="2.840"> -0.27478637E+00  0.32082089E+00  0.19901621E+00 -0.12850972E+00</point>
      <point r="2.860"> -0.27117954E+00  0.31795379E+00  0.19890209E+00 -0.12620867E+00</point>
      <point r="2.880"> -0.26758266E+00  0.31506070E+00  0.19872597E+00 -0.12394067E+00</point>
      <point r="2.900"> -0.26399705E+00  0.31214362E+00  0.19849016E+00 -0.12170546E+00</point>
      <point r="2.920"> -0.26042400E+00  0.30920449E+00  0.19819688E+00 -0.11950277E+00</point>
      <point r="2.940"> -0.25686474E+00  0.30624518E+00  0.19784830E+00 -0.11733232E+00</point>
      <point r="2.960"> -0.25332044E+00  0.30326752E+00  0.19744652E+00 -0.11519383E+00</point>
      <point r="2.980"> -0.24979223E+00  0.30027327E+00  0.19699354E+00 -0.11308700E+00</point>
      <point r="3.000"> -0.24628116E+00  0.29726415E+00  0.19649132E+00 -0.11101155E+00</point>
      <point r="3.020"> -0.24278827E+00  0.29424181E+00  0.19594173E+00 -0.10896716E+00</point>
      <point r="3.040"> -0.23931451E+00  0.29120785E+00  0.19534657E+00 -0.10695353E+00</point>
      <point r="3.060"> -0.23586081E+00  0.28816380E+00  0.19470755E+00 -0.10497036E+00</point>
      <point r="3.080"> -0.23242804E+00  0.28511116E+00  0.19402636E+00 -0.10301734E+00</point>
      <point r="3.100"> -0.22901701E+00  0.28205134E+00  0.19330456E+00 -0.10109413E+00</point>
      <point r="3.120"> -0.22562851E+00  0.27898573E+00  0.19254369E+00 -0.99200444E-01</point>
      <point r="3.140"> -0.22226327E+00  0.27591566E+00  0.19174520E+00 -0.97335946E-01</point>
      <point r="3.160"> -0.21892198E+00  0.27284240E+00  0.19091050E+00 -0.95500322E-01</point>
      <point r="3.180"> -0.21560530E+00  0.26976718E+00  0.19004092E+00 -0.93693253E-01</point>
      <point r="3.200"> -0.21231383E+00  0.26669121E+00  0.18913776E+00 -0.91914418E-01</point>
      <point r="3.220"> -0.20904817E+00  0.26361563E+00  0.18820226E+00 -0.90163498E-01</point>
      <point r="3.240"> -0.20580884E+00  0.26054154E+00  0.18723560E+00 -0.88440173E-01</point>
      <point r="3.260"> -0.20259637E+00  0.25747001E+00  0.18623895E+00 -0.86744125E-01</point>
      <point r="3.280"> -0.19941122E+00  0.25440208E+00  0.18521340E+00 -0.85075037E-01</point>
      <point r="3.300"> -0.19625385E+00  0.25133875E+00  0.18416004E+00 -0.83432590E-01</point>
      <point r="3.320"> -0.19312466E+00  0.24828097E+00  0.18307991E+00 -0.81816470E-01</point>
      <point r="3.340"> -0.19002361E+00  0.24522865E+00  0.18196839E+00 -0.80226083E-01</point>
      <point r="3.360"> -0.18695138E+00  0.24218360E+00  0.18083508E+00 -0.78661691E-01</point>
      <point r="3.380"> -0.18390846E+00  0.23914685E+00  0.17967818E+00 -0.77122688E-01</point>
      <point r="3.400"> -0.18089516E+00  0.23611926E+00  0.17849866E+00 -0.75608762E-01</point>
      <point r="3.420"> -0.17791177E+00  0.23310164E+00  0.17729747E+00 -0.74119601E-01</point>
      <point r="3.440"> -0.17495853E+00  0.23009477E+00  0.17607552E+00 -0.72654898E-01</point>
      <point r="3.460"> -0.17203569E+00  0.22709940E+00  0.17483370E+00 -0.71214346E-01</point>
      <point r="3.480"> -0.16914345E+00  0.22411626E+00  0.17357287E+00 -0.69797640E-01</point>
      <point r="3.500"> -0.16628202E+00  0.22114605E+00  0.17229389E+00 -0.68404477E-01</point>
      <point r="3.520"> -0.16345155E+00  0.21818945E+00  0.17099757E+00 -0.67034556E-01</point>
      <point r="3.540"> -0.16065221E+00  0.21524709E+00  0.16968471E+00 -0.65687580E-01</point>
      <point r="3.560"> -0.15788411E+00  0.21231960E+00  0.16835610E+00 -0.64363252E-01</point>
      <point r="3.580"> -0.15514738E+00  0.20940758E+00  0.16701250E+00 -0.63061277E-01</point>
      <point r="3.600"> -0.15244211E+00  0.20651160E+00  0.16565467E+00 -0.61781365E-01</point>
      <point r="3.620"> -0.14976837E+00  0.20363222E+00  0.16428333E+00 -0.60523225E-01</point>
      <point r="3.640"> -0.14712623E+00  0.20076997E+00  0.16289921E+00 -0.59286571E-01</point>
      <point r="3.660"> -0.14451572E+00  0.19792535E+00  0.16150300E+00 -0.58071119E-01</point>
      <point r="3.680"> -0.14193689E+00  0.19509886E+00  0.16009539E+00 -0.56876584E-01</point>
      <point r="3.700"> -0.13938974E+00  0.19229095E+00  0.15867706E+00 -0.55702689E-01</point>
      <point r="3.720"> -0.13687428E+00  0.18950208E+00  0.15724865E+00 -0.54549155E-01</point>
      <point r="3.740"> -0.13439048E+00  0.18673267E+00  0.15581082E+00 -0.53415707E-01</point>
      <point r="3.760"> -0.13193832E+00  0.18398313E+00  0.15436420E+00 -0.52302072E-01</point>
      <point r="3.780"> -0.12951775E+00  0.18125384E+00  0.15290940E+00 -0.51207981E-01</point>
      <point r="3.800"> -0.12712873E+00  0.17854517E+00  0.15144703E+00 -0.50133165E-01</point>
      <point r="3.820"> -0.12477118E+00  0.17585747E+00  0.14997768E+00 -0.49077360E-01</point>
      <point r="3.840"> -0.12244502E+00  0.17319108E+00  0.14850194E+00 -0.48040301E-01</point>
      <point r="3.860"> -0.12015016E+00  0.17054630E+00  0.14702036E+00 -0.47021729E-01</point>
      <point r="3.880"> -0.11788651E+00  0.16792343E+00  0.14553351E+00 -0.46021386E-01</point>
      <point r="3.900"> -0.11565394E+00  0.16532277E+00  0.14404193E+00 -0.45039016E-01</point>
      <point r="3.920"> -0.11345233E+00  0.16274456E+00  0.14254616E+00 -0.44074366E-01</point>
      <point r="3.940"> -0.11128157E+00  0.16018906E+00  0.14104673E+00 -0.43127185E-01</point>
      <point r="3.960"> -0.10914149E+00  0.15765651E+00  0.13954414E+00 -0.42197225E-01</point>
      <point r="3.980"> -0.10703196E+00  0.15514712E+00  0.13803890E+00 -0.41284241E-01</point>
      <point r="4.000"> -0.10495281E+00  0.15266111E+00  0.13653152E+00 -0.40387988E-01</point>
      <point r="4.020"> -0.10290388E+00  0.15019866E+00  0.13502248E+00 -0.39508226E-01</point>
      <point r="4.040"> -0.10088500E+00  0.14775996E+00  0.13351225E+00 -0.38644715E-01</point>
      <point r="4.060"> -0.98895994E-01  0.14534518E+00  0.13200131E+00 -0.37797221E-01</point>
      <point r="4.080"> -0.96936668E-01  0.14295445E+00  0.13049012E+00 -0.36965508E-01</point>
      <point r="4.100"> -0.95006832E-01  0.14058794E+00  0.12897913E+00 -0.36149345E-01</point>
      <point r="4.120"> -0.93106288E-01  0.13824576E+00  0.12746878E+00 -0.35348503E-01</point>
      <point r="4.140"> -0.91234832E-01  0.13592804E+00  0.12595951E+00 -0.34562755E-01</point>
      <point r="4.160"> -0.89392254E-01  0.13363488E+00  0.12445176E+00 -0.33791876E-01</point>
      <point r="4.180"> -0.87578341E-01  0.13136637E+00  0.12294593E+00 -0.33035644E-01</point>
      <point r="4.200"> -0.85792872E-01  0.12912261E+00  0.12144244E+00 -0.32293838E-01</point>
      <point r="4.220"> -0.84035622E-01  0.12690366E+00  0.11994169E+00 -0.31566241E-01</point>
      <point r="4.240"> -0.82306363E-01  0.12470959E+00  0.11844408E+00 -0.30852637E-01</point>
      <point r="4.260"> -0.80604862E-01  0.12254044E+00  0.11694998E+00 -0.30152813E-01</point>
      <point r="4.280"> -0.78930881E-01  0.12039626E+00  0.11545978E+00 -0.29466556E-01</point>
      <point r="4.300"> -0.77284180E-01  0.11827708E+00  0.11397384E+00 -0.28793659E-01</point>
      <point r="4.320"> -0.75664512E-01  0.11618291E+00  0.11249252E+00 -0.28133915E-01</point>
      <point r="4.340"> -0.74071630E-01  0.11411379E+00  0.11101617E+00 -0.27487118E-01</point>
      <point r="4.360"> -0.72505281E-01  0.11206969E+00  0.10954513E+00 -0.26853067E-01</point>
      <point r="4.380"> -0.70965211E-01  0.11005062E+00  0.10807974E+00 -0.26231560E-01</point>
      <point r="4.400"> -0.69451163E-01  0.10805656E+00  0.10662032E+00 -0.25622400E-01</point>
      <point r="4.420"> -0.67962875E-01  0.10608748E+00  0.10516717E+00 -0.25025391E-01</point>
      <point r="4.440"> -0.66500085E-01  0.10414336E+00  0.10372062E+00 -0.24440339E-01</point>
      <point r="4.460"> -0.65062527E-01  0.10222414E+00  0.10228095E+00 -0.23867052E-01</point>
      <point r="4.480"> -0.63649934E-01  0.10032978E+00  0.10084846E+00 -0.23305341E-01</point>
      <point r="4.500"> -0.62262036E-01  0.98460220E-01  0.99423431E-01 -0.22755017E-01</point>
      <point r="4.520"> -0.60898563E-01  0.96615394E-01  0.98006128E-01 -0.22215895E-01</point>
      <point r="4.540"> -0.59559242E-01  0.94795229E-01  0.96596820E-01 -0.21687792E-01</point>
      <point r="4.560"> -0.58243799E-01  0.92999646E-01  0.95195764E-01 -0.21170527E-01</point>
      <point r="4.580"> -0.56951959E-01  0.91228556E-01  0.93803208E-01 -0.20663920E-01</point>
      <point r="4.600"> -0.55683446E-01  0.89481866E-01  0.92419393E-01 -0.20167794E-01</point>
      <point r="4.620"> -0.54437983E-01  0.87759476E-01  0.91044553E-01 -0.19681973E-01</point>
      <point r="4.640"> -0.53215293E-01  0.86061280E-01  0.89678913E-01 -0.19206285E-01</point>
      <point r="4.660"> -0.52015096E-01  0.84387166E-01  0.88322691E-01 -0.18740559E-01</point>
      <point r="4.680"> -0.50837116E-01  0.82737016E-01  0.86976098E-01 -0.18284624E-01</point>
      <point r="4.700"> -0.49681073E-01  0.81110707E-01  0.85639337E-01 -0.17838315E-01</point>
      <point r="4.720"> -0.48546689E-01  0.79508110E-01  0.84312604E-01 -0.17401465E-01</point>
      <point r="4.740"> -0.47433684E-01  0.77929091E-01  0.82996087E-01 -0.16973913E-01</point>
      <point r="4.760"> -0.46341781E-01  0.76373511E-01  0.81689967E-01 -0.16555495E-01</point>
      <point r="4.780"> -0.45270702E-01  0.74841228E-01  0.80394419E-01 -0.16146054E-01</point>
      <point r="4.800"> -0.44220168E-01  0.73332093E-01  0.79109609E-01 -0.15745431E-01</point>
      <point r="4.820"> -0.43189903E-01  0.71845953E-01  0.77835696E-01 -0.15353473E-01</point>
      <point r="4.840"> -0.42179629E-01  0.70382653E-01  0.76572834E-01 -0.14970024E-01</point>
      <point r="4.860"> -0.41189071E-01  0.68942030E-01  0.75321167E-01 -0.14594933E-01</point>
      <point r="4.880"> -0.40217953E-01  0.67523921E-01  0.74080834E-01 -0.14228052E-01</point>
      <point r="4.900"> -0.39266002E-01  0.66128157E-01  0.72851967E-01 -0.13869232E-01</point>
      <point r="4.920"> -0.38332943E-01  0.64754565E-01  0.71634689E-01 -0.13518327E-01</point>
      <point r="4.940"> -0.37418505E-01  0.63402971E-01  0.70429117E-01 -0.13175194E-01</point>
      <point r="4.960"> -0.36522417E-01  0.62073195E-01  0.69235363E-01 -0.12839691E-01</point>
      <point r="4.980"> -0.35644409E-01  0.60765055E-01  0.68053529E-01 -0.12511677E-01</point>
      <point r="5.000"> -0.34784211E-01  0.59478365E-01  0.66883713E-01 -0.12191013E-01</point>
      <point r="5.020"> -0.33941557E-01  0.58212938E-01  0.65726004E-01 -0.11877565E-01</point>
      <point r="5.040"> -0.33116181E-01  0.56968582E-01  0.64580485E-01 -0.11571196E-01</point>
      <point r="5.060"> -0.32307818E-01  0.55745105E-01  0.63447234E-01 -0.11271774E-01</point>
      <point r="5.080"> -0.31516205E-01  0.54542308E-01  0.62326319E-01 -0.10979169E-01</point>
      <point r="5.100"> -0.30741081E-01  0.53359996E-01  0.61217806E-01 -0.10693250E-01</point>
      <point r="5.120"> -0.29982187E-01  0.52197965E-01  0.60121751E-01 -0.10413890E-01</point>
      <point r="5.140"> -0.29239265E-01  0.51056014E-01  0.59038205E-01 -0.10140965E-01</point>
      <point r="5.160"> -0.28512058E-01  0.49933938E-01  0.57967212E-01 -0.98743491E-02</point>
      <point r="5.180"> -0.27800313E-01  0.48831530E-01  0.56908813E-01 -0.96139214E-02</point>
      <point r="5.200"> -0.27103777E-01  0.47748581E-01  0.55863039E-01 -0.93595611E-02</point>
      <point r="5.220"> -0.26422199E-01  0.46684883E-01  0.54829917E-01 -0.91111498E-02</point>
      <point r="5.240"> -0.25755332E-01  0.45640223E-01  0.53809469E-01 -0.88685706E-02</point>
      <point r="5.260"> -0.25102929E-01  0.44614388E-01  0.52801711E-01 -0.86317082E-02</point>
      <point r="5.280"> -0.24464746E-01  0.43607166E-01  0.51806651E-01 -0.84004493E-02</point>
      <point r="5.300"> -0.23840541E-01  0.42618341E-01  0.50824296E-01 -0.81746820E-02</point>
      <point r="5.320"> -0.23230073E-01  0.41647697E-01  0.49854644E-01 -0.79542963E-02</point>
      <point r="5.340"> -0.22633106E-01  0.40695019E-01  0.48897690E-01 -0.77391837E-02</point>
      <point r="5.360"> -0.22049404E-01  0.39760087E-01  0.47953423E-01 -0.75292376E-02</point>
      <point r="5.380"> -0.21478734E-01  0.38842685E-01  0.47021828E-01 -0.73243527E-02</point>
      <point r="5.400"> -0.20920865E-01  0.37942594E-01  0.46102885E-01 -0.71244256E-02</point>
      <point r="5.420"> -0.20375569E-01  0.37059596E-01  0.45196567E-01 -0.69293545E-02</point>
      <point r="5.440"> -0.19842620E-01  0.36193471E-01  0.44302847E-01 -0.67390390E-02</point>
      <point r="5.460"> -0.19321793E-01  0.35344001E-01  0.43421690E-01 -0.65533804E-02</point>
      <point r="5.480"> -0.18812869E-01  0.34510967E-01  0.42553058E-01 -0.63722818E-02</point>
      <point r="5.500"> -0.18315628E-01  0.33694148E-01  0.41696909E-01 -0.61956475E-02</point>
      <point r="5.520"> -0.17829854E-01  0.32893328E-01  0.40853197E-01 -0.60233835E-02</point>
      <point r="5.540"> -0.17355333E-01  0.32108286E-01  0.40021872E-01 -0.58553974E-02</point>
      <point r="5.560"> -0.16891854E-01  0.31338804E-01  0.39202880E-01 -0.56915983E-02</point>
      <point r="5.580"> -0.16439209E-01  0.30584666E-01  0.38396163E-01 -0.55318968E-02</point>
      <point r="5.600"> -0.15997192E-01  0.29845653E-01  0.37601662E-01 -0.53762049E-02</point>
      <point r="5.620"> -0.15565598E-01  0.29121548E-01  0.36819310E-01 -0.52244362E-02</point>
      <point r="5.640"> -0.15144227E-01  0.28412136E-01  0.36049041E-01 -0.50765057E-02</point>
      <point r="5.660"> -0.14732880E-01  0.27717202E-01  0.35290784E-01 -0.49323299E-02</point>
      <point r="5.680"> -0.14331362E-01  0.27036530E-01  0.34544466E-01 -0.47918267E-02</point>
      <point r="5.700"> -0.13939478E-01  0.26369907E-01  0.33810008E-01 -0.46549154E-02</point>
      <point r="5.720"> -0.13557040E-01  0.25717121E-01  0.33087333E-01 -0.45215167E-02</point>
      <point r="5.740"> -0.13183857E-01  0.25077960E-01  0.32376357E-01 -0.43915528E-02</point>
      <point r="5.760"> -0.12819746E-01  0.24452214E-01  0.31676995E-01 -0.42649472E-02</point>
      <point r="5.780"> -0.12464522E-01  0.23839672E-01  0.30989161E-01 -0.41416247E-02</point>
      <point r="5.800"> -0.12118005E-01  0.23240128E-01  0.30312765E-01 -0.40215115E-02</point>
      <point r="5.820"> -0.11780019E-01  0.22653374E-01  0.29647714E-01 -0.39045353E-02</point>
      <point r="5.840"> -0.11450386E-01  0.22079205E-01  0.28993914E-01 -0.37906248E-02</point>
      <point r="5.860"> -0.11128936E-01  0.21517417E-01  0.28351270E-01 -0.36797101E-02</point>
      <point r="5.880"> -0.10815496E-01  0.20967807E-01  0.27719683E-01 -0.35717229E-02</point>
      <point r="5.900"> -0.10509901E-01  0.20430173E-01  0.27099053E-01 -0.34665956E-02</point>
      <point r="5.920"> -0.10211985E-01  0.19904317E-01  0.26489278E-01 -0.33642624E-02</point>
      <point r="5.940"> -0.99215860E-02  0.19390040E-01  0.25890255E-01 -0.32646584E-02</point>
      <point r="5.960"> -0.96385431E-02  0.18887145E-01  0.25301878E-01 -0.31677201E-02</point>
      <point r="5.980"> -0.93626994E-02  0.18395438E-01  0.24724042E-01 -0.30733851E-02</point>
      <point r="6.000"> -0.90938999E-02  0.17914726E-01  0.24156638E-01 -0.29815921E-02</point>
      <point r="6.020"> -0.88319920E-02  0.17444817E-01  0.23599558E-01 -0.28922814E-02</point>
      <point r="6.040"> -0.85768258E-02  0.16985521E-01  0.23052690E-01 -0.28053939E-02</point>
      <point r="6.060"> -0.83282537E-02  0.16536651E-01  0.22515923E-01 -0.27208720E-02</point>
      <point r="6.080"> -0.80861307E-02  0.16098021E-01  0.21989145E-01 -0.26386592E-02</point>
      <point r="6.100"> -0.78503139E-02  0.15669446E-01  0.21472242E-01 -0.25587000E-02</point>
      <point r="6.120"> -0.76206631E-02  0.15250744E-01  0.20965101E-01 -0.24809402E-02</point>
      <point r="6.140"> -0.73970404E-02  0.14841735E-01  0.20467605E-01 -0.24053264E-02</point>
      <point r="6.160"> -0.71793102E-02  0.14442239E-01  0.19979639E-01 -0.23318064E-02</point>
      <point r="6.180"> -0.69673394E-02  0.14052082E-01  0.19501087E-01 -0.22603292E-02</point>
      <point r="6.200"> -0.67609970E-02  0.13671087E-01  0.19031832E-01 -0.21908447E-02</point>
      <point r="6.220"> -0.65601545E-02  0.13299082E-01  0.18571755E-01 -0.21233038E-02</point>
      <point r="6.240"> -0.63646856E-02  0.12935896E-01  0.18120741E-01 -0.20576584E-02</point>
      <point r="6.260"> -0.61744662E-02  0.12581362E-01  0.17678669E-01 -0.19938616E-02</point>
      <point r="6.280"> -0.59893747E-02  0.12235312E-01  0.17245422E-01 -0.19318672E-02</point>
      <point r="6.300"> -0.58092914E-02  0.11897582E-01  0.16820881E-01 -0.18716301E-02</point>
      <point r="6.320"> -0.56340991E-02  0.11568009E-01  0.16404928E-01 -0.18131061E-02</point>
      <point r="6.340"> -0.54636824E-02  0.11246432E-01  0.15997443E-01 -0.17562521E-02</point>
      <point r="6.360"> -0.52979285E-02  0.10932694E-01  0.15598308E-01 -0.17010258E-02</point>
      <point r="6.380"> -0.51367264E-02  0.10626638E-01  0.15207404E-01 -0.16473856E-02</point>
      <point r="6.400"> -0.49799674E-02  0.10328110E-01  0.14824613E-01 -0.15952912E-02</point>
      <point r="6.420"> -0.48275446E-02  0.10036957E-01  0.14449815E-01 -0.15447029E-02</point>
      <point r="6.440"> -0.46793536E-02  0.97530295E-02  0.14082893E-01 -0.14955820E-02</point>
      <point r="6.460"> -0.45352917E-02  0.94761795E-02  0.13723730E-01 -0.14478904E-02</point>
      <point r="6.480"> -0.43952583E-02  0.92062609E-02  0.13372206E-01 -0.14015912E-02</point>
      <point r="6.500"> -0.42591548E-02  0.89431299E-02  0.13028206E-01 -0.13566481E-02</point>
      <point r="6.520"> -0.41268846E-02  0.86866449E-02  0.12691611E-01 -0.13130256E-02</point>
      <point r="6.540"> -0.39983529E-02  0.84366660E-02  0.12362308E-01 -0.12706890E-02</point>
      <point r="6.560"> -0.38734670E-02  0.81930557E-02  0.12040178E-01 -0.12296045E-02</point>
      <point r="6.580"> -0.37521360E-02  0.79556786E-02  0.11725108E-01 -0.11897390E-02</point>
      <point r="6.600"> -0.36342707E-02  0.77244011E-02  0.11416982E-01 -0.11510601E-02</point>
      <point r="6.620"> -0.35197841E-02  0.74990919E-02  0.11115688E-01 -0.11135361E-02</point>
      <point r="6.640"> -0.34085907E-02  0.72796216E-02  0.10821111E-01 -0.10771362E-02</point>
      <point r="6.660"> -0.33006068E-02  0.70658630E-02  0.10533139E-01 -0.10418302E-02</point>
      <point r="6.680"> -0.31933020E-02  0.68517920E-02  0.10240680E-01 -0.10071522E-02</point>
      <point r="6.700"> -0.30915104E-02  0.66491081E-02  0.99656210E-02 -0.97395531E-03</point>
      <point r="6.720"> -0.29926912E-02  0.64517743E-02  0.96968497E-02 -0.94176639E-03</point>
      <point r="6.740"> -0.28967678E-02  0.62596717E-02  0.94342574E-02 -0.91055799E-03</point>
      <point r="6.760"> -0.28036651E-02  0.60726828E-02  0.91777358E-02 -0.88030326E-03</point>
      <point r="6.780"> -0.27133096E-02  0.58906927E-02  0.89271777E-02 -0.85097598E-03</point>
      <point r="6.800"> -0.26256294E-02  0.57135878E-02  0.86824766E-02 -0.82255057E-03</point>
      <point r="6.820"> -0.25405543E-02  0.55412570E-02  0.84435271E-02 -0.79500207E-03</point>
      <point r="6.840"> -0.24580154E-02  0.53735908E-02  0.82102246E-02 -0.76830608E-03</point>
      <point r="6.860"> -0.23779456E-02  0.52104817E-02  0.79824655E-02 -0.74243883E-03</point>
      <point r="6.880"> -0.23002791E-02  0.50518239E-02  0.77601473E-02 -0.71737711E-03</point>
      <point r="6.900"> -0.22249516E-02  0.48975138E-02  0.75431683E-02 -0.69309828E-03</point>
      <point r="6.920"> -0.21519003E-02  0.47474494E-02  0.73314281E-02 -0.66958023E-03</point>
      <point r="6.940"> -0.20810639E-02  0.46015306E-02  0.71248270E-02 -0.64680144E-03</point>
      <point r="6.960"> -0.20123823E-02  0.44596592E-02  0.69232666E-02 -0.62474087E-03</point>
      <point r="6.980"> -0.19457969E-02  0.43217387E-02  0.67266495E-02 -0.60337804E-03</point>
      <point r="7.000"> -0.18812505E-02  0.41876745E-02  0.65348795E-02 -0.58269297E-03</point>
      <point r="7.020"> -0.18186870E-02  0.40573735E-02  0.63478612E-02 -0.56266617E-03</point>
      <point r="7.040"> -0.17580519E-02  0.39307447E-02  0.61655007E-02 -0.54327865E-03</point>
      <point r="7.060"> -0.16992918E-02  0.38076986E-02  0.59877048E-02 -0.52451189E-03</point>
      <point r="7.080"> -0.16423546E-02  0.36881475E-02  0.58143819E-02 -0.50634787E-03</point>
      <point r="7.100"> -0.15871894E-02  0.35720053E-02  0.56454410E-02 -0.48876898E-03</point>
      <point r="7.120"> -0.15337465E-02  0.34591878E-02  0.54807927E-02 -0.47175811E-03</point>
      <point r="7.140"> -0.14819776E-02  0.33496121E-02  0.53203486E-02 -0.45529857E-03</point>
      <point r="7.160"> -0.14318352E-02  0.32431972E-02  0.51640213E-02 -0.43937408E-03</point>
      <point r="7.180"> -0.13832733E-02  0.31398637E-02  0.50117248E-02 -0.42396883E-03</point>
      <point r="7.200"> -0.13362468E-02  0.30395336E-02  0.48633741E-02 -0.40906737E-03</point>
      <point r="7.220"> -0.12907117E-02  0.29421308E-02  0.47188856E-02 -0.39465470E-03</point>
      <point r="7.240"> -0.12466253E-02  0.28475804E-02  0.45781765E-02 -0.38071620E-03</point>
      <point r="7.260"> -0.12039457E-02  0.27558094E-02  0.44411656E-02 -0.36723761E-03</point>
      <point r="7.280"> -0.11626322E-02  0.26667459E-02  0.43077727E-02 -0.35420509E-03</point>
      <point r="7.300"> -0.11226450E-02  0.25803198E-02  0.41779186E-02 -0.34160515E-03</point>
      <point r="7.320"> -0.10839454E-02  0.24964624E-02  0.40515256E-02 -0.32942465E-03</point>
      <point r="7.340"> -0.10464958E-02  0.24151065E-02  0.39285170E-02 -0.31765083E-03</point>
      <point r="7.360"> -0.10102592E-02  0.23361863E-02  0.38088174E-02 -0.30627125E-03</point>
      <point r="7.380"> -0.97519987E-03  0.22596372E-02  0.36923524E-02 -0.29527382E-03</point>
      <point r="7.400"> -0.94128284E-03  0.21853965E-02  0.35790490E-02 -0.28464679E-03</point>
      <point r="7.420"> -0.90847410E-03  0.21134024E-02  0.34688352E-02 -0.27437870E-03</point>
      <point r="7.440"> -0.87674050E-03  0.20435946E-02  0.33616404E-02 -0.26445844E-03</point>
      <point r="7.460"> -0.84604977E-03  0.19759142E-02  0.32573949E-02 -0.25487518E-03</point>
      <point r="7.480"> -0.81637049E-03  0.19103037E-02  0.31560304E-02 -0.24561842E-03</point>
      <point r="7.500"> -0.78767205E-03  0.18467067E-02  0.30574796E-02 -0.23667792E-03</point>
      <point r="7.520"> -0.75992468E-03  0.17850682E-02  0.29616765E-02 -0.22804375E-03</point>
      <point r="7.540"> -0.73309937E-03  0.17253343E-02  0.28685562E-02 -0.21970624E-03</point>
      <point r="7.560"> -0.70716793E-03  0.16674526E-02  0.27780550E-02 -0.21165602E-03</point>
      <point r="7.580"> -0.68210290E-03  0.16113717E-02  0.26901101E-02 -0.20388395E-03</point>
      <point r="7.600"> -0.65787757E-03  0.15570415E-02  0.26046603E-02 -0.19638118E-03</point>
      <point r="7.620"> -0.63446597E-03  0.15044130E-02  0.25216450E-02 -0.18913910E-03</point>
      <point r="7.640"> -0.61184283E-03  0.14534385E-02  0.24410051E-02 -0.18214935E-03</point>
      <point r="7.660"> -0.58998357E-03  0.14040712E-02  0.23626826E-02 -0.17540381E-03</point>
      <point r="7.680"> -0.56886432E-03  0.13562657E-02  0.22866204E-02 -0.16889459E-03</point>
      <point r="7.700"> -0.54846183E-03  0.13099775E-02  0.22127626E-02 -0.16261403E-03</point>
      <point r="7.720"> -0.52875354E-03  0.12651632E-02  0.21410544E-02 -0.15655471E-03</point>
      <point r="7.740"> -0.50971751E-03  0.12217806E-02  0.20714420E-02 -0.15070940E-03</point>
      <point r="7.760"> -0.49133240E-03  0.11797884E-02  0.20038728E-02 -0.14507111E-03</point>
      <point r="7.780"> -0.47357751E-03  0.11391464E-02  0.19382952E-02 -0.13963302E-03</point>
      <point r="7.800"> -0.45643270E-03  0.10998153E-02  0.18746586E-02 -0.13438856E-03</point>
      <point r="7.820"> -0.43987843E-03  0.10617569E-02  0.18129135E-02 -0.12933132E-03</point>
      <point r="7.840"> -0.42389570E-03  0.10249339E-02  0.17530113E-02 -0.12445508E-03</point>
      <point r="7.860"> -0.40846608E-03  0.98930991E-03  0.16949044E-02 -0.11975384E-03</point>
      <point r="7.880"> -0.39357167E-03  0.95484959E-03  0.16385465E-02 -0.11522174E-03</point>
      <point r="7.900"> -0.37919507E-03  0.92151840E-03  0.15838919E-02 -0.11085312E-03</point>
      <point r="7.920"> -0.36531942E-03  0.88928273E-03  0.15308961E-02 -0.10664248E-03</point>
      <point r="7.940"> -0.35192835E-03  0.85810984E-03  0.14795154E-02 -0.10258450E-03</point>
      <point r="7.960"> -0.33900595E-03  0.82796784E-03  0.14297071E-02 -0.98674024E-04</point>
      <point r="7.980"> -0.32653681E-03  0.79882566E-03  0.13814296E-02 -0.94906034E-04</point>
      <point r="8.000"> -0.31450597E-03  0.77065307E-03  0.13346419E-02 -0.91275680E-04</point>
      <point r="8.020"> -0.30289891E-03  0.74342063E-03  0.12893041E-02 -0.87778263E-04</point>
      <point r="8.040"> -0.29170155E-03  0.71709968E-03  0.12453772E-02 -0.84409224E-04</point>
      <point r="8.060"> -0.28090023E-03  0.69166233E-03  0.12028229E-02 -0.81164150E-04</point>
      <point r="8.080"> -0.27048173E-03  0.66708143E-03  0.11616039E-02 -0.78038763E-04</point>
      <point r="8.100"> -0.26043319E-03  0.64333057E-03  0.11216838E-02 -0.75028920E-04</point>
      <point r="8.120"> -0.25074217E-03  0.62038406E-03  0.10830269E-02 -0.72130609E-04</point>
      <point r="8.140"> -0.24139660E-03  0.59821688E-03  0.10455982E-02 -0.69339943E-04</point>
      <point r="8.160"> -0.23238479E-03  0.57680472E-03  0.10093639E-02 -0.66653160E-04</point>
      <point r="8.180"> -0.22369540E-03  0.55612392E-03  0.97429047E-03 -0.64066615E-04</point>
      <point r="8.200"> -0.21531744E-03  0.53615149E-03  0.94034557E-03 -0.61576782E-04</point>
      <point r="8.220"> -0.20724027E-03  0.51686507E-03  0.90749740E-03 -0.59180247E-04</point>
      <point r="8.240"> -0.19945359E-03  0.49824290E-03  0.87571498E-03 -0.56873704E-04</point>
      <point r="8.260"> -0.19194739E-03  0.48026384E-03  0.84496802E-03 -0.54653956E-04</point>
      <point r="8.280"> -0.18471200E-03  0.46290735E-03  0.81522697E-03 -0.52517908E-04</point>
      <point r="8.300"> -0.17773805E-03  0.44615345E-03  0.78646296E-03 -0.50462566E-04</point>
      <point r="8.320"> -0.17101647E-03  0.42998275E-03  0.75864784E-03 -0.48485035E-04</point>
      <point r="8.340"> -0.16453846E-03  0.41437638E-03  0.73175409E-03 -0.46582512E-04</point>
      <point r="8.360"> -0.15829553E-03  0.39931602E-03  0.70575490E-03 -0.44752288E-04</point>
      <point r="8.380"> -0.15227943E-03  0.38478387E-03  0.68062408E-03 -0.42991744E-04</point>
      <point r="8.400"> -0.14648219E-03  0.37076264E-03  0.65633608E-03 -0.41298346E-04</point>
      <point r="8.420"> -0.14089610E-03  0.35723555E-03  0.63286598E-03 -0.39669645E-04</point>
      <point r="8.440"> -0.13551369E-03  0.34418630E-03  0.61018946E-03 -0.38103274E-04</point>
      <point r="8.460"> -0.13032774E-03  0.33159904E-03  0.58828283E-03 -0.36596944E-04</point>
      <point r="8.480"> -0.12533126E-03  0.31945842E-03  0.56712295E-03 -0.35148445E-04</point>
      <point r="8.500"> -0.12051749E-03  0.30774952E-03  0.54668729E-03 -0.33755640E-04</point>
      <point r="8.520"> -0.11587990E-03  0.29645785E-03  0.52695386E-03 -0.32416465E-04</point>
      <point r="8.540"> -0.11141215E-03  0.28556937E-03  0.50790123E-03 -0.31128927E-04</point>
      <point r="8.560"> -0.10710814E-03  0.27507044E-03  0.48950854E-03 -0.29891099E-04</point>
      <point r="8.580"> -0.10296196E-03  0.26494785E-03  0.47175542E-03 -0.28701122E-04</point>
      <point r="8.600"> -0.98967883E-04  0.25518875E-03  0.45462206E-03 -0.27557200E-04</point>
      <point r="8.620"> -0.95120402E-04  0.24578072E-03  0.43808913E-03 -0.26457600E-04</point>
      <point r="8.640"> -0.91414175E-04  0.23671168E-03  0.42213782E-03 -0.25400648E-04</point>
      <point r="8.660"> -0.87844046E-04  0.22796993E-03  0.40674982E-03 -0.24384730E-04</point>
      <point r="8.680"> -0.84405034E-04  0.21954414E-03  0.39190727E-03 -0.23408287E-04</point>
      <point r="8.700"> -0.81092329E-04  0.21142332E-03  0.37759280E-03 -0.22469814E-04</point>
      <point r="8.720"> -0.77901281E-04  0.20359680E-03  0.36378951E-03 -0.21567862E-04</point>
      <point r="8.740"> -0.74827405E-04  0.19605428E-03  0.35048094E-03 -0.20701031E-04</point>
      <point r="8.760"> -0.71866369E-04  0.18878576E-03  0.33765108E-03 -0.19867971E-04</point>
      <point r="8.780"> -0.69013989E-04  0.18178154E-03  0.32528434E-03 -0.19067382E-04</point>
      <point r="8.800"> -0.66266229E-04  0.17503225E-03  0.31336557E-03 -0.18298008E-04</point>
      <point r="8.820"> -0.63619194E-04  0.16852882E-03  0.30188003E-03 -0.17558641E-04</point>
      <point r="8.840"> -0.61069125E-04  0.16226244E-03  0.29081339E-03 -0.16848114E-04</point>
      <point r="8.860"> -0.58612397E-04  0.15622461E-03  0.28015171E-03 -0.16165304E-04</point>
      <point r="8.880"> -0.56245513E-04  0.15040709E-03  0.26988145E-03 -0.15509129E-04</point>
      <point r="8.900"> -0.53965099E-04  0.14480192E-03  0.25998944E-03 -0.14878547E-04</point>
      <point r="8.920"> -0.51767904E-04  0.13940139E-03  0.25046290E-03 -0.14272553E-04</point>
      <point r="8.940"> -0.49650793E-04  0.13419804E-03  0.24128939E-03 -0.13690179E-04</point>
      <point r="8.960"> -0.47610744E-04  0.12918466E-03  0.23245686E-03 -0.13130495E-04</point>
      <point r="8.980"> -0.45644845E-04  0.12435428E-03  0.22395357E-03 -0.12592604E-04</point>
      <point r="9.000"> -0.43750290E-04  0.11970017E-03  0.21576816E-03 -0.12075643E-04</point>
      <point r="9.020"> -0.41924377E-04  0.11521581E-03  0.20788957E-03 -0.11578781E-04</point>
      <point r="9.040"> -0.40164501E-04  0.11089490E-03  0.20030708E-03 -0.11101220E-04</point>
      <point r="9.060"> -0.38468155E-04  0.10673138E-03  0.19301031E-03 -0.10642189E-04</point>
      <point r="9.080"> -0.36832926E-04  0.10271937E-03  0.18598914E-03 -0.10200949E-04</point>
      <point r="9.100"> -0.35256488E-04  0.98853194E-04  0.17923381E-03 -0.97767900E-05</point>
      <point r="9.120"> -0.33736605E-04  0.95127382E-04  0.17273481E-03 -0.93690264E-05</point>
      <point r="9.140"> -0.32271124E-04  0.91536649E-04  0.16648294E-03 -0.89770009E-05</point>
      <point r="9.160"> -0.30857974E-04  0.88075894E-04  0.16046930E-03 -0.86000813E-05</point>
      <point r="9.180"> -0.29495162E-04  0.84740195E-04  0.15468522E-03 -0.82376598E-05</point>
      <point r="9.200"> -0.28180771E-04  0.81524806E-04  0.14912235E-03 -0.78891525E-05</point>
      <point r="9.220"> -0.26912960E-04  0.78425146E-04  0.14377257E-03 -0.75539983E-05</point>
      <point r="9.240"> -0.25689955E-04  0.75436799E-04  0.13862801E-03 -0.72316582E-05</point>
      <point r="9.260"> -0.24510053E-04  0.72555510E-04  0.13368108E-03 -0.69216144E-05</point>
      <point r="9.280"> -0.23371617E-04  0.69777172E-04  0.12892441E-03 -0.66233700E-05</point>
      <point r="9.300"> -0.22273074E-04  0.67097833E-04  0.12435087E-03 -0.63364475E-05</point>
      <point r="9.320"> -0.21212911E-04  0.64513681E-04  0.11995357E-03 -0.60603888E-05</point>
      <point r="9.340"> -0.20189677E-04  0.62021047E-04  0.11572583E-03 -0.57947540E-05</point>
      <point r="9.360"> -0.19201976E-04  0.59616394E-04  0.11166120E-03 -0.55391211E-05</point>
      <point r="9.380"> -0.18248469E-04  0.57296320E-04  0.10775344E-03 -0.52930852E-05</point>
      <point r="9.400"> -0.17327870E-04  0.55057548E-04  0.10399651E-03 -0.50562577E-05</point>
      <point r="9.420"> -0.16438945E-04  0.52896923E-04  0.10038459E-03 -0.48282660E-05</point>
      <point r="9.440"> -0.15580508E-04  0.50811413E-04  0.96912033E-04 -0.46087527E-05</point>
      <point r="9.460"> -0.14751422E-04  0.48798097E-04  0.93573405E-04 -0.43973754E-05</point>
      <point r="9.480"> -0.13950596E-04  0.46854167E-04  0.90363444E-04 -0.41938054E-05</point>
      <point r="9.500"> -0.13176982E-04  0.44976923E-04  0.87277072E-04 -0.39977279E-05</point>
      <point r="9.520"> -0.12429578E-04  0.43163770E-04  0.84309386E-04 -0.38088413E-05</point>
      <point r="9.540"> -0.11707418E-04  0.41412212E-04  0.81455655E-04 -0.36268563E-05</point>
      <point r="9.560"> -0.11009581E-04  0.39719852E-04  0.78711312E-04 -0.34514959E-05</point>
      <point r="9.580"> -0.10335180E-04  0.38084386E-04  0.76071953E-04 -0.32824949E-05</point>
      <point r="9.600"> -0.96833666E-05  0.36503602E-04  0.73533328E-04 -0.31195991E-05</point>
      <point r="9.620"> -0.90533269E-05  0.34975375E-04  0.71091342E-04 -0.29625652E-05</point>
      <point r="9.640"> -0.84442811E-05  0.33497664E-04  0.68742044E-04 -0.28111601E-05</point>
      <point r="9.660"> -0.78554818E-05  0.32068511E-04  0.66481627E-04 -0.26651609E-05</point>
      <point r="9.680"> -0.72862129E-05  0.30686037E-04  0.64306425E-04 -0.25243540E-05</point>
      <point r="9.700"> -0.67357884E-05  0.29348439E-04  0.62212904E-04 -0.23885352E-05</point>
      <point r="9.720"> -0.62035512E-05  0.28053988E-04  0.60197662E-04 -0.22575090E-05</point>
      <point r="9.740"> -0.56888716E-05  0.26801024E-04  0.58257423E-04 -0.21310883E-05</point>
      <point r="9.760"> -0.51911468E-05  0.25587957E-04  0.56389034E-04 -0.20090942E-05</point>
      <point r="9.780"> -0.47097993E-05  0.24413263E-04  0.54589461E-04 -0.18913556E-05</point>
      <point r="9.800"> -0.42442763E-05  0.23275480E-04  0.52855786E-04 -0.17777088E-05</point>
      <point r="9.820"> -0.37940482E-05  0.22173209E-04  0.51185201E-04 -0.16679974E-05</point>
      <point r="9.840"> -0.33586082E-05  0.21105108E-04  0.49575007E-04 -0.15620716E-05</point>
      <point r="9.860"> -0.29374709E-05  0.20069894E-04  0.48022611E-04 -0.14597884E-05</point>
      <point r="9.880"> -0.25301716E-05  0.19066336E-04  0.46525520E-04 -0.13610111E-05</point>
      <point r="9.900"> -0.21362654E-05  0.18093257E-04  0.45081339E-04 -0.12656090E-05</point>
      <point r="9.920"> -0.17553265E-05  0.17149530E-04  0.43687770E-04 -0.11734572E-05</point>
      <point r="9.940"> -0.13869471E-05  0.16234078E-04  0.42342604E-04 -0.10844362E-05</point>
      <point r="9.960"> -0.10307370E-05  0.15345867E-04  0.41043724E-04 -0.99843198E-06</point>
      <point r="9.980"> -0.68632232E-06  0.14483913E-04  0.39789095E-04 -0.91533559E-06</point>
      <point r="10.000"> -0.35334542E-06  0.13647270E-04  0.38576770E-04 -0.83504286E-06</point>
      <point r="10.020"> -0.31463686E-07  0.12835038E-04  0.37404879E-04 -0.75745428E-06</point>
      <point r="10.040">  0.27965087E-06  0.12046353E-04  0.36271629E-04 -0.68247479E-06</point>
      <point r="10.060">  0.58031235E-06  0.11280392E-04  0.35175305E-04 -0.61001360E-06</point>
      <point r="10.080">  0.87082160E-06  0.10536367E-04  0.34114262E-04 -0.53998393E-06</point>
      <point r="10.100">  0.11514668E-05  0.98135256E-05  0.33086925E-04 -0.47230291E-06</point>
      <point r="10.120">  0.14225238E-05  0.91111494E-05  0.32091789E-04 -0.40689137E-06</point>
      <point r="10.140">  0.16842574E-05  0.84285521E-05  0.31127411E-04 -0.34367366E-06</point>
      <point r="10.160">  0.19369208E-05  0.77650784E-05  0.30192411E-04 -0.28257752E-06</point>
      <point r="10.180">  0.21807570E-05  0.71201024E-05  0.29285473E-04 -0.22353390E-06</point>
      <point r="10.200">  0.24159990E-05  0.64930271E-05  0.28405335E-04 -0.16647682E-06</point>
      <point r="10.220">  0.26428701E-05  0.58832824E-05  0.27550794E-04 -0.11134326E-06</point>
      <point r="10.240">  0.28615847E-05  0.52903244E-05  0.26720701E-04 -0.58072946E-07</point>
      <point r="10.260">  0.30723484E-05  0.47136339E-05  0.25913958E-04 -0.66082972E-08</point>
      <point r="10.280">  0.32753585E-05  0.41527158E-05  0.25129519E-04  0.43105738E-07</point>
      <point r="10.300">  0.34708046E-05  0.36070976E-05  0.24366385E-04  0.91121797E-07</point>
      <point r="10.320">  0.36588688E-05  0.30763286E-05  0.23623605E-04  0.13749021E-06</point>
      <point r="10.340">  0.38397261E-05  0.25599787E-05  0.22900272E-04  0.18225914E-06</point>
      <point r="10.360">  0.40135447E-05  0.20576381E-05  0.22195521E-04  0.22547463E-06</point>
      <point r="10.380">  0.41804866E-05  0.15689154E-05  0.21508532E-04  0.26718077E-06</point>
      <point r="10.400">  0.00000000E+00</point>
    </H_spline>
    <S_spline>
      <point r="0.000">  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.44000000E+01</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400">  0.79377352E+00 -0.26075042E+00  0.70399492E+00  0.71171912E+00</point>
      <point r="0.420">  0.79065131E+00 -0.26976300E+00  0.70164744E+00  0.71129666E+00</point>
      <point r="0.440">  0.78737534E+00 -0.27846013E+00  0.69897853E+00  0.71079864E+00</point>
      <point r="0.460">  0.78396009E+00 -0.28686030E+00  0.69597616E+00  0.71022068E+00</point>
      <point r="0.480">  0.78042009E+00 -0.29498185E+00  0.69262985E+00  0.70955866E+00</point>
      <point r="0.500">  0.77676972E+00 -0.30284277E+00  0.68893061E+00  0.70880871E+00</point>
      <point r="0.520">  0.77302306E+00 -0.31046056E+00  0.68487096E+00  0.70796726E+00</point>
      <point r="0.540">  0.76919373E+00 -0.31785209E+00  0.68044485E+00  0.70703102E+00</point>
      <point r="0.560">  0.76529482E+00 -0.32503350E+00  0.67564764E+00  0.70599696E+00</point>
      <point r="0.580">  0.76133874E+00 -0.33202014E+00  0.67047609E+00  0.70486235E+00</point>
      <point r="0.600">  0.75733721E+00 -0.33882649E+00  0.66492823E+00  0.70362472E+00</point>
      <point r="0.620">  0.75330118E+00 -0.34546610E+00  0.65900336E+00  0.70228190E+00</point>
      <point r="0.640">  0.74924081E+00 -0.35195162E+00  0.65270198E+00  0.70083196E+00</point>
      <point r="0.660">  0.74516544E+00 -0.35829475E+00  0.64602574E+00  0.69927327E+00</point>
      <point r="0.680">  0.74108360E+00 -0.36450623E+00  0.63897734E+00  0.69760442E+00</point>
      <point r="0.700">  0.73700299E+00 -0.37059589E+00  0.63156051E+00  0.69582428E+00</point>
      <point r="0.720">  0.73293054E+00 -0.37657261E+00  0.62377994E+00  0.69393194E+00</point>
      <point r="0.740">  0.72887235E+00 -0.38244440E+00  0.61564120E+00  0.69192676E+00</point>
      <point r="0.760">  0.72483380E+00 -0.38821839E+00  0.60715069E+00  0.68980829E+00</point>
      <point r="0.780">  0.72081950E+00 -0.39390087E+00  0.59831560E+00  0.68757632E+00</point>
      <point r="0.800">  0.71683340E+00 -0.39949734E+00  0.58914382E+00  0.68523085E+00</point>
      <point r="0.820">  0.71287877E+00 -0.40501253E+00  0.57964391E+00  0.68277207E+00</point>
      <point r="0.840">  0.70895824E+00 -0.41045045E+00  0.56982503E+00  0.68020037E+00</point>
      <point r="0.860">  0.70507387E+00 -0.41581443E+00  0.55969690E+00  0.67751632E+00</point>
      <point r="0.880">  0.70122717E+00 -0.42110715E+00  0.54926972E+00  0.67472067E+00</point>
      <point r="0.900">  0.69741913E+00 -0.42633071E+00  0.53855416E+00  0.67181434E+00</point>
      <point r="0.920">  0.69365028E+00 -0.43148665E+00  0.52756130E+00  0.66879838E+00</point>
      <point r="0.940">  0.68992072E+00 -0.43657598E+00  0.51630258E+00  0.66567402E+00</point>
      <point r="0.960">  0.68623013E+00 -0.44159925E+00  0.50478973E+00  0.66244262E+00</point>
      <point r="0.980">  0.68257788E+00 -0.44655657E+00  0.49303479E+00  0.65910567E+00</point>
      <point r="1.000">  0.67896296E+00 -0.45144765E+00  0.48105001E+00  0.65566478E+00</point>
      <point r="1.020">  0.67538413E+00 -0.45627185E+00  0.46884786E+00  0.65212168E+00</point>
      <point r="1.040">  0.67183985E+00 -0.46102819E+00  0.45644096E+00  0.64847820E+00</point>
      <point r="1.060">  0.66832838E+00 -0.46571540E+00  0.44384205E+00  0.64473629E+00</point>
      <point r="1.080">  0.66484777E+00 -0.47033196E+00  0.43106397E+00  0.64089798E+00</point>
      <point r="1.100">  0.66139591E+00 -0.47487613E+00  0.41811963E+00  0.63696538E+00</point>
      <point r="1.120">  0.65797056E+00 -0.47934596E+00  0.40502198E+00  0.63294070E+00</point>
      <point r="1.140">  0.65456935E+00 -0.48373931E+00  0.39178397E+00  0.62882619E+00</point>
      <point r="1.160">  0.65118981E+00 -0.48805394E+00  0.37841851E+00  0.62462421E+00</point>
      <point r="1.180">  0.64782942E+00 -0.49228745E+00  0.36493851E+00  0.62033715E+00</point>
      <point r="1.200">  0.64448559E+00 -0.49643736E+00  0.35135679E+00  0.61596746E+00</point>
      <point r="1.220">  0.64115571E+00 -0.50050112E+00  0.33768606E+00  0.61151764E+00</point>
      <point r="1.240">  0.63783714E+00 -0.50447612E+00  0.32393897E+00  0.60699023E+00</point>
      <point r="1.260">  0.63452724E+00 -0.50835970E+00  0.31012799E+00  0.60238782E+00</point>
      <point r="1.280">  0.63122340E+00 -0.51214922E+00  0.29626549E+00  0.59771302E+00</point>
      <point r="1.300">  0.62792301E+00 -0.51584199E+00  0.28236365E+00  0.59296847E+00</point>
      <point r="1.320">  0.62462353E+00 -0.51943536E+00  0.26843448E+00  0.58815684E+00</point>
      <point r="1.340">  0.62132244E+00 -0.52292671E+00  0.25448979E+00  0.58328081E+00</point>
      <point r="1.360">  0.61801730E+00 -0.52631345E+00  0.24054119E+00  0.57834308E+00</point>
      <point r="1.380">  0.61470572E+00 -0.52959305E+00  0.22660007E+00  0.57334636E+00</point>
      <point r="1.400">  0.61138539E+00 -0.53276302E+00  0.21267761E+00  0.56829337E+00</point>
      <point r="1.420">  0.60805408E+00 -0.53582098E+00  0.19878471E+00  0.56318682E+00</point>
      <point r="1.440">  0.60470965E+00 -0.53876459E+00  0.18493205E+00  0.55802943E+00</point>
      <point r="1.460">  0.60135006E+00 -0.54159162E+00  0.17113005E+00  0.55282391E+00</point>
      <point r="1.480">  0.59797335E+00 -0.54429994E+00  0.15738885E+00  0.54757296E+00</point>
      <point r="1.500">  0.59457767E+00 -0.54688752E+00  0.14371833E+00  0.54227928E+00</point>
      <point r="1.520">  0.59116127E+00 -0.54935241E+00  0.13012809E+00  0.53694556E+00</point>
      <point r="1.540">  0.58772251E+00 -0.55169281E+00  0.11662745E+00  0.53157445E+00</point>
      <point r="1.560">  0.58425984E+00 -0.55390701E+00  0.10322542E+00  0.52616859E+00</point>
      <point r="1.580">  0.58077183E+00 -0.55599344E+00  0.89930741E-01  0.52073062E+00</point>
      <point r="1.600">  0.57725718E+00 -0.55795063E+00  0.76751849E-01  0.51526313E+00</point>
      <point r="1.620">  0.57371466E+00 -0.55977727E+00  0.63696880E-01  0.50976869E+00</point>
      <point r="1.640">  0.57014317E+00 -0.56147213E+00  0.50773672E-01  0.50424985E+00</point>
      <point r="1.660">  0.56654173E+00 -0.56303415E+00  0.37989759E-01  0.49870913E+00</point>
      <point r="1.680">  0.56290945E+00 -0.56446239E+00  0.25352371E-01  0.49314900E+00</point>
      <point r="1.700">  0.55924555E+00 -0.56575602E+00  0.12868436E-01  0.48757193E+00</point>
      <point r="1.720">  0.55554936E+00 -0.56691437E+00  0.54457569E-03  0.48198031E+00</point>
      <point r="1.740">  0.55182031E+00 -0.56793687E+00 -0.11612890E-01  0.47637653E+00</point>
      <point r="1.760">  0.54805795E+00 -0.56882311E+00 -0.23597946E-01  0.47076294E+00</point>
      <point r="1.780">  0.54426191E+00 -0.56957279E+00 -0.35404875E-01  0.46514182E+00</point>
      <point r="1.800">  0.54043192E+00 -0.57018575E+00 -0.47028263E-01  0.45951544E+00</point>
      <point r="1.820">  0.53656782E+00 -0.57066194E+00 -0.58462993E-01  0.45388602E+00</point>
      <point r="1.840">  0.53266954E+00 -0.57100146E+00 -0.69704243E-01  0.44825574E+00</point>
      <point r="1.860">  0.52873708E+00 -0.57120451E+00 -0.80747482E-01  0.44262673E+00</point>
      <point r="1.880">  0.52477056E+00 -0.57127143E+00 -0.91588470E-01  0.43700107E+00</point>
      <point r="1.900">  0.52077017E+00 -0.57120266E+00 -0.10222326E+00  0.43138082E+00</point>
      <point r="1.920">  0.51673617E+00 -0.57099877E+00 -0.11264817E+00  0.42576796E+00</point>
      <point r="1.940">  0.51266893E+00 -0.57066044E+00 -0.12285981E+00  0.42016446E+00</point>
      <point r="1.960">  0.50856886E+00 -0.57018847E+00 -0.13285507E+00  0.41457221E+00</point>
      <point r="1.980">  0.50443648E+00 -0.56958376E+00 -0.14263110E+00  0.40899308E+00</point>
      <point r="2.000">  0.50027236E+00 -0.56884732E+00 -0.15218531E+00  0.40342888E+00</point>
      <point r="2.020">  0.49607713E+00 -0.56798025E+00 -0.16151540E+00  0.39788138E+00</point>
      <point r="2.040">  0.49185150E+00 -0.56698376E+00 -0.17061929E+00  0.39235228E+00</point>
      <point r="2.060">  0.48759625E+00 -0.56585916E+00 -0.17949517E+00  0.38684327E+00</point>
      <point r="2.080">  0.48331218E+00 -0.56460786E+00 -0.18814148E+00  0.38135595E+00</point>
      <point r="2.100">  0.47900018E+00 -0.56323133E+00 -0.19655689E+00  0.37589191E+00</point>
      <point r="2.120">  0.47466118E+00 -0.56173115E+00 -0.20474031E+00  0.37045267E+00</point>
      <point r="2.140">  0.47029615E+00 -0.56010898E+00 -0.21269088E+00  0.36503970E+00</point>
      <point r="2.160">  0.46590612E+00 -0.55836656E+00 -0.22040797E+00  0.35965443E+00</point>
      <point r="2.180">  0.46149215E+00 -0.55650569E+00 -0.22789115E+00  0.35429825E+00</point>
      <point r="2.200">  0.45705536E+00 -0.55452827E+00 -0.23514023E+00  0.34897250E+00</point>
      <point r="2.220">  0.45259687E+00 -0.55243623E+00 -0.24215520E+00  0.34367844E+00</point>
      <point r="2.240">  0.44811787E+00 -0.55023161E+00 -0.24893626E+00  0.33841734E+00</point>
      <point r="2.260">  0.44361955E+00 -0.54791647E+00 -0.25548382E+00  0.33319038E+00</point>
      <point r="2.280">  0.43910316E+00 -0.54549297E+00 -0.26179845E+00  0.32799870E+00</point>
      <point r="2.300">  0.43456994E+00 -0.54296328E+00 -0.26788093E+00  0.32284342E+00</point>
      <point r="2.320">  0.43002118E+00 -0.54032965E+00 -0.27373220E+00  0.31772558E+00</point>
      <point r="2.340">  0.42545818E+00 -0.53759437E+00 -0.27935338E+00  0.31264620E+00</point>
      <point r="2.360">  0.42088226E+00 -0.53475977E+00 -0.28474574E+00  0.30760624E+00</point>
      <point r="2.380">  0.41629474E+00 -0.53182823E+00 -0.28991073E+00  0.30260664E+00</point>
      <point r="2.400">  0.41169697E+00 -0.52880216E+00 -0.29484993E+00  0.29764826E+00</point>
      <point r="2.420">  0.40709030E+00 -0.52568401E+00 -0.29956508E+00  0.29273194E+00</point>
      <point r="2.440">  0.40247609E+00 -0.52247625E+00 -0.30405806E+00  0.28785848E+00</point>
      <point r="2.460">  0.39785570E+00 -0.51918138E+00 -0.30833089E+00  0.28302863E+00</point>
      <point r="2.480">  0.39323052E+00 -0.51580194E+00 -0.31238569E+00  0.27824311E+00</point>
      <point r="2.500">  0.38860190E+00 -0.51234047E+00 -0.31622473E+00  0.27350258E+00</point>
      <point r="2.520">  0.38397122E+00 -0.50879955E+00 -0.31985040E+00  0.26880768E+00</point>
      <point r="2.540">  0.37933983E+00 -0.50518175E+00 -0.32326518E+00  0.26415899E+00</point>
      <point r="2.560">  0.37470911E+00 -0.50148966E+00 -0.32647167E+00  0.25955708E+00</point>
      <point r="2.580">  0.37008041E+00 -0.49772591E+00 -0.32947256E+00  0.25500246E+00</point>
      <point r="2.600">  0.36545507E+00 -0.49389308E+00 -0.33227065E+00  0.25049561E+00</point>
      <point r="2.620">  0.36083444E+00 -0.48999381E+00 -0.33486881E+00  0.24603697E+00</point>
      <point r="2.640">  0.35621983E+00 -0.48603071E+00 -0.33727000E+00  0.24162695E+00</point>
      <point r="2.660">  0.35161256E+00 -0.48200639E+00 -0.33947727E+00  0.23726592E+00</point>
      <point r="2.680">  0.34701392E+00 -0.47792347E+00 -0.34149373E+00  0.23295422E+00</point>
      <point r="2.700">  0.34242521E+00 -0.47378455E+00 -0.34332256E+00  0.22869216E+00</point>
      <point r="2.720">  0.33784768E+00 -0.46959223E+00 -0.34496699E+00  0.22448002E+00</point>
      <point r="2.740">  0.33328259E+00 -0.46534910E+00 -0.34643034E+00  0.22031802E+00</point>
      <point r="2.760">  0.32873117E+00 -0.46105773E+00 -0.34771596E+00  0.21620638E+00</point>
      <point r="2.780">  0.32419461E+00 -0.45672068E+00 -0.34882724E+00  0.21214529E+00</point>
      <point r="2.800">  0.31967411E+00 -0.45234050E+00 -0.34976763E+00  0.20813488E+00</point>
      <point r="2.820">  0.31517084E+00 -0.44791970E+00 -0.35054063E+00  0.20417529E+00</point>
      <point r="2.840">  0.31068594E+00 -0.44346079E+00 -0.35114974E+00  0.20026660E+00</point>
      <point r="2.860">  0.30622051E+00 -0.43896625E+00 -0.35159851E+00  0.19640888E+00</point>
      <point r="2.880">  0.30177567E+00 -0.43443854E+00 -0.35189053E+00  0.19260217E+00</point>
      <point r="2.900">  0.29735248E+00 -0.42988008E+00 -0.35202939E+00  0.18884648E+00</point>
      <point r="2.920">  0.29295197E+00 -0.42529329E+00 -0.35201870E+00  0.18514180E+00</point>
      <point r="2.940">  0.28857517E+00 -0.42068053E+00 -0.35186210E+00  0.18148810E+00</point>
      <point r="2.960">  0.28422307E+00 -0.41604415E+00 -0.35156322E+00  0.17788531E+00</point>
      <point r="2.980">  0.27989662E+00 -0.41138646E+00 -0.35112572E+00  0.17433336E+00</point>
      <point r="3.000">  0.27559677E+00 -0.40670974E+00 -0.35055324E+00  0.17083213E+00</point>
      <point r="3.020">  0.27132441E+00 -0.40201624E+00 -0.34984944E+00  0.16738151E+00</point>
      <point r="3.040">  0.26708043E+00 -0.39730816E+00 -0.34901796E+00  0.16398134E+00</point>
      <point r="3.060">  0.26286567E+00 -0.39258769E+00 -0.34806244E+00  0.16063146E+00</point>
      <point r="3.080">  0.25868095E+00 -0.38785696E+00 -0.34698652E+00  0.15733168E+00</point>
      <point r="3.100">  0.25452707E+00 -0.38311807E+00 -0.34579381E+00  0.15408180E+00</point>
      <point r="3.120">  0.25040479E+00 -0.37837307E+00 -0.34448793E+00  0.15088159E+00</point>
      <point r="3.140">  0.24631484E+00 -0.37362400E+00 -0.34307244E+00  0.14773083E+00</point>
      <point r="3.160">  0.24225792E+00 -0.36887283E+00 -0.34155092E+00  0.14462924E+00</point>
      <point r="3.180">  0.23823471E+00 -0.36412150E+00 -0.33992692E+00  0.14157655E+00</point>
      <point r="3.200">  0.23424585E+00 -0.35937191E+00 -0.33820393E+00  0.13857249E+00</point>
      <point r="3.220">  0.23029195E+00 -0.35462591E+00 -0.33638547E+00  0.13561674E+00</point>
      <point r="3.240">  0.22637361E+00 -0.34988533E+00 -0.33447498E+00  0.13270900E+00</point>
      <point r="3.260">  0.22249138E+00 -0.34515193E+00 -0.33247589E+00  0.12984893E+00</point>
      <point r="3.280">  0.21864580E+00 -0.34042743E+00 -0.33039159E+00  0.12703619E+00</point>
      <point r="3.300">  0.21483735E+00 -0.33571354E+00 -0.32822544E+00  0.12427042E+00</point>
      <point r="3.320">  0.21106652E+00 -0.33101187E+00 -0.32598076E+00  0.12155126E+00</point>
      <point r="3.340">  0.20733375E+00 -0.32632404E+00 -0.32366083E+00  0.11887834E+00</point>
      <point r="3.360">  0.20363944E+00 -0.32165159E+00 -0.32126888E+00  0.11625126E+00</point>
      <point r="3.380">  0.19998401E+00 -0.31699604E+00 -0.31880811E+00  0.11366963E+00</point>
      <point r="3.400">  0.19636779E+00 -0.31235884E+00 -0.31628168E+00  0.11113305E+00</point>
      <point r="3.420">  0.19279113E+00 -0.30774143E+00 -0.31369269E+00  0.10864109E+00</point>
      <point r="3.440">  0.18925435E+00 -0.30314517E+00 -0.31104421E+00  0.10619333E+00</point>
      <point r="3.460">  0.18575771E+00 -0.29857139E+00 -0.30833924E+00  0.10378934E+00</point>
      <point r="3.480">  0.18230148E+00 -0.29402140E+00 -0.30558075E+00  0.10142869E+00</point>
      <point r="3.500">  0.17888588E+00 -0.28949642E+00 -0.30277166E+00  0.99110923E-01</point>
      <point r="3.520">  0.17551113E+00 -0.28499767E+00 -0.29991483E+00  0.96835593E-01</point>
      <point r="3.540">  0.17217741E+00 -0.28052630E+00 -0.29701308E+00  0.94602242E-01</point>
      <point r="3.560">  0.16888487E+00 -0.27608343E+00 -0.29406917E+00  0.92410406E-01</point>
      <point r="3.580">  0.16563365E+00 -0.27167013E+00 -0.29108581E+00  0.90259618E-01</point>
      <point r="3.600">  0.16242387E+00 -0.26728743E+00 -0.28806565E+00  0.88149405E-01</point>
      <point r="3.620">  0.15925561E+00 -0.26293632E+00 -0.28501131E+00  0.86079291E-01</point>
      <point r="3.640">  0.15612893E+00 -0.25861775E+00 -0.28192533E+00  0.84048796E-01</point>
      <point r="3.660">  0.15304389E+00 -0.25433262E+00 -0.27881020E+00  0.82057436E-01</point>
      <point r="3.680">  0.15000052E+00 -0.25008180E+00 -0.27566837E+00  0.80104725E-01</point>
      <point r="3.700">  0.14699880E+00 -0.24586611E+00 -0.27250223E+00  0.78190174E-01</point>
      <point r="3.720">  0.14403874E+00 -0.24168635E+00 -0.26931409E+00  0.76313292E-01</point>
      <point r="3.740">  0.14112030E+00 -0.23754325E+00 -0.26610626E+00  0.74473586E-01</point>
      <point r="3.760">  0.13824341E+00 -0.23343752E+00 -0.26288093E+00  0.72670562E-01</point>
      <point r="3.780">  0.13540802E+00 -0.22936984E+00 -0.25964028E+00  0.70903724E-01</point>
      <point r="3.800">  0.13261402E+00 -0.22534084E+00 -0.25638642E+00  0.69172576E-01</point>
      <point r="3.820">  0.12986132E+00 -0.22135111E+00 -0.25312141E+00  0.67476621E-01</point>
      <point r="3.840">  0.12714979E+00 -0.21740121E+00 -0.24984724E+00  0.65815363E-01</point>
      <point r="3.860">  0.12447929E+00 -0.21349167E+00 -0.24656586E+00  0.64188304E-01</point>
      <point r="3.880">  0.12184967E+00 -0.20962298E+00 -0.24327917E+00  0.62594947E-01</point>
      <point r="3.900">  0.11926075E+00 -0.20579558E+00 -0.23998901E+00  0.61034797E-01</point>
      <point r="3.920">  0.11671234E+00 -0.20200991E+00 -0.23669715E+00  0.59507358E-01</point>
      <point r="3.940">  0.11420426E+00 -0.19826635E+00 -0.23340533E+00  0.58012137E-01</point>
      <point r="3.960">  0.11173628E+00 -0.19456525E+00 -0.23011522E+00  0.56548642E-01</point>
      <point r="3.980">  0.10930817E+00 -0.19090694E+00 -0.22682845E+00  0.55116380E-01</point>
      <point r="4.000">  0.10691971E+00 -0.18729171E+00 -0.22354660E+00  0.53714862E-01</point>
      <point r="4.020">  0.10457064E+00 -0.18371983E+00 -0.22027119E+00  0.52343603E-01</point>
      <point r="4.040">  0.10226069E+00 -0.18019152E+00 -0.21700368E+00  0.51002116E-01</point>
      <point r="4.060">  0.99989602E-01 -0.17670700E+00 -0.21374549E+00  0.49689920E-01</point>
      <point r="4.080">  0.97757079E-01 -0.17326643E+00 -0.21049800E+00  0.48406534E-01</point>
      <point r="4.100">  0.95562829E-01 -0.16986997E+00 -0.20726252E+00  0.47151482E-01</point>
      <point r="4.120">  0.93406550E-01 -0.16651774E+00 -0.20404033E+00  0.45924289E-01</point>
      <point r="4.140">  0.91287928E-01 -0.16320983E+00 -0.20083265E+00  0.44724484E-01</point>
      <point r="4.160">  0.89206642E-01 -0.15994632E+00 -0.19764066E+00  0.43551601E-01</point>
      <point r="4.180">  0.87162360E-01 -0.15672724E+00 -0.19446548E+00  0.42405174E-01</point>
      <point r="4.200">  0.85154746E-01 -0.15355262E+00 -0.19130819E+00  0.41284743E-01</point>
      <point r="4.220">  0.83183453E-01 -0.15042247E+00 -0.18816985E+00  0.40189851E-01</point>
      <point r="4.240">  0.81248130E-01 -0.14733674E+00 -0.18505143E+00  0.39120046E-01</point>
      <point r="4.260">  0.79348418E-01 -0.14429541E+00 -0.18195389E+00  0.38074877E-01</point>
      <point r="4.280">  0.77483951E-01 -0.14129839E+00 -0.17887814E+00  0.37053900E-01</point>
      <point r="4.300">  0.75654359E-01 -0.13834561E+00 -0.17582503E+00  0.36056674E-01</point>
      <point r="4.320">  0.73859266E-01 -0.13543695E+00 -0.17279540E+00  0.35082763E-01</point>
      <point r="4.340">  0.72098292E-01 -0.13257229E+00 -0.16979002E+00  0.34131733E-01</point>
      <point r="4.360">  0.70371050E-01 -0.12975148E+00 -0.16680963E+00  0.33203158E-01</point>
      <point r="4.380">  0.68677151E-01 -0.12697436E+00 -0.16385495E+00  0.32296613E-01</point>
      <point r="4.400">  0.67016202E-01 -0.12424075E+00 -0.16092662E+00  0.31411680E-01</point>
      <point r="4.420">  0.65387806E-01 -0.12155044E+00 -0.15802528E+00  0.30547946E-01</point>
      <point r="4.440">  0.63791562E-01 -0.11890324E+00 -0.15515152E+00  0.29705000E-01</point>
      <point r="4.460">  0.62227069E-01 -0.11629890E+00 -0.15230589E+00  0.28882438E-01</point>
      <point r="4.480">  0.60693920E-01 -0.11373719E+00 -0.14948892E+00  0.28079860E-01</point>
      <point r="4.500">  0.59191709E-01 -0.11121784E+00 -0.14670108E+00  0.27296871E-01</point>
      <point r="4.520">  0.57720026E-01 -0.10874059E+00 -0.14394283E+00  0.26533082E-01</point>
      <point r="4.540">  0.56278460E-01 -0.10630515E+00 -0.14121458E+00  0.25788107E-01</point>
      <point r="4.560">  0.54866600E-01 -0.10391123E+00 -0.13851673E+00  0.25061567E-01</point>
      <point r="4.580">  0.53484033E-01 -0.10155852E+00 -0.13584963E+00  0.24353086E-01</point>
      <point r="4.600">  0.52130345E-01 -0.99246698E-01 -0.13321360E+00  0.23662294E-01</point>
      <point r="4.620">  0.50805122E-01 -0.96975436E-01 -0.13060895E+00  0.22988827E-01</point>
      <point r="4.640">  0.49507950E-01 -0.94744394E-01 -0.12803593E+00  0.22332324E-01</point>
      <point r="4.660">  0.48238415E-01 -0.92553223E-01 -0.12549480E+00  0.21692431E-01</point>
      <point r="4.680">  0.46996104E-01 -0.90401564E-01 -0.12298575E+00  0.21068799E-01</point>
      <point r="4.700">  0.45780602E-01 -0.88289050E-01 -0.12050899E+00  0.20461082E-01</point>
      <point r="4.720">  0.44591498E-01 -0.86215306E-01 -0.11806467E+00  0.19868941E-01</point>
      <point r="4.740">  0.43428380E-01 -0.84179951E-01 -0.11565292E+00  0.19292042E-01</point>
      <point r="4.760">  0.42290838E-01 -0.82182593E-01 -0.11327387E+00  0.18730056E-01</point>
      <point r="4.780">  0.41178463E-01 -0.80222837E-01 -0.11092759E+00  0.18182658E-01</point>
      <point r="4.800">  0.40090848E-01 -0.78300280E-01 -0.10861417E+00  0.17649530E-01</point>
      <point r="4.820">  0.39027587E-01 -0.76414515E-01 -0.10633363E+00  0.17130357E-01</point>
      <point r="4.840">  0.37988278E-01 -0.74565128E-01 -0.10408602E+00  0.16624831E-01</point>
      <point r="4.860">  0.36972519E-01 -0.72751700E-01 -0.10187133E+00  0.16132646E-01</point>
      <point r="4.880">  0.35979909E-01 -0.70973809E-01 -0.99689539E-01  0.15653506E-01</point>
      <point r="4.900">  0.35010054E-01 -0.69231028E-01 -0.97540627E-01  0.15187114E-01</point>
      <point r="4.920">  0.34062559E-01 -0.67522928E-01 -0.95424536E-01  0.14733184E-01</point>
      <point r="4.940">  0.33137032E-01 -0.65849073E-01 -0.93341197E-01  0.14291429E-01</point>
      <point r="4.960">  0.32233085E-01 -0.64209028E-01 -0.91290524E-01  0.13861572E-01</point>
      <point r="4.980">  0.31350332E-01 -0.62602354E-01 -0.89272414E-01  0.13443338E-01</point>
      <point r="5.000">  0.30488391E-01 -0.61028609E-01 -0.87286751E-01  0.13036457E-01</point>
      <point r="5.020">  0.29646883E-01 -0.59487350E-01 -0.85333400E-01  0.12640665E-01</point>
      <point r="5.040">  0.28825431E-01 -0.57978134E-01 -0.83412215E-01  0.12255701E-01</point>
      <point r="5.060">  0.28023663E-01 -0.56500512E-01 -0.81523036E-01  0.11881311E-01</point>
      <point r="5.080">  0.27241210E-01 -0.55054040E-01 -0.79665689E-01  0.11517244E-01</point>
      <point r="5.100">  0.26477707E-01 -0.53638270E-01 -0.77839988E-01  0.11163254E-01</point>
      <point r="5.120">  0.25732792E-01 -0.52252754E-01 -0.76045734E-01  0.10819100E-01</point>
      <point r="5.140">  0.25006107E-01 -0.50897044E-01 -0.74282719E-01  0.10484545E-01</point>
      <point r="5.160">  0.24297297E-01 -0.49570693E-01 -0.72550723E-01  0.10159356E-01</point>
      <point r="5.180">  0.23606012E-01 -0.48273254E-01 -0.70849515E-01  0.98433060E-02</point>
      <point r="5.200">  0.22931907E-01 -0.47004281E-01 -0.69178856E-01  0.95361715E-02</point>
      <point r="5.220">  0.22274638E-01 -0.45763328E-01 -0.67538496E-01  0.92377335E-02</point>
      <point r="5.240">  0.21633867E-01 -0.44549952E-01 -0.65928177E-01  0.89477771E-02</point>
      <point r="5.260">  0.21009261E-01 -0.43363709E-01 -0.64347634E-01  0.86660919E-02</point>
      <point r="5.280">  0.20400490E-01 -0.42204159E-01 -0.62796592E-01  0.83924716E-02</point>
      <point r="5.300">  0.19807226E-01 -0.41070862E-01 -0.61274771E-01  0.81267141E-02</point>
      <point r="5.320">  0.19229150E-01 -0.39963382E-01 -0.59781882E-01  0.78686215E-02</point>
      <point r="5.340">  0.18665944E-01 -0.38881284E-01 -0.58317632E-01  0.76179997E-02</point>
      <point r="5.360">  0.18117294E-01 -0.37824134E-01 -0.56881721E-01  0.73746589E-02</point>
      <point r="5.380">  0.17582892E-01 -0.36791504E-01 -0.55473841E-01  0.71384132E-02</point>
      <point r="5.400">  0.17062434E-01 -0.35782966E-01 -0.54093684E-01  0.69090803E-02</point>
      <point r="5.420">  0.16555619E-01 -0.34798095E-01 -0.52740933E-01  0.66864821E-02</point>
      <point r="5.440">  0.16062153E-01 -0.33836471E-01 -0.51415268E-01  0.64704440E-02</point>
      <point r="5.460">  0.15581742E-01 -0.32897674E-01 -0.50116365E-01  0.62607955E-02</point>
      <point r="5.480">  0.15114102E-01 -0.31981292E-01 -0.48843898E-01  0.60573695E-02</point>
      <point r="5.500">  0.14658949E-01 -0.31086912E-01 -0.47597534E-01  0.58600025E-02</point>
      <point r="5.520">  0.14216004E-01 -0.30214126E-01 -0.46376942E-01  0.56685347E-02</point>
      <point r="5.540">  0.13784995E-01 -0.29362531E-01 -0.45181784E-01  0.54828098E-02</point>
      <point r="5.560">  0.13365651E-01 -0.28531725E-01 -0.44011722E-01  0.53026750E-02</point>
      <point r="5.580">  0.12957707E-01 -0.27721313E-01 -0.42866415E-01  0.51279807E-02</point>
      <point r="5.600">  0.12560904E-01 -0.26930903E-01 -0.41745522E-01  0.49585810E-02</point>
      <point r="5.620">  0.12174983E-01 -0.26160105E-01 -0.40648699E-01  0.47943329E-02</point>
      <point r="5.640">  0.11799693E-01 -0.25408535E-01 -0.39575601E-01  0.46350971E-02</point>
      <point r="5.660">  0.11434787E-01 -0.24675813E-01 -0.38525883E-01  0.44807370E-02</point>
      <point r="5.680">  0.11080020E-01 -0.23961564E-01 -0.37499198E-01  0.43311195E-02</point>
      <point r="5.700">  0.10735153E-01 -0.23265416E-01 -0.36495201E-01  0.41861144E-02</point>
      <point r="5.720">  0.10399951E-01 -0.22587002E-01 -0.35513545E-01  0.40455947E-02</point>
      <point r="5.740">  0.10074183E-01 -0.21925959E-01 -0.34553884E-01  0.39094362E-02</point>
      <point r="5.760">  0.97576221E-02 -0.21281929E-01 -0.33615871E-01  0.37775177E-02</point>
      <point r="5.780">  0.94500460E-02 -0.20654559E-01 -0.32699162E-01  0.36497209E-02</point>
      <point r="5.800">  0.91512360E-02 -0.20043499E-01 -0.31803411E-01  0.35259303E-02</point>
      <point r="5.820">  0.88609776E-02 -0.19448406E-01 -0.30928276E-01  0.34060331E-02</point>
      <point r="5.840">  0.85790602E-02 -0.18868939E-01 -0.30073413E-01  0.32899194E-02</point>
      <point r="5.860">  0.83052776E-02 -0.18304763E-01 -0.29238482E-01  0.31774818E-02</point>
      <point r="5.880">  0.80394272E-02 -0.17755548E-01 -0.28423142E-01  0.30686156E-02</point>
      <point r="5.900">  0.77813104E-02 -0.17220969E-01 -0.27627056E-01  0.29632186E-02</point>
      <point r="5.920">  0.75307326E-02 -0.16700705E-01 -0.26849889E-01  0.28611911E-02</point>
      <point r="5.940">  0.72875030E-02 -0.16194439E-01 -0.26091304E-01  0.27624360E-02</point>
      <point r="5.960">  0.70514346E-02 -0.15701860E-01 -0.25350971E-01  0.26668586E-02</point>
      <point r="5.980">  0.68223441E-02 -0.15222661E-01 -0.24628560E-01  0.25743665E-02</point>
      <point r="6.000">  0.66000519E-02 -0.14756541E-01 -0.23923743E-01  0.24848696E-02</point>
      <point r="6.020">  0.63843823E-02 -0.14303202E-01 -0.23236195E-01  0.23982802E-02</point>
      <point r="6.040">  0.61751630E-02 -0.13862352E-01 -0.22565594E-01  0.23145127E-02</point>
      <point r="6.060">  0.59722255E-02 -0.13433703E-01 -0.21911620E-01  0.22334839E-02</point>
      <point r="6.080">  0.57754045E-02 -0.13016972E-01 -0.21273956E-01  0.21551124E-02</point>
      <point r="6.100">  0.55845386E-02 -0.12611881E-01 -0.20652289E-01  0.20793192E-02</point>
      <point r="6.120">  0.53994696E-02 -0.12218157E-01 -0.20046308E-01  0.20060273E-02</point>
      <point r="6.140">  0.52200428E-02 -0.11835530E-01 -0.19455704E-01  0.19351617E-02</point>
      <point r="6.160">  0.50461068E-02 -0.11463737E-01 -0.18880173E-01  0.18666491E-02</point>
      <point r="6.180">  0.48775135E-02 -0.11102517E-01 -0.18319414E-01  0.18004186E-02</point>
      <point r="6.200">  0.47141182E-02 -0.10751616E-01 -0.17773128E-01  0.17364009E-02</point>
      <point r="6.220">  0.45557794E-02 -0.10410783E-01 -0.17241021E-01  0.16745285E-02</point>
      <point r="6.240">  0.44023585E-02 -0.10079773E-01 -0.16722801E-01  0.16147359E-02</point>
      <point r="6.260">  0.42537203E-02 -0.97583435E-02 -0.16218180E-01  0.15569592E-02</point>
      <point r="6.280">  0.41097326E-02 -0.94462580E-02 -0.15726874E-01  0.15011364E-02</point>
      <point r="6.300">  0.39702662E-02 -0.91432841E-02 -0.15248603E-01  0.14472068E-02</point>
      <point r="6.320">  0.38351949E-02 -0.88491937E-02 -0.14783088E-01  0.13951119E-02</point>
      <point r="6.340">  0.37043954E-02 -0.85637632E-02 -0.14330057E-01  0.13447944E-02</point>
      <point r="6.360">  0.35777474E-02 -0.82867730E-02 -0.13889239E-01  0.12961987E-02</point>
      <point r="6.380">  0.34551332E-02 -0.80180083E-02 -0.13460368E-01  0.12492707E-02</point>
      <point r="6.400">  0.33364380E-02 -0.77572580E-02 -0.13043182E-01  0.12039579E-02</point>
      <point r="6.420">  0.32215499E-02 -0.75043156E-02 -0.12637422E-01  0.11602093E-02</point>
      <point r="6.440">  0.31103594E-02 -0.72589786E-02 -0.12242833E-01  0.11179751E-02</point>
      <point r="6.460">  0.30027598E-02 -0.70210486E-02 -0.11859164E-01  0.10772070E-02</point>
      <point r="6.480">  0.28986471E-02 -0.67903313E-02 -0.11486167E-01  0.10378583E-02</point>
      <point r="6.500">  0.27979196E-02 -0.65666365E-02 -0.11123599E-01  0.99988328E-03</point>
      <point r="6.520">  0.27004783E-02 -0.63497777E-02 -0.10771220E-01  0.96323772E-03</point>
      <point r="6.540">  0.26062266E-02 -0.61395727E-02 -0.10428793E-01  0.92787862E-03</point>
      <point r="6.560">  0.25150703E-02 -0.59358429E-02 -0.10096087E-01  0.89376423E-03</point>
      <point r="6.580">  0.24269175E-02 -0.57384136E-02 -0.97728729E-02  0.86085401E-03</point>
      <point r="6.600">  0.23416787E-02 -0.55471139E-02 -0.94589260E-02  0.82910859E-03</point>
      <point r="6.620">  0.22592668E-02 -0.53617766E-02 -0.91540255E-02  0.79848974E-03</point>
      <point r="6.640">  0.21795966E-02 -0.51822382E-02 -0.88579541E-02  0.76896036E-03</point>
      <point r="6.660">  0.21025853E-02 -0.50083387E-02 -0.85704984E-02  0.74048445E-03</point>
      <point r="6.680">  0.20281524E-02 -0.48399219E-02 -0.82914487E-02  0.71302707E-03</point>
      <point r="6.700">  0.19562192E-02 -0.46768349E-02 -0.80205990E-02  0.68655430E-03</point>
      <point r="6.720">  0.18867091E-02 -0.45189284E-02 -0.77577471E-02  0.66103327E-03</point>
      <point r="6.740">  0.18195478E-02 -0.43660564E-02 -0.75026942E-02  0.63643204E-03</point>
      <point r="6.760">  0.17546626E-02 -0.42180764E-02 -0.72552454E-02  0.61271969E-03</point>
      <point r="6.780">  0.16919830E-02 -0.40748491E-02 -0.70152094E-02  0.58986618E-03</point>
      <point r="6.800">  0.16314403E-02 -0.39362386E-02 -0.67823984E-02  0.56784243E-03</point>
      <point r="6.820">  0.15729676E-02 -0.38021119E-02 -0.65566283E-02  0.54662021E-03</point>
      <point r="6.840">  0.15164999E-02 -0.36723394E-02 -0.63377183E-02  0.52617217E-03</point>
      <point r="6.860">  0.14619739E-02 -0.35467947E-02 -0.61254914E-02  0.50647180E-03</point>
      <point r="6.880">  0.14093280E-02 -0.34253541E-02 -0.59197739E-02  0.48749340E-03</point>
      <point r="6.900">  0.13585026E-02 -0.33078972E-02 -0.57203956E-02  0.46921208E-03</point>
      <point r="6.920">  0.13094393E-02 -0.31943064E-02 -0.55271897E-02  0.45160371E-03</point>
      <point r="6.940">  0.12620817E-02 -0.30844671E-02 -0.53399929E-02  0.43464491E-03</point>
      <point r="6.960">  0.12163747E-02 -0.29782673E-02 -0.51586450E-02  0.41831307E-03</point>
      <point r="6.980">  0.11722650E-02 -0.28755981E-02 -0.49829893E-02  0.40258624E-03</point>
      <point r="7.000">  0.11297006E-02 -0.27763531E-02 -0.48128723E-02  0.38744321E-03</point>
      <point r="7.020">  0.10886311E-02 -0.26804287E-02 -0.46481437E-02  0.37286341E-03</point>
      <point r="7.040">  0.10490076E-02 -0.25877239E-02 -0.44886565E-02  0.35882695E-03</point>
      <point r="7.060">  0.10107824E-02 -0.24981403E-02 -0.43342669E-02  0.34531458E-03</point>
      <point r="7.080">  0.97390943E-03 -0.24115821E-02 -0.41848339E-02  0.33230763E-03</point>
      <point r="7.100">  0.93834377E-03 -0.23279559E-02 -0.40402199E-02  0.31978809E-03</point>
      <point r="7.120">  0.90404191E-03 -0.22471707E-02 -0.39002902E-02  0.30773849E-03</point>
      <point r="7.140">  0.87096157E-03 -0.21691381E-02 -0.37649131E-02  0.29614194E-03</point>
      <point r="7.160">  0.83906177E-03 -0.20937717E-02 -0.36339598E-02  0.28498212E-03</point>
      <point r="7.180">  0.80830270E-03 -0.20209877E-02 -0.35073046E-02  0.27424323E-03</point>
      <point r="7.200">  0.77864578E-03 -0.19507045E-02 -0.33848244E-02  0.26391000E-03</point>
      <point r="7.220">  0.75005357E-03 -0.18828426E-02 -0.32663992E-02  0.25396765E-03</point>
      <point r="7.240">  0.72248974E-03 -0.18173247E-02 -0.31519116E-02  0.24440190E-03</point>
      <point r="7.260">  0.69591911E-03 -0.17540756E-02 -0.30412469E-02  0.23519896E-03</point>
      <point r="7.280">  0.67030753E-03 -0.16930222E-02 -0.29342932E-02  0.22634549E-03</point>
      <point r="7.300">  0.64562192E-03 -0.16340934E-02 -0.28309413E-02  0.21782859E-03</point>
      <point r="7.320">  0.62183021E-03 -0.15772201E-02 -0.27310845E-02  0.20963582E-03</point>
      <point r="7.340">  0.59890132E-03 -0.15223351E-02 -0.26346188E-02  0.20175514E-03</point>
      <point r="7.360">  0.57680516E-03 -0.14693732E-02 -0.25414425E-02  0.19417493E-03</point>
      <point r="7.380">  0.55551256E-03 -0.14182708E-02 -0.24514566E-02  0.18688398E-03</point>
      <point r="7.400">  0.53499527E-03 -0.13689664E-02 -0.23645645E-02  0.17987144E-03</point>
      <point r="7.420">  0.51522594E-03 -0.13214002E-02 -0.22806719E-02  0.17312684E-03</point>
      <point r="7.440">  0.49617808E-03 -0.12755139E-02 -0.21996869E-02  0.16664010E-03</point>
      <point r="7.460">  0.47782606E-03 -0.12312512E-02 -0.21215200E-02  0.16040146E-03</point>
      <point r="7.480">  0.46014506E-03 -0.11885572E-02 -0.20460839E-02  0.15440152E-03</point>
      <point r="7.500">  0.44311107E-03 -0.11473788E-02 -0.19732935E-02  0.14863119E-03</point>
      <point r="7.520">  0.42670084E-03 -0.11076644E-02 -0.19030659E-02  0.14308172E-03</point>
      <point r="7.540">  0.41089190E-03 -0.10693638E-02 -0.18353204E-02  0.13774466E-03</point>
      <point r="7.560">  0.39566251E-03 -0.10324286E-02 -0.17699784E-02  0.13261186E-03</point>
      <point r="7.580">  0.38099163E-03 -0.99681160E-03 -0.17069634E-02  0.12767547E-03</point>
      <point r="7.600">  0.36685895E-03 -0.96246710E-03 -0.16462007E-02  0.12292791E-03</point>
      <point r="7.620">  0.35324479E-03 -0.92935082E-03 -0.15876178E-02  0.11836188E-03</point>
      <point r="7.640">  0.34013017E-03 -0.89741984E-03 -0.15311442E-02  0.11397034E-03</point>
      <point r="7.660">  0.32749672E-03 -0.86663252E-03 -0.14767112E-02  0.10974650E-03</point>
      <point r="7.680">  0.31532670E-03 -0.83694855E-03 -0.14242519E-02  0.10568382E-03</point>
      <point r="7.700">  0.30360297E-03 -0.80832887E-03 -0.13737012E-02  0.10177602E-03</point>
      <point r="7.720">  0.29230897E-03 -0.78073562E-03 -0.13249961E-02  0.98017023E-04</point>
      <point r="7.740">  0.28142872E-03 -0.75413217E-03 -0.12780751E-02  0.94400985E-04</point>
      <point r="7.760">  0.27094678E-03 -0.72848302E-03 -0.12328784E-02  0.90922283E-04</point>
      <point r="7.780">  0.26084825E-03 -0.70375384E-03 -0.11893479E-02  0.87575499E-04</point>
      <point r="7.800">  0.25111873E-03 -0.67991136E-03 -0.11474273E-02  0.84355419E-04</point>
      <point r="7.820">  0.24174435E-03 -0.65692341E-03 -0.11070618E-02  0.81257023E-04</point>
      <point r="7.840">  0.23271172E-03 -0.63475887E-03 -0.10681980E-02  0.78275483E-04</point>
      <point r="7.860">  0.22400791E-03 -0.61338762E-03 -0.10307843E-02  0.75406150E-04</point>
      <point r="7.880">  0.21562045E-03 -0.59278054E-03 -0.99477058E-03  0.72644553E-04</point>
      <point r="7.900">  0.20753734E-03 -0.57290947E-03 -0.96010799E-03  0.69986393E-04</point>
      <point r="7.920">  0.19974697E-03 -0.55374719E-03 -0.92674927E-03  0.67427534E-04</point>
      <point r="7.940">  0.19223818E-03 -0.53526738E-03 -0.89464853E-03  0.64963998E-04</point>
      <point r="7.960">  0.18500020E-03 -0.51744461E-03 -0.86376124E-03  0.62591964E-04</point>
      <point r="7.980">  0.17802265E-03 -0.50025432E-03 -0.83404420E-03  0.60307757E-04</point>
      <point r="8.000">  0.17129555E-03 -0.48367278E-03 -0.80545553E-03  0.58107843E-04</point>
      <point r="8.020">  0.16480925E-03 -0.46767707E-03 -0.77795462E-03  0.55988831E-04</point>
      <point r="8.040">  0.15855450E-03 -0.45224507E-03 -0.75150209E-03  0.53947459E-04</point>
      <point r="8.060">  0.15252234E-03 -0.43735540E-03 -0.72605978E-03  0.51980596E-04</point>
      <point r="8.080">  0.14670420E-03 -0.42298745E-03 -0.70159071E-03  0.50085234E-04</point>
      <point r="8.100">  0.14109179E-03 -0.40912133E-03 -0.67805907E-03  0.48258485E-04</point>
      <point r="8.120">  0.13567716E-03 -0.39573785E-03 -0.65543015E-03  0.46497575E-04</point>
      <point r="8.140">  0.13045263E-03 -0.38281849E-03 -0.63367035E-03  0.44799844E-04</point>
      <point r="8.160">  0.12541085E-03 -0.37034541E-03 -0.61274714E-03  0.43162736E-04</point>
      <point r="8.180">  0.12054471E-03 -0.35830138E-03 -0.59262902E-03  0.41583799E-04</point>
      <point r="8.200">  0.11584741E-03 -0.34666983E-03 -0.57328552E-03  0.40060683E-04</point>
      <point r="8.220">  0.11131239E-03 -0.33543476E-03 -0.55468716E-03  0.38591130E-04</point>
      <point r="8.240">  0.10693335E-03 -0.32458078E-03 -0.53680539E-03  0.37172977E-04</point>
      <point r="8.260">  0.10270424E-03 -0.31409305E-03 -0.51961265E-03  0.35804150E-04</point>
      <point r="8.280">  0.98619239E-04 -0.30395728E-03 -0.50308224E-03  0.34482659E-04</point>
      <point r="8.300">  0.94672759E-04 -0.29415972E-03 -0.48718839E-03  0.33206596E-04</point>
      <point r="8.320">  0.90859435E-04 -0.28468713E-03 -0.47190617E-03  0.31974135E-04</point>
      <point r="8.340">  0.87174111E-04 -0.27552676E-03 -0.45721149E-03  0.30783523E-04</point>
      <point r="8.360">  0.83611840E-04 -0.26666636E-03 -0.44308110E-03  0.29633084E-04</point>
      <point r="8.380">  0.80167872E-04 -0.25809414E-03 -0.42949253E-03  0.28521208E-04</point>
      <point r="8.400">  0.76837648E-04 -0.24979875E-03 -0.41642407E-03  0.27446358E-04</point>
      <point r="8.420">  0.73616796E-04 -0.24176930E-03 -0.40385480E-03  0.26407059E-04</point>
      <point r="8.440">  0.70501118E-04 -0.23399531E-03 -0.39176449E-03  0.25401899E-04</point>
      <point r="8.460">  0.67486591E-04 -0.22646670E-03 -0.38013365E-03  0.24429526E-04</point>
      <point r="8.480">  0.64569356E-04 -0.21917381E-03 -0.36894347E-03  0.23488649E-04</point>
      <point r="8.500">  0.61745711E-04 -0.21210735E-03 -0.35817580E-03  0.22578028E-04</point>
      <point r="8.520">  0.59012112E-04 -0.20525839E-03 -0.34781315E-03  0.21696480E-04</point>
      <point r="8.540">  0.56365157E-04 -0.19861837E-03 -0.33783866E-03  0.20842872E-04</point>
      <point r="8.560">  0.53801591E-04 -0.19217907E-03 -0.32823610E-03  0.20016120E-04</point>
      <point r="8.580">  0.51318292E-04 -0.18593262E-03 -0.31898980E-03  0.19215187E-04</point>
      <point r="8.600">  0.48912273E-04 -0.17987143E-03 -0.31008471E-03  0.18439083E-04</point>
      <point r="8.620">  0.46580670E-04 -0.17398826E-03 -0.30150629E-03  0.17686861E-04</point>
      <point r="8.640">  0.44320743E-04 -0.16827616E-03 -0.29324060E-03  0.16957613E-04</point>
      <point r="8.660">  0.42129869E-04 -0.16272846E-03 -0.28527419E-03  0.16250476E-04</point>
      <point r="8.680">  0.40005536E-04 -0.15733877E-03 -0.27759412E-03  0.15564620E-04</point>
      <point r="8.700">  0.37945340E-04 -0.15210098E-03 -0.27018796E-03  0.14899256E-04</point>
      <point r="8.720">  0.35946983E-04 -0.14700923E-03 -0.26304375E-03  0.14253628E-04</point>
      <point r="8.740">  0.34008265E-04 -0.14205791E-03 -0.25615001E-03  0.13627015E-04</point>
      <point r="8.760">  0.32127081E-04 -0.13724165E-03 -0.24949569E-03  0.13018726E-04</point>
      <point r="8.780">  0.30301419E-04 -0.13255532E-03 -0.24307019E-03  0.12428104E-04</point>
      <point r="8.800">  0.28529355E-04 -0.12799399E-03 -0.23686332E-03  0.11854520E-04</point>
      <point r="8.820">  0.26809050E-04 -0.12355297E-03 -0.23086531E-03  0.11297373E-04</point>
      <point r="8.840">  0.25138744E-04 -0.11922777E-03 -0.22506677E-03  0.10756090E-04</point>
      <point r="8.860">  0.23516757E-04 -0.11501410E-03 -0.21945872E-03  0.10230123E-04</point>
      <point r="8.880">  0.21941482E-04 -0.11090784E-03 -0.21403252E-03  0.97189499E-05</point>
      <point r="8.900">  0.20411384E-04 -0.10690509E-03 -0.20877989E-03  0.92220712E-05</point>
      <point r="8.920">  0.18924995E-04 -0.10300209E-03 -0.20369292E-03  0.87390106E-05</point>
      <point r="8.940">  0.17480915E-04 -0.99195267E-04 -0.19876402E-03  0.82693135E-05</point>
      <point r="8.960">  0.16077804E-04 -0.95481220E-04 -0.19398590E-03  0.78125457E-05</point>
      <point r="8.980">  0.14714383E-04 -0.91856689E-04 -0.18935162E-03  0.73682926E-05</point>
      <point r="9.000">  0.13389431E-04 -0.88318569E-04 -0.18485452E-03  0.69361589E-05</point>
      <point r="9.020">  0.12101780E-04 -0.84863900E-04 -0.18048823E-03  0.65157666E-05</point>
      <point r="9.040">  0.10850315E-04 -0.81489859E-04 -0.17624667E-03  0.61067554E-05</point>
      <point r="9.060">  0.96339720E-05 -0.78193754E-04 -0.17212401E-03  0.57087810E-05</point>
      <point r="9.080">  0.84517332E-05 -0.74973023E-04 -0.16811471E-03  0.53215146E-05</point>
      <point r="9.100">  0.73026273E-05 -0.71825222E-04 -0.16421345E-03  0.49446425E-05</point>
      <point r="9.120">  0.61857262E-05 -0.68748026E-04 -0.16041518E-03  0.45778650E-05</point>
      <point r="9.140">  0.51001434E-05 -0.65739219E-04 -0.15671507E-03  0.42208956E-05</point>
      <point r="9.160">  0.40450317E-05 -0.62796694E-04 -0.15310852E-03  0.38734609E-05</point>
      <point r="9.180">  0.30195818E-05 -0.59918444E-04 -0.14959113E-03  0.35352996E-05</point>
      <point r="9.200">  0.20230202E-05 -0.57102560E-04 -0.14615873E-03  0.32061620E-05</point>
      <point r="9.220">  0.10546076E-05 -0.54347226E-04 -0.14280735E-03  0.28858092E-05</point>
      <point r="9.240">  0.11363728E-06 -0.51650715E-04 -0.13953319E-03  0.25740132E-05</point>
      <point r="9.260"> -0.80056659E-06 -0.49011384E-04 -0.13633267E-03  0.22705554E-05</point>
      <point r="9.280"> -0.16886504E-05 -0.46427672E-04 -0.13320235E-03  0.19752270E-05</point>
      <point r="9.300"> -0.25512325E-05 -0.43898095E-04 -0.13013900E-03  0.16878282E-05</point>
      <point r="9.320"> -0.33889047E-05 -0.41421242E-04 -0.12713953E-03  0.14081675E-05</point>
      <point r="9.340"> -0.42022336E-05 -0.38995772E-04 -0.12420101E-03  0.11360616E-05</point>
      <point r="9.360"> -0.49917618E-05 -0.36620413E-04 -0.12132067E-03  0.87133471E-06</point>
      <point r="9.380"> -0.57580092E-05 -0.34293952E-04 -0.11849589E-03  0.61381853E-06</point>
      <point r="9.400"> -0.65014743E-05 -0.32015242E-04 -0.11572417E-03  0.36335145E-06</point>
      <point r="9.420"> -0.72226349E-05 -0.29783189E-04 -0.11300317E-03  0.11977843E-06</point>
      <point r="9.440"> -0.79219497E-05 -0.27596755E-04 -0.11033066E-03 -0.11704942E-06</point>
      <point r="9.460"> -0.85998590E-05 -0.25454954E-04 -0.10770454E-03 -0.34727516E-06</point>
      <point r="9.480"> -0.92567860E-05 -0.23356849E-04 -0.10512283E-03 -0.57103632E-06</point>
      <point r="9.500"> -0.98931373E-05 -0.21301551E-04 -0.10258366E-03 -0.78846524E-06</point>
      <point r="9.520"> -0.10509304E-04 -0.19288212E-04 -0.10008526E-03 -0.99968933E-06</point>
      <point r="9.540"> -0.11105663E-04 -0.17316028E-04 -0.97625974E-04 -0.12048313E-05</point>
      <point r="9.560"> -0.11682577E-04 -0.15384234E-04 -0.95204248E-04 -0.14040096E-05</point>
      <point r="9.580"> -0.12240397E-04 -0.13492103E-04 -0.92818614E-04 -0.15973383E-05</point>
      <point r="9.600"> -0.12779459E-04 -0.11638942E-04 -0.90467697E-04 -0.17849279E-05</point>
      <point r="9.620"> -0.13300090E-04 -0.98240912E-05 -0.88150208E-04 -0.19668848E-05</point>
      <point r="9.640"> -0.13802606E-04 -0.80469242E-05 -0.85864938E-04 -0.21433123E-05</point>
      <point r="9.660"> -0.14287312E-04 -0.63068421E-05 -0.83610760E-04 -0.23143102E-05</point>
      <point r="9.680"> -0.14754503E-04 -0.46032748E-05 -0.81386618E-04 -0.24799754E-05</point>
      <point r="9.700"> -0.15204466E-04 -0.29356779E-05 -0.79191528E-04 -0.26404016E-05</point>
      <point r="9.720"> -0.15637479E-04 -0.13035316E-05 -0.77024574E-04 -0.27956800E-05</point>
      <point r="9.740"> -0.16053813E-04  0.29366067E-06 -0.74884903E-04 -0.29458994E-05</point>
      <point r="9.760"> -0.16453732E-04  0.18563742E-05 -0.72771724E-04 -0.30911458E-05</point>
      <point r="9.780"> -0.16837490E-04  0.33850640E-05 -0.70684302E-04 -0.32315031E-05</point>
      <point r="9.800"> -0.17205337E-04  0.48801664E-05 -0.68621960E-04 -0.33670533E-05</point>
      <point r="9.820"> -0.17557519E-04  0.63420998E-05 -0.66584069E-04 -0.34978760E-05</point>
      <point r="9.840"> -0.17894271E-04  0.77712668E-05 -0.64570052E-04 -0.36240492E-05</point>
      <point r="9.860"> -0.18215828E-04  0.91680544E-05 -0.62579377E-04 -0.37456491E-05</point>
      <point r="9.880"> -0.18522417E-04  0.10532836E-04 -0.60611557E-04 -0.38627501E-05</point>
      <point r="9.900"> -0.18814261E-04  0.11865970E-04 -0.58666145E-04 -0.39754254E-05</point>
      <point r="9.920"> -0.19091580E-04  0.13167806E-04 -0.56742734E-04 -0.40837464E-05</point>
      <point r="9.940"> -0.19354590E-04  0.14438681E-04 -0.54840955E-04 -0.41877832E-05</point>
      <point r="9.960"> -0.19603501E-04  0.15678920E-04 -0.52960471E-04 -0.42876048E-05</point>
      <point r="9.980"> -0.19838522E-04  0.16888840E-04 -0.51100981E-04 -0.43832789E-05</point>
      <point r="10.000"> -0.20059858E-04  0.18068751E-04 -0.49262210E-04 -0.44748721E-05</point>
      <point r="10.020"> -0.20267711E-04  0.19218953E-04 -0.47443916E-04 -0.45624500E-05</point>
      <point r="10.040"> -0.20462282E-04  0.20339738E-04 -0.45645882E-04 -0.46460770E-05</point>
      <point r="10.060"> -0.20643768E-04  0.21431395E-04 -0.43867916E-04 -0.47258169E-05</point>
      <point r="10.080"> -0.20812363E-04  0.22494204E-04 -0.42109849E-04 -0.48017324E-05</point>
      <point r="10.100"> -0.20968260E-04  0.23528442E-04 -0.40371534E-04 -0.48738855E-05</point>
      <point r="10.120"> -0.21111650E-04  0.24534379E-04 -0.38652844E-04 -0.49423374E-05</point>
      <point r="10.140"> -0.21242724E-04  0.25512282E-04 -0.36953672E-04 -0.50071486E-05</point>
      <point r="10.160"> -0.21361667E-04  0.26462416E-04 -0.35273926E-04 -0.50683790E-05</point>
      <point r="10.180"> -0.21468667E-04  0.27385040E-04 -0.33613531E-04 -0.51260878E-05</point>
      <point r="10.200"> -0.21563908E-04  0.28280412E-04 -0.31972428E-04 -0.51803334E-05</point>
      <point r="10.220"> -0.21647573E-04  0.29148787E-04 -0.30350568E-04 -0.52311740E-05</point>
      <point r="10.240"> -0.21719846E-04  0.29990419E-04 -0.28747917E-04 -0.52786671E-05</point>
      <point r="10.260"> -0.21780908E-04  0.30805559E-04 -0.27164452E-04 -0.53228697E-05</point>
      <point r="10.280"> -0.21830938E-04  0.31594457E-04 -0.25600158E-04 -0.53638383E-05</point>
      <point r="10.300"> -0.21870116E-04  0.32357363E-04 -0.24055031E-04 -0.54016289E-05</point>
      <point r="10.320"> -0.21898622E-04  0.33094526E-04 -0.22529073E-04 -0.54362971E-05</point>
      <point r="10.340"> -0.21916632E-04  0.33806194E-04 -0.21022295E-04 -0.54678983E-05</point>
      <point r="10.360"> -0.21924325E-04  0.34492615E-04 -0.19534714E-04 -0.54964871E-05</point>
      <point r="10.380"> -0.21921876E-04  0.35154037E-04 -0.18066351E-04 -0.55221180E-05</point>
      <point r="10.400">  0.00000000E+00</point>
    </S_spline>
    <Vrep_spline>
      <point r="0.000">  0.37750294E+03</point>
      <point r="0.100">  0.28193386E+03</point>
      <point r="0.200">  0.21033306E+03</point>
      <point r="0.300">  0.15668940E+03</point>
      <point r="0.400">  0.11649933E+03</point>
      <point r="0.500">  0.86388733E+02</point>
      <point r="0.600">  0.63829735E+02</point>
      <point r="0.700">  0.46928429E+02</point>
      <point r="0.800">  0.34265891E+02</point>
      <point r="0.900">  0.24779059E+02</point>
      <point r="1.000">  0.17671480E+02</point>
      <point r="1.058">  0.14395043E+02</point>
      <point r="1.117">  0.11723444E+02</point>
      <point r="1.175">  0.95534825E+01</point>
      <point r="1.234">  0.77972440E+01</point>
      <point r="1.292">  0.63802484E+01</point>
      <point r="1.351">  0.52397649E+01</point>
      <point r="1.409">  0.43232896E+01</point>
      <point r="1.468">  0.35871732E+01</point>
      <point r="1.526">  0.29953890E+01</point>
      <point r="1.585">  0.25184317E+01</point>
      <point r="1.643">  0.21323373E+01</point>
      <point r="1.701">  0.18178159E+01</point>
      <point r="1.760">  0.15594878E+01</point>
      <point r="1.818">  0.13452147E+01</point>
      <point r="1.877">  0.11655188E+01</point>
      <point r="1.935">  0.10130810E+01</point>
      <point r="1.994">  0.88231199E+00</point>
      <point r="2.052">  0.76898840E+00</point>
      <point r="2.111">  0.66994751E+00</point>
      <point r="2.169">  0.58283470E+00</point>
      <point r="2.228">  0.50589712E+00</point>
      <point r="2.286">  0.43781821E+00</point>
      <point r="2.344">  0.37758773E+00</point>
      <point r="2.403">  0.32440228E+00</point>
      <point r="2.461">  0.27759177E+00</point>
      <point r="2.520">  0.23656748E+00</point>
      <point r="2.578">  0.20078787E+00</point>
      <point r="2.637">  0.16973837E+00</point>
      <point r="2.695">  0.14292205E+00</point>
      <point r="2.754">  0.11985815E+00</point>
      <point r="2.812">  0.10008580E+00</point>
      <point r="2.870">  0.83170890E-01</point>
      <point r="2.929">  0.68713934E-01</point>
      <point r="2.987">  0.56357484E-01</point>
      <point r="3.046">  0.45791817E-01</point>
      <point r="3.104">  0.36757980E-01</point>
      <point r="3.163">  0.29047634E-01</point>
      <point r="3.221">  0.22499458E-01</point>
      <point r="3.280">  0.16992229E-01</point>
      <point r="3.338">  0.12435013E-01</point>
      <point r="3.397">  0.87552610E-02</point>
      <point r="3.455">  0.58859192E-02</point>
      <point r="3.513">  0.37530168E-02</point>
      <point r="3.572">  0.22655281E-02</point>
      <point r="3.630">  0.13096400E-02</point>
      <point r="3.689">  0.74989755E-03</point>
      <point r="3.747">  0.44003621E-03</point>
      <point r="3.806">  0.24664847E-03</point>
      <point r="3.864">  0.89168817E-04</point>
      <point r="3.923">  0.00000000E+00</point>
    </Vrep_spline>
  </per_pair_data>
  <per_pair_data type1="2" type2="1" SK_cutoff="10.4000000000000004"
    Vrep_cutoff="3.9226500000000000" SK_npts="521" Vrep_npts="61">
    <H_spline>
      <point r="0.000">  0.39226500E+01  0.00000000E+00  0.11682400E+01 -0.42589700E+00  0.61469600E-01  0.10886800E+01 -0.15677300E+01  0.00000000E+00  0.39891300E-01 -0.29448900E+00</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400"> -0.87819926E+00 -0.42277802E+00 -0.41515211E+00 -0.60566109E-01</point>
      <point r="0.420"> -0.86277382E+00 -0.44076289E+00 -0.44333612E+00 -0.78141717E-01</point>
      <point r="0.440"> -0.84534046E+00 -0.45688711E+00 -0.46938649E+00 -0.95351112E-01</point>
      <point r="0.460"> -0.82636805E+00 -0.47112007E+00 -0.49327114E+00 -0.11214844E+00</point>
      <point r="0.480"> -0.80627264E+00 -0.48345184E+00 -0.51497917E+00 -0.12849461E+00</point>
      <point r="0.500"> -0.78541999E+00 -0.49389051E+00 -0.53451735E+00 -0.14435665E+00</point>
      <point r="0.520"> -0.76412878E+00 -0.50245944E+00 -0.55190696E+00 -0.15970719E+00</point>
      <point r="0.540"> -0.74267422E+00 -0.50919562E+00 -0.56718169E+00 -0.17452396E+00</point>
      <point r="0.560"> -0.72129134E+00 -0.51414727E+00 -0.58038549E+00 -0.18878925E+00</point>
      <point r="0.580"> -0.70017849E+00 -0.51737108E+00 -0.59157070E+00 -0.20248945E+00</point>
      <point r="0.600"> -0.67950210E+00 -0.51893157E+00 -0.60079893E+00 -0.21561490E+00</point>
      <point r="0.620"> -0.65939994E+00 -0.51889933E+00 -0.60813930E+00 -0.22815919E+00</point>
      <point r="0.640"> -0.63998286E+00 -0.51734777E+00 -0.61366366E+00 -0.24011863E+00</point>
      <point r="0.660"> -0.62133803E+00 -0.51435304E+00 -0.61744606E+00 -0.25149234E+00</point>
      <point r="0.680"> -0.60353313E+00 -0.50999568E+00 -0.61956525E+00 -0.26228216E+00</point>
      <point r="0.700"> -0.58661858E+00 -0.50435932E+00 -0.62010387E+00 -0.27249184E+00</point>
      <point r="0.720"> -0.57062861E+00 -0.49752729E+00 -0.61914437E+00 -0.28212669E+00</point>
      <point r="0.740"> -0.55558297E+00 -0.48958095E+00 -0.61676712E+00 -0.29119372E+00</point>
      <point r="0.760"> -0.54148947E+00 -0.48060045E+00 -0.61305190E+00 -0.29970171E+00</point>
      <point r="0.780"> -0.52834646E+00 -0.47066609E+00 -0.60808030E+00 -0.30766070E+00</point>
      <point r="0.800"> -0.51614466E+00 -0.45985829E+00 -0.60193605E+00 -0.31508170E+00</point>
      <point r="0.820"> -0.50487097E+00 -0.44825825E+00 -0.59471334E+00 -0.32197673E+00</point>
      <point r="0.840"> -0.49449886E+00 -0.43593661E+00 -0.58647387E+00 -0.32835774E+00</point>
      <point r="0.860"> -0.48500150E+00 -0.42296851E+00 -0.57730004E+00 -0.33423841E+00</point>
      <point r="0.880"> -0.47634752E+00 -0.40942418E+00 -0.56726776E+00 -0.33963257E+00</point>
      <point r="0.900"> -0.46850230E+00 -0.39537137E+00 -0.55645056E+00 -0.34455445E+00</point>
      <point r="0.920"> -0.46142861E+00 -0.38087538E+00 -0.54491913E+00 -0.34901858E+00</point>
      <point r="0.940"> -0.45508739E+00 -0.36599881E+00 -0.53274138E+00 -0.35303975E+00</point>
      <point r="0.960"> -0.44943844E+00 -0.35080124E+00 -0.51998274E+00 -0.35663294E+00</point>
      <point r="0.980"> -0.44444103E+00 -0.33533901E+00 -0.50670631E+00 -0.35981321E+00</point>
      <point r="1.000"> -0.44005424E+00 -0.31966534E+00 -0.49297274E+00 -0.36259567E+00</point>
      <point r="1.020"> -0.43623719E+00 -0.30383050E+00 -0.47883993E+00 -0.36499538E+00</point>
      <point r="1.040"> -0.43294933E+00 -0.28788190E+00 -0.46436279E+00 -0.36702730E+00</point>
      <point r="1.060"> -0.43015064E+00 -0.27186416E+00 -0.44959325E+00 -0.36870627E+00</point>
      <point r="1.080"> -0.42780194E+00 -0.25581905E+00 -0.43458038E+00 -0.37004696E+00</point>
      <point r="1.100"> -0.42586515E+00 -0.23978554E+00 -0.41937061E+00 -0.37106387E+00</point>
      <point r="1.120"> -0.42430341E+00 -0.22379986E+00 -0.40400787E+00 -0.37177128E+00</point>
      <point r="1.140"> -0.42308120E+00 -0.20789555E+00 -0.38853363E+00 -0.37218317E+00</point>
      <point r="1.160"> -0.42216440E+00 -0.19210367E+00 -0.37298687E+00 -0.37231328E+00</point>
      <point r="1.180"> -0.42152034E+00 -0.17645281E+00 -0.35740408E+00 -0.37217501E+00</point>
      <point r="1.200"> -0.42111783E+00 -0.16096921E+00 -0.34181929E+00 -0.37178146E+00</point>
      <point r="1.220"> -0.42092723E+00 -0.14567677E+00 -0.32626406E+00 -0.37114540E+00</point>
      <point r="1.240"> -0.42092043E+00 -0.13059721E+00 -0.31076765E+00 -0.37027929E+00</point>
      <point r="1.260"> -0.42107087E+00 -0.11575014E+00 -0.29535707E+00 -0.36919526E+00</point>
      <point r="1.280"> -0.42135355E+00 -0.10115322E+00 -0.28005730E+00 -0.36790506E+00</point>
      <point r="1.300"> -0.42174495E+00 -0.86822325E-01 -0.26489131E+00 -0.36642007E+00</point>
      <point r="1.320"> -0.42222308E+00 -0.72771600E-01 -0.24988027E+00 -0.36475130E+00</point>
      <point r="1.340"> -0.42276743E+00 -0.59013595E-01 -0.23504353E+00 -0.36290942E+00</point>
      <point r="1.360"> -0.42335894E+00 -0.45559310E-01 -0.22039877E+00 -0.36090473E+00</point>
      <point r="1.380"> -0.42397994E+00 -0.32418256E-01 -0.20596196E+00 -0.35874718E+00</point>
      <point r="1.400"> -0.42461409E+00 -0.19598516E-01 -0.19174745E+00 -0.35644639E+00</point>
      <point r="1.420"> -0.42524635E+00 -0.71068188E-02 -0.17776803E+00 -0.35401161E+00</point>
      <point r="1.440"> -0.42586282E+00  0.50513810E-02 -0.16403500E+00 -0.35145175E+00</point>
      <point r="1.460"> -0.42645079E+00  0.16871809E-01 -0.15055828E+00 -0.34877535E+00</point>
      <point r="1.480"> -0.42699864E+00  0.28351282E-01 -0.13734652E+00 -0.34599059E+00</point>
      <point r="1.500"> -0.42749584E+00  0.39487618E-01 -0.12440715E+00 -0.34310537E+00</point>
      <point r="1.520"> -0.42793291E+00  0.50279555E-01 -0.11174654E+00 -0.34012725E+00</point>
      <point r="1.540"> -0.42830143E+00  0.60726670E-01 -0.99370034E-01 -0.33706349E+00</point>
      <point r="1.560"> -0.42859393E+00  0.70829311E-01 -0.87282032E-01 -0.33392108E+00</point>
      <point r="1.580"> -0.42880385E+00  0.80588533E-01 -0.75486051E-01 -0.33070671E+00</point>
      <point r="1.600"> -0.42892550E+00  0.90006043E-01 -0.63984773E-01 -0.32742677E+00</point>
      <point r="1.620"> -0.42895394E+00  0.99084146E-01 -0.52780092E-01 -0.32408735E+00</point>
      <point r="1.640"> -0.42888492E+00  0.10782570E+00 -0.41873153E-01 -0.32069427E+00</point>
      <point r="1.660"> -0.42870293E+00  0.11620790E+00 -0.31308214E-01 -0.31725282E+00</point>
      <point r="1.680"> -0.42843228E+00  0.12429451E+00 -0.20988447E-01 -0.31376869E+00</point>
      <point r="1.700"> -0.42805468E+00  0.13205523E+00 -0.10965758E-01 -0.31024664E+00</point>
      <point r="1.720"> -0.42756813E+00  0.13949444E+00 -0.12391277E-02 -0.30669146E+00</point>
      <point r="1.740"> -0.42697121E+00  0.14661696E+00  0.81928854E-02 -0.30310770E+00</point>
      <point r="1.760"> -0.42626294E+00  0.15342793E+00  0.17332137E-01 -0.29949969E+00</point>
      <point r="1.780"> -0.42544284E+00  0.15993277E+00  0.26180898E-01 -0.29587153E+00</point>
      <point r="1.800"> -0.42451079E+00  0.16613716E+00  0.34741861E-01 -0.29222712E+00</point>
      <point r="1.820"> -0.42346704E+00  0.17204693E+00  0.43018125E-01 -0.28857014E+00</point>
      <point r="1.840"> -0.42231218E+00  0.17766811E+00  0.51013185E-01 -0.28490408E+00</point>
      <point r="1.860"> -0.42104709E+00  0.18300681E+00  0.58730899E-01 -0.28123223E+00</point>
      <point r="1.880"> -0.41967291E+00  0.18806927E+00  0.66175452E-01 -0.27755772E+00</point>
      <point r="1.900"> -0.41819106E+00  0.19286182E+00  0.73351305E-01 -0.27388347E+00</point>
      <point r="1.920"> -0.41660317E+00  0.19739086E+00  0.80263150E-01 -0.27021225E+00</point>
      <point r="1.940"> -0.41491110E+00  0.20166287E+00  0.86915855E-01 -0.26654668E+00</point>
      <point r="1.960"> -0.41311692E+00  0.20568437E+00  0.93314422E-01 -0.26288920E+00</point>
      <point r="1.980"> -0.41122285E+00  0.20946191E+00  0.99463941E-01 -0.25924213E+00</point>
      <point r="2.000"> -0.40923129E+00  0.21300204E+00  0.10536956E+00 -0.25560762E+00</point>
      <point r="2.020"> -0.40714477E+00  0.21631132E+00  0.11103647E+00 -0.25198771E+00</point>
      <point r="2.040"> -0.40496592E+00  0.21939627E+00  0.11646987E+00 -0.24838429E+00</point>
      <point r="2.060"> -0.40269748E+00  0.22226336E+00  0.12167497E+00 -0.24479914E+00</point>
      <point r="2.080"> -0.40034224E+00  0.22491902E+00  0.12665700E+00 -0.24123391E+00</point>
      <point r="2.100"> -0.39790307E+00  0.22736958E+00  0.13142118E+00 -0.23769015E+00</point>
      <point r="2.120"> -0.39538288E+00  0.22962131E+00  0.13597274E+00 -0.23416927E+00</point>
      <point r="2.140"> -0.39278461E+00  0.23168042E+00  0.14031693E+00 -0.23067260E+00</point>
      <point r="2.160"> -0.39011125E+00  0.23355301E+00  0.14445898E+00 -0.22720137E+00</point>
      <point r="2.180"> -0.38736579E+00  0.23524510E+00  0.14840415E+00 -0.22375672E+00</point>
      <point r="2.200"> -0.38455125E+00  0.23676265E+00  0.15215768E+00 -0.22033968E+00</point>
      <point r="2.220"> -0.38167065E+00  0.23811153E+00  0.15572479E+00 -0.21695123E+00</point>
      <point r="2.240"> -0.37872704E+00  0.23929750E+00  0.15911069E+00 -0.21359223E+00</point>
      <point r="2.260"> -0.37572344E+00  0.24032628E+00  0.16232059E+00 -0.21026349E+00</point>
      <point r="2.280"> -0.37266288E+00  0.24120346E+00  0.16535961E+00 -0.20696574E+00</point>
      <point r="2.300"> -0.36954838E+00  0.24193456E+00  0.16823285E+00 -0.20369966E+00</point>
      <point r="2.320"> -0.36638292E+00  0.24252497E+00  0.17094533E+00 -0.20046582E+00</point>
      <point r="2.340"> -0.36316948E+00  0.24298000E+00  0.17350202E+00 -0.19726476E+00</point>
      <point r="2.360"> -0.35991097E+00  0.24330482E+00  0.17590776E+00 -0.19409696E+00</point>
      <point r="2.380"> -0.35661029E+00  0.24350448E+00  0.17816733E+00 -0.19096283E+00</point>
      <point r="2.400"> -0.35327026E+00  0.24358391E+00  0.18028538E+00 -0.18786274E+00</point>
      <point r="2.420"> -0.34989368E+00  0.24354790E+00  0.18226645E+00 -0.18479699E+00</point>
      <point r="2.440"> -0.34648325E+00  0.24340109E+00  0.18411496E+00 -0.18176586E+00</point>
      <point r="2.460"> -0.34304166E+00  0.24314800E+00  0.18583521E+00 -0.17876956E+00</point>
      <point r="2.480"> -0.33957149E+00  0.24279302E+00  0.18743136E+00 -0.17580828E+00</point>
      <point r="2.500"> -0.33607527E+00  0.24234037E+00  0.18890747E+00 -0.17288215E+00</point>
      <point r="2.520"> -0.33255549E+00  0.24179418E+00  0.19026744E+00 -0.16999129E+00</point>
      <point r="2.540"> -0.32901452E+00  0.24115841E+00  0.19151508E+00 -0.16713576E+00</point>
      <point r="2.560"> -0.32545471E+00  0.24043693E+00  0.19265407E+00 -0.16431560E+00</point>
      <point r="2.580"> -0.32187831E+00  0.23963347E+00  0.19368799E+00 -0.16153083E+00</point>
      <point r="2.600"> -0.31828752E+00  0.23875164E+00  0.19462028E+00 -0.15878143E+00</point>
      <point r="2.620"> -0.31468445E+00  0.23779495E+00  0.19545432E+00 -0.15606736E+00</point>
      <point r="2.640"> -0.31107119E+00  0.23676679E+00  0.19619337E+00 -0.15338853E+00</point>
      <point r="2.660"> -0.30744971E+00  0.23567046E+00  0.19684059E+00 -0.15074487E+00</point>
      <point r="2.680"> -0.30382195E+00  0.23450914E+00  0.19739906E+00 -0.14813626E+00</point>
      <point r="2.700"> -0.30018978E+00  0.23328594E+00  0.19787180E+00 -0.14556256E+00</point>
      <point r="2.720"> -0.29655501E+00  0.23200385E+00  0.19826171E+00 -0.14302363E+00</point>
      <point r="2.740"> -0.29291941E+00  0.23066579E+00  0.19857164E+00 -0.14051930E+00</point>
      <point r="2.760"> -0.28928466E+00  0.22927458E+00  0.19880436E+00 -0.13804938E+00</point>
      <point r="2.780"> -0.28565242E+00  0.22783298E+00  0.19896258E+00 -0.13561368E+00</point>
      <point r="2.800"> -0.28202427E+00  0.22634364E+00  0.19904893E+00 -0.13321199E+00</point>
      <point r="2.820"> -0.27840176E+00  0.22480913E+00  0.19906597E+00 -0.13084408E+00</point>
      <point r="2.840"> -0.27478637E+00  0.22323198E+00  0.19901621E+00 -0.12850972E+00</point>
      <point r="2.860"> -0.27117954E+00  0.22161460E+00  0.19890209E+00 -0.12620867E+00</point>
      <point r="2.880"> -0.26758266E+00  0.21995934E+00  0.19872597E+00 -0.12394067E+00</point>
      <point r="2.900"> -0.26399705E+00  0.21826848E+00  0.19849016E+00 -0.12170546E+00</point>
      <point r="2.920"> -0.26042400E+00  0.21654421E+00  0.19819688E+00 -0.11950277E+00</point>
      <point r="2.940"> -0.25686474E+00  0.21478866E+00  0.19784830E+00 -0.11733232E+00</point>
      <point r="2.960"> -0.25332044E+00  0.21300387E+00  0.19744652E+00 -0.11519383E+00</point>
      <point r="2.980"> -0.24979223E+00  0.21119182E+00  0.19699354E+00 -0.11308700E+00</point>
      <point r="3.000"> -0.24628116E+00  0.20935440E+00  0.19649132E+00 -0.11101155E+00</point>
      <point r="3.020"> -0.24278827E+00  0.20749344E+00  0.19594173E+00 -0.10896716E+00</point>
      <point r="3.040"> -0.23931451E+00  0.20561068E+00  0.19534657E+00 -0.10695353E+00</point>
      <point r="3.060"> -0.23586081E+00  0.20370782E+00  0.19470755E+00 -0.10497036E+00</point>
      <point r="3.080"> -0.23242804E+00  0.20178646E+00  0.19402636E+00 -0.10301734E+00</point>
      <point r="3.100"> -0.22901701E+00  0.19984815E+00  0.19330456E+00 -0.10109413E+00</point>
      <point r="3.120"> -0.22562851E+00  0.19789438E+00  0.19254369E+00 -0.99200444E-01</point>
      <point r="3.140"> -0.22226327E+00  0.19592655E+00  0.19174520E+00 -0.97335946E-01</point>
      <point r="3.160"> -0.21892198E+00  0.19394605E+00  0.19091050E+00 -0.95500322E-01</point>
      <point r="3.180"> -0.21560530E+00  0.19195417E+00  0.19004092E+00 -0.93693253E-01</point>
      <point r="3.200"> -0.21231383E+00  0.18995217E+00  0.18913776E+00 -0.91914418E-01</point>
      <point r="3.220"> -0.20904817E+00  0.18794124E+00  0.18820226E+00 -0.90163498E-01</point>
      <point r="3.240"> -0.20580884E+00  0.18592256E+00  0.18723560E+00 -0.88440173E-01</point>
      <point r="3.260"> -0.20259637E+00  0.18389721E+00  0.18623895E+00 -0.86744125E-01</point>
      <point r="3.280"> -0.19941122E+00  0.18186627E+00  0.18521340E+00 -0.85075037E-01</point>
      <point r="3.300"> -0.19625385E+00  0.17983077E+00  0.18416004E+00 -0.83432590E-01</point>
      <point r="3.320"> -0.19312466E+00  0.17779168E+00  0.18307991E+00 -0.81816470E-01</point>
      <point r="3.340"> -0.19002361E+00  0.17574736E+00  0.18196839E+00 -0.80226083E-01</point>
      <point r="3.360"> -0.18695138E+00  0.17370264E+00  0.18083508E+00 -0.78661691E-01</point>
      <point r="3.380"> -0.18390846E+00  0.17165721E+00  0.17967818E+00 -0.77122688E-01</point>
      <point r="3.400"> -0.18089516E+00  0.16961192E+00  0.17849866E+00 -0.75608762E-01</point>
      <point r="3.420"> -0.17791177E+00  0.16756761E+00  0.17729747E+00 -0.74119601E-01</point>
      <point r="3.440"> -0.17495853E+00  0.16552506E+00  0.17607552E+00 -0.72654898E-01</point>
      <point r="3.460"> -0.17203569E+00  0.16348505E+00  0.17483370E+00 -0.71214346E-01</point>
      <point r="3.480"> -0.16914345E+00  0.16144831E+00  0.17357287E+00 -0.69797640E-01</point>
      <point r="3.500"> -0.16628202E+00  0.15941555E+00  0.17229389E+00 -0.68404477E-01</point>
      <point r="3.520"> -0.16345155E+00  0.15738744E+00  0.17099757E+00 -0.67034556E-01</point>
      <point r="3.540"> -0.16065221E+00  0.15536463E+00  0.16968471E+00 -0.65687580E-01</point>
      <point r="3.560"> -0.15788411E+00  0.15334776E+00  0.16835610E+00 -0.64363252E-01</point>
      <point r="3.580"> -0.15514738E+00  0.15133743E+00  0.16701250E+00 -0.63061277E-01</point>
      <point r="3.600"> -0.15244211E+00  0.14933422E+00  0.16565467E+00 -0.61781365E-01</point>
      <point r="3.620"> -0.14976837E+00  0.14733869E+00  0.16428333E+00 -0.60523225E-01</point>
      <point r="3.640"> -0.14712623E+00  0.14535138E+00  0.16289921E+00 -0.59286571E-01</point>
      <point r="3.660"> -0.14451572E+00  0.14337280E+00  0.16150300E+00 -0.58071119E-01</point>
      <point r="3.680"> -0.14193689E+00  0.14140344E+00  0.16009539E+00 -0.56876584E-01</point>
      <point r="3.700"> -0.13938974E+00  0.13944379E+00  0.15867706E+00 -0.55702689E-01</point>
      <point r="3.720"> -0.13687428E+00  0.13749430E+00  0.15724865E+00 -0.54549155E-01</point>
      <point r="3.740"> -0.13439048E+00  0.13555540E+00  0.15581082E+00 -0.53415707E-01</point>
      <point r="3.760"> -0.13193832E+00  0.13362752E+00  0.15436420E+00 -0.52302072E-01</point>
      <point r="3.780"> -0.12951775E+00  0.13171104E+00  0.15290940E+00 -0.51207981E-01</point>
      <point r="3.800"> -0.12712873E+00  0.12980636E+00  0.15144703E+00 -0.50133165E-01</point>
      <point r="3.820"> -0.12477118E+00  0.12791384E+00  0.14997768E+00 -0.49077360E-01</point>
      <point r="3.840"> -0.12244502E+00  0.12603382E+00  0.14850194E+00 -0.48040301E-01</point>
      <point r="3.860"> -0.12015016E+00  0.12416663E+00  0.14702036E+00 -0.47021729E-01</point>
      <point r="3.880"> -0.11788651E+00  0.12231259E+00  0.14553351E+00 -0.46021386E-01</point>
      <point r="3.900"> -0.11565394E+00  0.12047200E+00  0.14404193E+00 -0.45039016E-01</point>
      <point r="3.920"> -0.11345233E+00  0.11864513E+00  0.14254616E+00 -0.44074366E-01</point>
      <point r="3.940"> -0.11128157E+00  0.11683227E+00  0.14104673E+00 -0.43127185E-01</point>
      <point r="3.960"> -0.10914149E+00  0.11503367E+00  0.13954414E+00 -0.42197225E-01</point>
      <point r="3.980"> -0.10703196E+00  0.11324957E+00  0.13803890E+00 -0.41284241E-01</point>
      <point r="4.000"> -0.10495281E+00  0.11148020E+00  0.13653152E+00 -0.40387988E-01</point>
      <point r="4.020"> -0.10290388E+00  0.10972578E+00  0.13502248E+00 -0.39508226E-01</point>
      <point r="4.040"> -0.10088500E+00  0.10798651E+00  0.13351225E+00 -0.38644715E-01</point>
      <point r="4.060"> -0.98895994E-01  0.10626259E+00  0.13200131E+00 -0.37797221E-01</point>
      <point r="4.080"> -0.96936668E-01  0.10455420E+00  0.13049012E+00 -0.36965508E-01</point>
      <point r="4.100"> -0.95006832E-01  0.10286152E+00  0.12897913E+00 -0.36149345E-01</point>
      <point r="4.120"> -0.93106288E-01  0.10118469E+00  0.12746878E+00 -0.35348503E-01</point>
      <point r="4.140"> -0.91234832E-01  0.99523886E-01  0.12595951E+00 -0.34562755E-01</point>
      <point r="4.160"> -0.89392254E-01  0.97879232E-01  0.12445176E+00 -0.33791876E-01</point>
      <point r="4.180"> -0.87578341E-01  0.96250865E-01  0.12294593E+00 -0.33035644E-01</point>
      <point r="4.200"> -0.85792872E-01  0.94638904E-01  0.12144244E+00 -0.32293838E-01</point>
      <point r="4.220"> -0.84035622E-01  0.93043461E-01  0.11994169E+00 -0.31566241E-01</point>
      <point r="4.240"> -0.82306363E-01  0.91464638E-01  0.11844408E+00 -0.30852637E-01</point>
      <point r="4.260"> -0.80604862E-01  0.89902528E-01  0.11694998E+00 -0.30152813E-01</point>
      <point r="4.280"> -0.78930881E-01  0.88357214E-01  0.11545978E+00 -0.29466556E-01</point>
      <point r="4.300"> -0.77284180E-01  0.86828771E-01  0.11397384E+00 -0.28793659E-01</point>
      <point r="4.320"> -0.75664512E-01  0.85317265E-01  0.11249252E+00 -0.28133915E-01</point>
      <point r="4.340"> -0.74071630E-01  0.83822754E-01  0.11101617E+00 -0.27487118E-01</point>
      <point r="4.360"> -0.72505281E-01  0.82345288E-01  0.10954513E+00 -0.26853067E-01</point>
      <point r="4.380"> -0.70965211E-01  0.80884907E-01  0.10807974E+00 -0.26231560E-01</point>
      <point r="4.400"> -0.69451163E-01  0.79441645E-01  0.10662032E+00 -0.25622400E-01</point>
      <point r="4.420"> -0.67962875E-01  0.78015527E-01  0.10516717E+00 -0.25025391E-01</point>
      <point r="4.440"> -0.66500085E-01  0.76606571E-01  0.10372062E+00 -0.24440339E-01</point>
      <point r="4.460"> -0.65062527E-01  0.75214788E-01  0.10228095E+00 -0.23867052E-01</point>
      <point r="4.480"> -0.63649934E-01  0.73840181E-01  0.10084846E+00 -0.23305341E-01</point>
      <point r="4.500"> -0.62262036E-01  0.72482748E-01  0.99423431E-01 -0.22755017E-01</point>
      <point r="4.520"> -0.60898563E-01  0.71142477E-01  0.98006128E-01 -0.22215895E-01</point>
      <point r="4.540"> -0.59559242E-01  0.69819354E-01  0.96596820E-01 -0.21687792E-01</point>
      <point r="4.560"> -0.58243799E-01  0.68513354E-01  0.95195764E-01 -0.21170527E-01</point>
      <point r="4.580"> -0.56951959E-01  0.67224450E-01  0.93803208E-01 -0.20663920E-01</point>
      <point r="4.600"> -0.55683446E-01  0.65952608E-01  0.92419393E-01 -0.20167794E-01</point>
      <point r="4.620"> -0.54437983E-01  0.64697786E-01  0.91044553E-01 -0.19681973E-01</point>
      <point r="4.640"> -0.53215293E-01  0.63459941E-01  0.89678913E-01 -0.19206285E-01</point>
      <point r="4.660"> -0.52015096E-01  0.62239020E-01  0.88322691E-01 -0.18740559E-01</point>
      <point r="4.680"> -0.50837116E-01  0.61034969E-01  0.86976098E-01 -0.18284624E-01</point>
      <point r="4.700"> -0.49681073E-01  0.59847727E-01  0.85639337E-01 -0.17838315E-01</point>
      <point r="4.720"> -0.48546689E-01  0.58677228E-01  0.84312604E-01 -0.17401465E-01</point>
      <point r="4.740"> -0.47433684E-01  0.57523404E-01  0.82996087E-01 -0.16973913E-01</point>
      <point r="4.760"> -0.46341781E-01  0.56386179E-01  0.81689967E-01 -0.16555495E-01</point>
      <point r="4.780"> -0.45270702E-01  0.55265476E-01  0.80394419E-01 -0.16146054E-01</point>
      <point r="4.800"> -0.44220168E-01  0.54161212E-01  0.79109609E-01 -0.15745431E-01</point>
      <point r="4.820"> -0.43189903E-01  0.53073300E-01  0.77835696E-01 -0.15353473E-01</point>
      <point r="4.840"> -0.42179629E-01  0.52001650E-01  0.76572834E-01 -0.14970024E-01</point>
      <point r="4.860"> -0.41189071E-01  0.50946168E-01  0.75321167E-01 -0.14594933E-01</point>
      <point r="4.880"> -0.40217953E-01  0.49906755E-01  0.74080834E-01 -0.14228052E-01</point>
      <point r="4.900"> -0.39266002E-01  0.48883312E-01  0.72851967E-01 -0.13869232E-01</point>
      <point r="4.920"> -0.38332943E-01  0.47875734E-01  0.71634689E-01 -0.13518327E-01</point>
      <point r="4.940"> -0.37418505E-01  0.46883911E-01  0.70429117E-01 -0.13175194E-01</point>
      <point r="4.960"> -0.36522417E-01  0.45907734E-01  0.69235363E-01 -0.12839691E-01</point>
      <point r="4.980"> -0.35644409E-01  0.44947089E-01  0.68053529E-01 -0.12511677E-01</point>
      <point r="5.000"> -0.34784211E-01  0.44001858E-01  0.66883713E-01 -0.12191013E-01</point>
      <point r="5.020"> -0.33941557E-01  0.43071922E-01  0.65726004E-01 -0.11877565E-01</point>
      <point r="5.040"> -0.33116181E-01  0.42157159E-01  0.64580485E-01 -0.11571196E-01</point>
      <point r="5.060"> -0.32307818E-01  0.41257442E-01  0.63447234E-01 -0.11271774E-01</point>
      <point r="5.080"> -0.31516205E-01  0.40372646E-01  0.62326319E-01 -0.10979169E-01</point>
      <point r="5.100"> -0.30741081E-01  0.39502640E-01  0.61217806E-01 -0.10693250E-01</point>
      <point r="5.120"> -0.29982187E-01  0.38647292E-01  0.60121751E-01 -0.10413890E-01</point>
      <point r="5.140"> -0.29239265E-01  0.37806467E-01  0.59038205E-01 -0.10140965E-01</point>
      <point r="5.160"> -0.28512058E-01  0.36980031E-01  0.57967212E-01 -0.98743491E-02</point>
      <point r="5.180"> -0.27800313E-01  0.36167844E-01  0.56908813E-01 -0.96139214E-02</point>
      <point r="5.200"> -0.27103777E-01  0.35369766E-01  0.55863039E-01 -0.93595611E-02</point>
      <point r="5.220"> -0.26422199E-01  0.34585658E-01  0.54829917E-01 -0.91111498E-02</point>
      <point r="5.240"> -0.25755332E-01  0.33815374E-01  0.53809469E-01 -0.88685706E-02</point>
      <point r="5.260"> -0.25102929E-01  0.33058772E-01  0.52801711E-01 -0.86317082E-02</point>
      <point r="5.280"> -0.24464746E-01  0.32315704E-01  0.51806651E-01 -0.84004493E-02</point>
      <point r="5.300"> -0.23840541E-01  0.31586025E-01  0.50824296E-01 -0.81746820E-02</point>
      <point r="5.320"> -0.23230073E-01  0.30869585E-01  0.49854644E-01 -0.79542963E-02</point>
      <point r="5.340"> -0.22633106E-01  0.30166236E-01  0.48897690E-01 -0.77391837E-02</point>
      <point r="5.360"> -0.22049404E-01  0.29475827E-01  0.47953423E-01 -0.75292376E-02</point>
      <point r="5.380"> -0.21478734E-01  0.28798208E-01  0.47021828E-01 -0.73243527E-02</point>
      <point r="5.400"> -0.20920865E-01  0.28133227E-01  0.46102885E-01 -0.71244256E-02</point>
      <point r="5.420"> -0.20375569E-01  0.27480732E-01  0.45196567E-01 -0.69293545E-02</point>
      <point r="5.440"> -0.19842620E-01  0.26840569E-01  0.44302847E-01 -0.67390390E-02</point>
      <point r="5.460"> -0.19321793E-01  0.26212586E-01  0.43421690E-01 -0.65533804E-02</point>
      <point r="5.480"> -0.18812869E-01  0.25596628E-01  0.42553058E-01 -0.63722818E-02</point>
      <point r="5.500"> -0.18315628E-01  0.24992542E-01  0.41696909E-01 -0.61956475E-02</point>
      <point r="5.520"> -0.17829854E-01  0.24400173E-01  0.40853197E-01 -0.60233835E-02</point>
      <point r="5.540"> -0.17355333E-01  0.23819368E-01  0.40021872E-01 -0.58553974E-02</point>
      <point r="5.560"> -0.16891854E-01  0.23249972E-01  0.39202880E-01 -0.56915983E-02</point>
      <point r="5.580"> -0.16439209E-01  0.22691829E-01  0.38396163E-01 -0.55318968E-02</point>
      <point r="5.600"> -0.15997192E-01  0.22144787E-01  0.37601662E-01 -0.53762049E-02</point>
      <point r="5.620"> -0.15565598E-01  0.21608691E-01  0.36819310E-01 -0.52244362E-02</point>
      <point r="5.640"> -0.15144227E-01  0.21083387E-01  0.36049041E-01 -0.50765057E-02</point>
      <point r="5.660"> -0.14732880E-01  0.20568721E-01  0.35290784E-01 -0.49323299E-02</point>
      <point r="5.680"> -0.14331362E-01  0.20064540E-01  0.34544466E-01 -0.47918267E-02</point>
      <point r="5.700"> -0.13939478E-01  0.19570692E-01  0.33810008E-01 -0.46549154E-02</point>
      <point r="5.720"> -0.13557040E-01  0.19087024E-01  0.33087333E-01 -0.45215167E-02</point>
      <point r="5.740"> -0.13183857E-01  0.18613384E-01  0.32376357E-01 -0.43915528E-02</point>
      <point r="5.760"> -0.12819746E-01  0.18149620E-01  0.31676995E-01 -0.42649472E-02</point>
      <point r="5.780"> -0.12464522E-01  0.17695583E-01  0.30989161E-01 -0.41416247E-02</point>
      <point r="5.800"> -0.12118005E-01  0.17251122E-01  0.30312765E-01 -0.40215115E-02</point>
      <point r="5.820"> -0.11780019E-01  0.16816089E-01  0.29647714E-01 -0.39045353E-02</point>
      <point r="5.840"> -0.11450386E-01  0.16390334E-01  0.28993914E-01 -0.37906248E-02</point>
      <point r="5.860"> -0.11128936E-01  0.15973710E-01  0.28351270E-01 -0.36797101E-02</point>
      <point r="5.880"> -0.10815496E-01  0.15566071E-01  0.27719683E-01 -0.35717229E-02</point>
      <point r="5.900"> -0.10509901E-01  0.15167270E-01  0.27099053E-01 -0.34665956E-02</point>
      <point r="5.920"> -0.10211985E-01  0.14777163E-01  0.26489278E-01 -0.33642624E-02</point>
      <point r="5.940"> -0.99215860E-02  0.14395607E-01  0.25890255E-01 -0.32646584E-02</point>
      <point r="5.960"> -0.96385431E-02  0.14022458E-01  0.25301878E-01 -0.31677201E-02</point>
      <point r="5.980"> -0.93626994E-02  0.13657574E-01  0.24724042E-01 -0.30733851E-02</point>
      <point r="6.000"> -0.90938999E-02  0.13300816E-01  0.24156638E-01 -0.29815921E-02</point>
      <point r="6.020"> -0.88319920E-02  0.12952043E-01  0.23599558E-01 -0.28922814E-02</point>
      <point r="6.040"> -0.85768258E-02  0.12611117E-01  0.23052690E-01 -0.28053939E-02</point>
      <point r="6.060"> -0.83282537E-02  0.12277902E-01  0.22515923E-01 -0.27208720E-02</point>
      <point r="6.080"> -0.80861307E-02  0.11952260E-01  0.21989145E-01 -0.26386592E-02</point>
      <point r="6.100"> -0.78503139E-02  0.11634059E-01  0.21472242E-01 -0.25587000E-02</point>
      <point r="6.120"> -0.76206631E-02  0.11323164E-01  0.20965101E-01 -0.24809402E-02</point>
      <point r="6.140"> -0.73970404E-02  0.11019443E-01  0.20467605E-01 -0.24053264E-02</point>
      <point r="6.160"> -0.71793102E-02  0.10722766E-01  0.19979639E-01 -0.23318064E-02</point>
      <point r="6.180"> -0.69673394E-02  0.10433002E-01  0.19501087E-01 -0.22603292E-02</point>
      <point r="6.200"> -0.67609970E-02  0.10150025E-01  0.19031832E-01 -0.21908447E-02</point>
      <point r="6.220"> -0.65601545E-02  0.98737074E-02  0.18571755E-01 -0.21233038E-02</point>
      <point r="6.240"> -0.63646856E-02  0.96039237E-02  0.18120741E-01 -0.20576584E-02</point>
      <point r="6.260"> -0.61744662E-02  0.93405504E-02  0.17678669E-01 -0.19938616E-02</point>
      <point r="6.280"> -0.59893747E-02  0.90834648E-02  0.17245422E-01 -0.19318672E-02</point>
      <point r="6.300"> -0.58092914E-02  0.88325461E-02  0.16820881E-01 -0.18716301E-02</point>
      <point r="6.320"> -0.56340991E-02  0.85876746E-02  0.16404928E-01 -0.18131061E-02</point>
      <point r="6.340"> -0.54636824E-02  0.83487322E-02  0.15997443E-01 -0.17562521E-02</point>
      <point r="6.360"> -0.52979285E-02  0.81156023E-02  0.15598308E-01 -0.17010258E-02</point>
      <point r="6.380"> -0.51367264E-02  0.78881697E-02  0.15207404E-01 -0.16473856E-02</point>
      <point r="6.400"> -0.49799674E-02  0.76663208E-02  0.14824613E-01 -0.15952912E-02</point>
      <point r="6.420"> -0.48275446E-02  0.74499433E-02  0.14449815E-01 -0.15447029E-02</point>
      <point r="6.440"> -0.46793536E-02  0.72389266E-02  0.14082893E-01 -0.14955820E-02</point>
      <point r="6.460"> -0.45352917E-02  0.70331613E-02  0.13723730E-01 -0.14478904E-02</point>
      <point r="6.480"> -0.43952583E-02  0.68325399E-02  0.13372206E-01 -0.14015912E-02</point>
      <point r="6.500"> -0.42591548E-02  0.66369562E-02  0.13028206E-01 -0.13566481E-02</point>
      <point r="6.520"> -0.41268846E-02  0.64463053E-02  0.12691611E-01 -0.13130256E-02</point>
      <point r="6.540"> -0.39983529E-02  0.62604843E-02  0.12362308E-01 -0.12706890E-02</point>
      <point r="6.560"> -0.38734670E-02  0.60793913E-02  0.12040178E-01 -0.12296045E-02</point>
      <point r="6.580"> -0.37521360E-02  0.59029262E-02  0.11725108E-01 -0.11897390E-02</point>
      <point r="6.600"> -0.36342707E-02  0.57309904E-02  0.11416982E-01 -0.11510601E-02</point>
      <point r="6.620"> -0.35197841E-02  0.55634866E-02  0.11115688E-01 -0.11135361E-02</point>
      <point r="6.640"> -0.34085907E-02  0.54003192E-02  0.10821111E-01 -0.10771362E-02</point>
      <point r="6.660"> -0.33006068E-02  0.52413941E-02  0.10533139E-01 -0.10418302E-02</point>
      <point r="6.680"> -0.31933020E-02  0.50820636E-02  0.10240680E-01 -0.10071522E-02</point>
      <point r="6.700"> -0.30915104E-02  0.49313732E-02  0.99656210E-02 -0.97395531E-03</point>
      <point r="6.720"> -0.29926912E-02  0.47846577E-02  0.96968497E-02 -0.94176639E-03</point>
      <point r="6.740"> -0.28967678E-02  0.46418287E-02  0.94342574E-02 -0.91055799E-03</point>
      <point r="6.760"> -0.28036651E-02  0.45027995E-02  0.91777358E-02 -0.88030326E-03</point>
      <point r="6.780"> -0.27133096E-02  0.43674847E-02  0.89271777E-02 -0.85097598E-03</point>
      <point r="6.800"> -0.26256294E-02  0.42358004E-02  0.86824766E-02 -0.82255057E-03</point>
      <point r="6.820"> -0.25405543E-02  0.41076639E-02  0.84435271E-02 -0.79500207E-03</point>
      <point r="6.840"> -0.24580154E-02  0.39829942E-02  0.82102246E-02 -0.76830608E-03</point>
      <point r="6.860"> -0.23779456E-02  0.38617116E-02  0.79824655E-02 -0.74243883E-03</point>
      <point r="6.880"> -0.23002791E-02  0.37437378E-02  0.77601473E-02 -0.71737711E-03</point>
      <point r="6.900"> -0.22249516E-02  0.36289957E-02  0.75431683E-02 -0.69309828E-03</point>
      <point r="6.920"> -0.21519003E-02  0.35174097E-02  0.73314281E-02 -0.66958023E-03</point>
      <point r="6.940"> -0.20810639E-02  0.34089058E-02  0.71248270E-02 -0.64680144E-03</point>
      <point r="6.960"> -0.20123823E-02  0.33034108E-02  0.69232666E-02 -0.62474087E-03</point>
      <point r="6.980"> -0.19457969E-02  0.32008533E-02  0.67266495E-02 -0.60337804E-03</point>
      <point r="7.000"> -0.18812505E-02  0.31011631E-02  0.65348795E-02 -0.58269297E-03</point>
      <point r="7.020"> -0.18186870E-02  0.30042710E-02  0.63478612E-02 -0.56266617E-03</point>
      <point r="7.040"> -0.17580519E-02  0.29101096E-02  0.61655007E-02 -0.54327865E-03</point>
      <point r="7.060"> -0.16992918E-02  0.28186124E-02  0.59877048E-02 -0.52451189E-03</point>
      <point r="7.080"> -0.16423546E-02  0.27297143E-02  0.58143819E-02 -0.50634787E-03</point>
      <point r="7.100"> -0.15871894E-02  0.26433515E-02  0.56454410E-02 -0.48876898E-03</point>
      <point r="7.120"> -0.15337465E-02  0.25594612E-02  0.54807927E-02 -0.47175811E-03</point>
      <point r="7.140"> -0.14819776E-02  0.24779822E-02  0.53203486E-02 -0.45529857E-03</point>
      <point r="7.160"> -0.14318352E-02  0.23988542E-02  0.51640213E-02 -0.43937408E-03</point>
      <point r="7.180"> -0.13832733E-02  0.23220182E-02  0.50117248E-02 -0.42396883E-03</point>
      <point r="7.200"> -0.13362468E-02  0.22474164E-02  0.48633741E-02 -0.40906737E-03</point>
      <point r="7.220"> -0.12907117E-02  0.21749922E-02  0.47188856E-02 -0.39465470E-03</point>
      <point r="7.240"> -0.12466253E-02  0.21046900E-02  0.45781765E-02 -0.38071620E-03</point>
      <point r="7.260"> -0.12039457E-02  0.20364555E-02  0.44411656E-02 -0.36723761E-03</point>
      <point r="7.280"> -0.11626322E-02  0.19702355E-02  0.43077727E-02 -0.35420509E-03</point>
      <point r="7.300"> -0.11226450E-02  0.19059779E-02  0.41779186E-02 -0.34160515E-03</point>
      <point r="7.320"> -0.10839454E-02  0.18436316E-02  0.40515256E-02 -0.32942465E-03</point>
      <point r="7.340"> -0.10464958E-02  0.17831466E-02  0.39285170E-02 -0.31765083E-03</point>
      <point r="7.360"> -0.10102592E-02  0.17244742E-02  0.38088174E-02 -0.30627125E-03</point>
      <point r="7.380"> -0.97519987E-03  0.16675665E-02  0.36923524E-02 -0.29527382E-03</point>
      <point r="7.400"> -0.94128284E-03  0.16123766E-02  0.35790490E-02 -0.28464679E-03</point>
      <point r="7.420"> -0.90847410E-03  0.15588589E-02  0.34688352E-02 -0.27437870E-03</point>
      <point r="7.440"> -0.87674050E-03  0.15069685E-02  0.33616404E-02 -0.26445844E-03</point>
      <point r="7.460"> -0.84604977E-03  0.14566616E-02  0.32573949E-02 -0.25487518E-03</point>
      <point r="7.480"> -0.81637049E-03  0.14078956E-02  0.31560304E-02 -0.24561842E-03</point>
      <point r="7.500"> -0.78767205E-03  0.13606286E-02  0.30574796E-02 -0.23667792E-03</point>
      <point r="7.520"> -0.75992468E-03  0.13148196E-02  0.29616765E-02 -0.22804375E-03</point>
      <point r="7.540"> -0.73309937E-03  0.12704287E-02  0.28685562E-02 -0.21970624E-03</point>
      <point r="7.560"> -0.70716793E-03  0.12274169E-02  0.27780550E-02 -0.21165602E-03</point>
      <point r="7.580"> -0.68210290E-03  0.11857461E-02  0.26901101E-02 -0.20388395E-03</point>
      <point r="7.600"> -0.65787757E-03  0.11453791E-02  0.26046603E-02 -0.19638118E-03</point>
      <point r="7.620"> -0.63446597E-03  0.11062793E-02  0.25216450E-02 -0.18913910E-03</point>
      <point r="7.640"> -0.61184283E-03  0.10684115E-02  0.24410051E-02 -0.18214935E-03</point>
      <point r="7.660"> -0.58998357E-03  0.10317409E-02  0.23626826E-02 -0.17540381E-03</point>
      <point r="7.680"> -0.56886432E-03  0.99623359E-03  0.22866204E-02 -0.16889459E-03</point>
      <point r="7.700"> -0.54846183E-03  0.96185668E-03  0.22127626E-02 -0.16261403E-03</point>
      <point r="7.720"> -0.52875354E-03  0.92857792E-03  0.21410544E-02 -0.15655471E-03</point>
      <point r="7.740"> -0.50971751E-03  0.89636588E-03  0.20714420E-02 -0.15070940E-03</point>
      <point r="7.760"> -0.49133240E-03  0.86518992E-03  0.20038728E-02 -0.14507111E-03</point>
      <point r="7.780"> -0.47357751E-03  0.83502013E-03  0.19382952E-02 -0.13963302E-03</point>
      <point r="7.800"> -0.45643270E-03  0.80582736E-03  0.18746586E-02 -0.13438856E-03</point>
      <point r="7.820"> -0.43987843E-03  0.77758319E-03  0.18129135E-02 -0.12933132E-03</point>
      <point r="7.840"> -0.42389570E-03  0.75025993E-03  0.17530113E-02 -0.12445508E-03</point>
      <point r="7.860"> -0.40846608E-03  0.72383055E-03  0.16949044E-02 -0.11975384E-03</point>
      <point r="7.880"> -0.39357167E-03  0.69826875E-03  0.16385465E-02 -0.11522174E-03</point>
      <point r="7.900"> -0.37919507E-03  0.67354888E-03  0.15838919E-02 -0.11085312E-03</point>
      <point r="7.920"> -0.36531942E-03  0.64964595E-03  0.15308961E-02 -0.10664248E-03</point>
      <point r="7.940"> -0.35192835E-03  0.62653563E-03  0.14795154E-02 -0.10258450E-03</point>
      <point r="7.960"> -0.33900595E-03  0.60419421E-03  0.14297071E-02 -0.98674024E-04</point>
      <point r="7.980"> -0.32653681E-03  0.58259860E-03  0.13814296E-02 -0.94906034E-04</point>
      <point r="8.000"> -0.31450597E-03  0.56172633E-03  0.13346419E-02 -0.91275680E-04</point>
      <point r="8.020"> -0.30289891E-03  0.54155549E-03  0.12893041E-02 -0.87778263E-04</point>
      <point r="8.040"> -0.29170155E-03  0.52206479E-03  0.12453772E-02 -0.84409224E-04</point>
      <point r="8.060"> -0.28090023E-03  0.50323347E-03  0.12028229E-02 -0.81164150E-04</point>
      <point r="8.080"> -0.27048173E-03  0.48504136E-03  0.11616039E-02 -0.78038763E-04</point>
      <point r="8.100"> -0.26043319E-03  0.46746879E-03  0.11216838E-02 -0.75028920E-04</point>
      <point r="8.120"> -0.25074217E-03  0.45049667E-03  0.10830269E-02 -0.72130609E-04</point>
      <point r="8.140"> -0.24139660E-03  0.43410638E-03  0.10455982E-02 -0.69339943E-04</point>
      <point r="8.160"> -0.23238479E-03  0.41827982E-03  0.10093639E-02 -0.66653160E-04</point>
      <point r="8.180"> -0.22369540E-03  0.40299941E-03  0.97429047E-03 -0.64066615E-04</point>
      <point r="8.200"> -0.21531744E-03  0.38824802E-03  0.94034557E-03 -0.61576782E-04</point>
      <point r="8.220"> -0.20724027E-03  0.37400899E-03  0.90749740E-03 -0.59180247E-04</point>
      <point r="8.240"> -0.19945359E-03  0.36026614E-03  0.87571498E-03 -0.56873704E-04</point>
      <point r="8.260"> -0.19194739E-03  0.34700373E-03  0.84496802E-03 -0.54653956E-04</point>
      <point r="8.280"> -0.18471200E-03  0.33420646E-03  0.81522697E-03 -0.52517908E-04</point>
      <point r="8.300"> -0.17773805E-03  0.32185944E-03  0.78646296E-03 -0.50462566E-04</point>
      <point r="8.320"> -0.17101647E-03  0.30994823E-03  0.75864784E-03 -0.48485035E-04</point>
      <point r="8.340"> -0.16453846E-03  0.29845876E-03  0.73175409E-03 -0.46582512E-04</point>
      <point r="8.360"> -0.15829553E-03  0.28737738E-03  0.70575490E-03 -0.44752288E-04</point>
      <point r="8.380"> -0.15227943E-03  0.27669083E-03  0.68062408E-03 -0.42991744E-04</point>
      <point r="8.400"> -0.14648219E-03  0.26638621E-03  0.65633608E-03 -0.41298346E-04</point>
      <point r="8.420"> -0.14089610E-03  0.25645101E-03  0.63286598E-03 -0.39669645E-04</point>
      <point r="8.440"> -0.13551369E-03  0.24687305E-03  0.61018946E-03 -0.38103274E-04</point>
      <point r="8.460"> -0.13032774E-03  0.23764054E-03  0.58828283E-03 -0.36596944E-04</point>
      <point r="8.480"> -0.12533126E-03  0.22874199E-03  0.56712295E-03 -0.35148445E-04</point>
      <point r="8.500"> -0.12051749E-03  0.22016628E-03  0.54668729E-03 -0.33755640E-04</point>
      <point r="8.520"> -0.11587990E-03  0.21190259E-03  0.52695386E-03 -0.32416465E-04</point>
      <point r="8.540"> -0.11141215E-03  0.20394043E-03  0.50790123E-03 -0.31128927E-04</point>
      <point r="8.560"> -0.10710814E-03  0.19626961E-03  0.48950854E-03 -0.29891099E-04</point>
      <point r="8.580"> -0.10296196E-03  0.18888025E-03  0.47175542E-03 -0.28701122E-04</point>
      <point r="8.600"> -0.98967883E-04  0.18176276E-03  0.45462206E-03 -0.27557200E-04</point>
      <point r="8.620"> -0.95120402E-04  0.17490783E-03  0.43808913E-03 -0.26457600E-04</point>
      <point r="8.640"> -0.91414175E-04  0.16830645E-03  0.42213782E-03 -0.25400648E-04</point>
      <point r="8.660"> -0.87844046E-04  0.16194985E-03  0.40674982E-03 -0.24384730E-04</point>
      <point r="8.680"> -0.84405034E-04  0.15582956E-03  0.39190727E-03 -0.23408287E-04</point>
      <point r="8.700"> -0.81092329E-04  0.14993734E-03  0.37759280E-03 -0.22469814E-04</point>
      <point r="8.720"> -0.77901281E-04  0.14426522E-03  0.36378951E-03 -0.21567862E-04</point>
      <point r="8.740"> -0.74827405E-04  0.13880545E-03  0.35048094E-03 -0.20701031E-04</point>
      <point r="8.760"> -0.71866369E-04  0.13355056E-03  0.33765108E-03 -0.19867971E-04</point>
      <point r="8.780"> -0.69013989E-04  0.12849328E-03  0.32528434E-03 -0.19067382E-04</point>
      <point r="8.800"> -0.66266229E-04  0.12362657E-03  0.31336557E-03 -0.18298008E-04</point>
      <point r="8.820"> -0.63619194E-04  0.11894362E-03  0.30188003E-03 -0.17558641E-04</point>
      <point r="8.840"> -0.61069125E-04  0.11443783E-03  0.29081339E-03 -0.16848114E-04</point>
      <point r="8.860"> -0.58612397E-04  0.11010282E-03  0.28015171E-03 -0.16165304E-04</point>
      <point r="8.880"> -0.56245513E-04  0.10593240E-03  0.26988145E-03 -0.15509129E-04</point>
      <point r="8.900"> -0.53965099E-04  0.10192058E-03  0.25998944E-03 -0.14878547E-04</point>
      <point r="8.920"> -0.51767904E-04  0.98061565E-04  0.25046290E-03 -0.14272553E-04</point>
      <point r="8.940"> -0.49650793E-04  0.94349750E-04  0.24128939E-03 -0.13690179E-04</point>
      <point r="8.960"> -0.47610744E-04  0.90779710E-04  0.23245686E-03 -0.13130495E-04</point>
      <point r="8.980"> -0.45644845E-04  0.87346200E-04  0.22395357E-03 -0.12592604E-04</point>
      <point r="9.000"> -0.43750290E-04  0.84044145E-04  0.21576816E-03 -0.12075643E-04</point>
      <point r="9.020"> -0.41924377E-04  0.80868640E-04  0.20788957E-03 -0.11578781E-04</point>
      <point r="9.040"> -0.40164501E-04  0.77814943E-04  0.20030708E-03 -0.11101220E-04</point>
      <point r="9.060"> -0.38468155E-04  0.74878470E-04  0.19301031E-03 -0.10642189E-04</point>
      <point r="9.080"> -0.36832926E-04  0.72054791E-04  0.18598914E-03 -0.10200949E-04</point>
      <point r="9.100"> -0.35256488E-04  0.69339628E-04  0.17923381E-03 -0.97767900E-05</point>
      <point r="9.120"> -0.33736605E-04  0.66728847E-04  0.17273481E-03 -0.93690264E-05</point>
      <point r="9.140"> -0.32271124E-04  0.64218456E-04  0.16648294E-03 -0.89770009E-05</point>
      <point r="9.160"> -0.30857974E-04  0.61804599E-04  0.16046930E-03 -0.86000813E-05</point>
      <point r="9.180"> -0.29495162E-04  0.59483557E-04  0.15468522E-03 -0.82376598E-05</point>
      <point r="9.200"> -0.28180771E-04  0.57251738E-04  0.14912235E-03 -0.78891525E-05</point>
      <point r="9.220"> -0.26912960E-04  0.55105675E-04  0.14377257E-03 -0.75539983E-05</point>
      <point r="9.240"> -0.25689955E-04  0.53042025E-04  0.13862801E-03 -0.72316582E-05</point>
      <point r="9.260"> -0.24510053E-04  0.51057562E-04  0.13368108E-03 -0.69216144E-05</point>
      <point r="9.280"> -0.23371617E-04  0.49149176E-04  0.12892441E-03 -0.66233700E-05</point>
      <point r="9.300"> -0.22273074E-04  0.47313868E-04  0.12435087E-03 -0.63364475E-05</point>
      <point r="9.320"> -0.21212911E-04  0.45548747E-04  0.11995357E-03 -0.60603888E-05</point>
      <point r="9.340"> -0.20189677E-04  0.43851026E-04  0.11572583E-03 -0.57947540E-05</point>
      <point r="9.360"> -0.19201976E-04  0.42218020E-04  0.11166120E-03 -0.55391211E-05</point>
      <point r="9.380"> -0.18248469E-04  0.40647144E-04  0.10775344E-03 -0.52930852E-05</point>
      <point r="9.400"> -0.17327870E-04  0.39135904E-04  0.10399651E-03 -0.50562577E-05</point>
      <point r="9.420"> -0.16438945E-04  0.37681904E-04  0.10038459E-03 -0.48282660E-05</point>
      <point r="9.440"> -0.15580508E-04  0.36282832E-04  0.96912033E-04 -0.46087527E-05</point>
      <point r="9.460"> -0.14751422E-04  0.34936467E-04  0.93573405E-04 -0.43973754E-05</point>
      <point r="9.480"> -0.13950596E-04  0.33640669E-04  0.90363444E-04 -0.41938054E-05</point>
      <point r="9.500"> -0.13176982E-04  0.32393380E-04  0.87277072E-04 -0.39977279E-05</point>
      <point r="9.520"> -0.12429578E-04  0.31192619E-04  0.84309386E-04 -0.38088413E-05</point>
      <point r="9.540"> -0.11707418E-04  0.30036485E-04  0.81455655E-04 -0.36268563E-05</point>
      <point r="9.560"> -0.11009581E-04  0.28923145E-04  0.78711312E-04 -0.34514959E-05</point>
      <point r="9.580"> -0.10335180E-04  0.27850842E-04  0.76071953E-04 -0.32824949E-05</point>
      <point r="9.600"> -0.96833666E-05  0.26817884E-04  0.73533328E-04 -0.31195991E-05</point>
      <point r="9.620"> -0.90533269E-05  0.25822647E-04  0.71091342E-04 -0.29625652E-05</point>
      <point r="9.640"> -0.84442811E-05  0.24863571E-04  0.68742044E-04 -0.28111601E-05</point>
      <point r="9.660"> -0.78554818E-05  0.23939157E-04  0.66481627E-04 -0.26651609E-05</point>
      <point r="9.680"> -0.72862129E-05  0.23047968E-04  0.64306425E-04 -0.25243540E-05</point>
      <point r="9.700"> -0.67357884E-05  0.22188623E-04  0.62212904E-04 -0.23885352E-05</point>
      <point r="9.720"> -0.62035512E-05  0.21359796E-04  0.60197662E-04 -0.22575090E-05</point>
      <point r="9.740"> -0.56888716E-05  0.20560218E-04  0.58257423E-04 -0.21310883E-05</point>
      <point r="9.760"> -0.51911468E-05  0.19788668E-04  0.56389034E-04 -0.20090942E-05</point>
      <point r="9.780"> -0.47097993E-05  0.19043979E-04  0.54589461E-04 -0.18913556E-05</point>
      <point r="9.800"> -0.42442763E-05  0.18325029E-04  0.52855786E-04 -0.17777088E-05</point>
      <point r="9.820"> -0.37940482E-05  0.17630745E-04  0.51185201E-04 -0.16679974E-05</point>
      <point r="9.840"> -0.33586082E-05  0.16960099E-04  0.49575007E-04 -0.15620716E-05</point>
      <point r="9.860"> -0.29374709E-05  0.16312104E-04  0.48022611E-04 -0.14597884E-05</point>
      <point r="9.880"> -0.25301716E-05  0.15685817E-04  0.46525520E-04 -0.13610111E-05</point>
      <point r="9.900"> -0.21362654E-05  0.15080335E-04  0.45081339E-04 -0.12656090E-05</point>
      <point r="9.920"> -0.17553265E-05  0.14494793E-04  0.43687770E-04 -0.11734572E-05</point>
      <point r="9.940"> -0.13869471E-05  0.13928364E-04  0.42342604E-04 -0.10844362E-05</point>
      <point r="9.960"> -0.10307370E-05  0.13380256E-04  0.41043724E-04 -0.99843198E-06</point>
      <point r="9.980"> -0.68632232E-06  0.12849713E-04  0.39789095E-04 -0.91533559E-06</point>
      <point r="10.000"> -0.35334542E-06  0.12336012E-04  0.38576770E-04 -0.83504286E-06</point>
      <point r="10.020"> -0.31463686E-07  0.11838460E-04  0.37404879E-04 -0.75745428E-06</point>
      <point r="10.040">  0.27965087E-06  0.11356397E-04  0.36271629E-04 -0.68247479E-06</point>
      <point r="10.060">  0.58031235E-06  0.10889192E-04  0.35175305E-04 -0.61001360E-06</point>
      <point r="10.080">  0.87082160E-06  0.10436244E-04  0.34114262E-04 -0.53998393E-06</point>
      <point r="10.100">  0.11514668E-05  0.99969752E-05  0.33086925E-04 -0.47230291E-06</point>
      <point r="10.120">  0.14225238E-05  0.95708384E-05  0.32091789E-04 -0.40689137E-06</point>
      <point r="10.140">  0.16842574E-05  0.91573094E-05  0.31127411E-04 -0.34367366E-06</point>
      <point r="10.160">  0.19369208E-05  0.87558887E-05  0.30192411E-04 -0.28257752E-06</point>
      <point r="10.180">  0.21807570E-05  0.83661001E-05  0.29285473E-04 -0.22353390E-06</point>
      <point r="10.200">  0.24159990E-05  0.79874896E-05  0.28405335E-04 -0.16647682E-06</point>
      <point r="10.220">  0.26428701E-05  0.76196247E-05  0.27550794E-04 -0.11134326E-06</point>
      <point r="10.240">  0.28615847E-05  0.72620933E-05  0.26720701E-04 -0.58072946E-07</point>
      <point r="10.260">  0.30723484E-05  0.69145031E-05  0.25913958E-04 -0.66082972E-08</point>
      <point r="10.280">  0.32753585E-05  0.65764803E-05  0.25129519E-04  0.43105738E-07</point>
      <point r="10.300">  0.34708046E-05  0.62476696E-05  0.24366385E-04  0.91121797E-07</point>
      <point r="10.320">  0.36588688E-05  0.59277325E-05  0.23623605E-04  0.13749021E-06</point>
      <point r="10.340">  0.38397261E-05  0.56163475E-05  0.22900272E-04  0.18225914E-06</point>
      <point r="10.360">  0.40135447E-05  0.53132085E-05  0.22195521E-04  0.22547463E-06</point>
      <point r="10.380">  0.41804866E-05  0.50180247E-05  0.21508532E-04  0.26718077E-06</point>
      <point r="10.400">  0.00000000E+00</point>
    </H_spline>
    <S_spline>
      <point r="0.000">  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.44000000E+01</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400">  0.79377352E+00  0.46598497E-01  0.70399492E+00  0.71171912E+00</point>
      <point r="0.420">  0.79065131E+00  0.47252568E-01  0.70164744E+00  0.71129666E+00</point>
      <point r="0.440">  0.78737534E+00  0.47656330E-01  0.69897853E+00  0.71079864E+00</point>
      <point r="0.460">  0.78396009E+00  0.47802761E-01  0.69597616E+00  0.71022068E+00</point>
      <point r="0.480">  0.78042009E+00  0.47686087E-01  0.69262985E+00  0.70955866E+00</point>
      <point r="0.500">  0.77676972E+00  0.47301774E-01  0.68893061E+00  0.70880871E+00</point>
      <point r="0.520">  0.77302306E+00  0.46646500E-01  0.68487096E+00  0.70796726E+00</point>
      <point r="0.540">  0.76919373E+00  0.45718120E-01  0.68044485E+00  0.70703102E+00</point>
      <point r="0.560">  0.76529482E+00  0.44515626E-01  0.67564764E+00  0.70599696E+00</point>
      <point r="0.580">  0.76133874E+00  0.43039096E-01  0.67047609E+00  0.70486235E+00</point>
      <point r="0.600">  0.75733721E+00  0.41289646E-01  0.66492823E+00  0.70362472E+00</point>
      <point r="0.620">  0.75330118E+00  0.39269367E-01  0.65900336E+00  0.70228190E+00</point>
      <point r="0.640">  0.74924081E+00  0.36981275E-01  0.65270198E+00  0.70083196E+00</point>
      <point r="0.660">  0.74516544E+00  0.34429241E-01  0.64602574E+00  0.69927327E+00</point>
      <point r="0.680">  0.74108360E+00  0.31617940E-01  0.63897734E+00  0.69760442E+00</point>
      <point r="0.700">  0.73700299E+00  0.28552784E-01  0.63156051E+00  0.69582428E+00</point>
      <point r="0.720">  0.73293054E+00  0.25239864E-01  0.62377994E+00  0.69393194E+00</point>
      <point r="0.740">  0.72887235E+00  0.21685891E-01  0.61564120E+00  0.69192676E+00</point>
      <point r="0.760">  0.72483380E+00  0.17898139E-01  0.60715069E+00  0.68980829E+00</point>
      <point r="0.780">  0.72081950E+00  0.13884385E-01  0.59831560E+00  0.68757632E+00</point>
      <point r="0.800">  0.71683340E+00  0.96528614E-02  0.58914382E+00  0.68523085E+00</point>
      <point r="0.820">  0.71287877E+00  0.52121982E-02  0.57964391E+00  0.68277207E+00</point>
      <point r="0.840">  0.70895824E+00  0.57137562E-03  0.56982503E+00  0.68020037E+00</point>
      <point r="0.860">  0.70507387E+00 -0.42603244E-02  0.55969690E+00  0.67751632E+00</point>
      <point r="0.880">  0.70122717E+00 -0.92733632E-02  0.54926972E+00  0.67472067E+00</point>
      <point r="0.900">  0.69741913E+00 -0.14457989E-01  0.53855416E+00  0.67181434E+00</point>
      <point r="0.920">  0.69365028E+00 -0.19804276E-01  0.52756130E+00  0.66879838E+00</point>
      <point r="0.940">  0.68992072E+00 -0.25302167E-01  0.51630258E+00  0.66567402E+00</point>
      <point r="0.960">  0.68623013E+00 -0.30941503E-01  0.50478973E+00  0.66244262E+00</point>
      <point r="0.980">  0.68257788E+00 -0.36712063E-01  0.49303479E+00  0.65910567E+00</point>
      <point r="1.000">  0.67896296E+00 -0.42603594E-01  0.48105001E+00  0.65566478E+00</point>
      <point r="1.020">  0.67538413E+00 -0.48605838E-01  0.46884786E+00  0.65212168E+00</point>
      <point r="1.040">  0.67183985E+00 -0.54708565E-01  0.45644096E+00  0.64847820E+00</point>
      <point r="1.060">  0.66832838E+00 -0.60901593E-01  0.44384205E+00  0.64473629E+00</point>
      <point r="1.080">  0.66484777E+00 -0.67174815E-01  0.43106397E+00  0.64089798E+00</point>
      <point r="1.100">  0.66139591E+00 -0.73518222E-01  0.41811963E+00  0.63696538E+00</point>
      <point r="1.120">  0.65797056E+00 -0.79921918E-01  0.40502198E+00  0.63294070E+00</point>
      <point r="1.140">  0.65456935E+00 -0.86376147E-01  0.39178397E+00  0.62882619E+00</point>
      <point r="1.160">  0.65118981E+00 -0.92871302E-01  0.37841851E+00  0.62462421E+00</point>
      <point r="1.180">  0.64782942E+00 -0.99397944E-01  0.36493851E+00  0.62033715E+00</point>
      <point r="1.200">  0.64448559E+00 -0.10594682E+00  0.35135679E+00  0.61596746E+00</point>
      <point r="1.220">  0.64115571E+00 -0.11250886E+00  0.33768606E+00  0.61151764E+00</point>
      <point r="1.240">  0.63783714E+00 -0.11907522E+00  0.32393897E+00  0.60699023E+00</point>
      <point r="1.260">  0.63452724E+00 -0.12563725E+00  0.31012799E+00  0.60238782E+00</point>
      <point r="1.280">  0.63122340E+00 -0.13218655E+00  0.29626549E+00  0.59771302E+00</point>
      <point r="1.300">  0.62792301E+00 -0.13871493E+00  0.28236365E+00  0.59296847E+00</point>
      <point r="1.320">  0.62462353E+00 -0.14521446E+00  0.26843448E+00  0.58815684E+00</point>
      <point r="1.340">  0.62132244E+00 -0.15167744E+00  0.25448979E+00  0.58328081E+00</point>
      <point r="1.360">  0.61801730E+00 -0.15809644E+00  0.24054119E+00  0.57834308E+00</point>
      <point r="1.380">  0.61470572E+00 -0.16446428E+00  0.22660007E+00  0.57334636E+00</point>
      <point r="1.400">  0.61138539E+00 -0.17077402E+00  0.21267761E+00  0.56829337E+00</point>
      <point r="1.420">  0.60805408E+00 -0.17701902E+00  0.19878471E+00  0.56318682E+00</point>
      <point r="1.440">  0.60470965E+00 -0.18319287E+00  0.18493205E+00  0.55802943E+00</point>
      <point r="1.460">  0.60135006E+00 -0.18928945E+00  0.17113005E+00  0.55282391E+00</point>
      <point r="1.480">  0.59797335E+00 -0.19530289E+00  0.15738885E+00  0.54757296E+00</point>
      <point r="1.500">  0.59457767E+00 -0.20122761E+00  0.14371833E+00  0.54227928E+00</point>
      <point r="1.520">  0.59116127E+00 -0.20705827E+00  0.13012809E+00  0.53694556E+00</point>
      <point r="1.540">  0.58772251E+00 -0.21278980E+00  0.11662745E+00  0.53157445E+00</point>
      <point r="1.560">  0.58425984E+00 -0.21841743E+00  0.10322542E+00  0.52616859E+00</point>
      <point r="1.580">  0.58077183E+00 -0.22393659E+00  0.89930741E-01  0.52073062E+00</point>
      <point r="1.600">  0.57725718E+00 -0.22934304E+00  0.76751849E-01  0.51526313E+00</point>
      <point r="1.620">  0.57371466E+00 -0.23463273E+00  0.63696880E-01  0.50976869E+00</point>
      <point r="1.640">  0.57014317E+00 -0.23980193E+00  0.50773672E-01  0.50424985E+00</point>
      <point r="1.660">  0.56654173E+00 -0.24484711E+00  0.37989759E-01  0.49870913E+00</point>
      <point r="1.680">  0.56290945E+00 -0.24976501E+00  0.25352371E-01  0.49314900E+00</point>
      <point r="1.700">  0.55924555E+00 -0.25455262E+00  0.12868436E-01  0.48757193E+00</point>
      <point r="1.720">  0.55554936E+00 -0.25920716E+00  0.54457569E-03  0.48198031E+00</point>
      <point r="1.740">  0.55182031E+00 -0.26372609E+00 -0.11612890E-01  0.47637653E+00</point>
      <point r="1.760">  0.54805795E+00 -0.26810710E+00 -0.23597946E-01  0.47076294E+00</point>
      <point r="1.780">  0.54426191E+00 -0.27234810E+00 -0.35404875E-01  0.46514182E+00</point>
      <point r="1.800">  0.54043192E+00 -0.27644725E+00 -0.47028263E-01  0.45951544E+00</point>
      <point r="1.820">  0.53656782E+00 -0.28040289E+00 -0.58462993E-01  0.45388602E+00</point>
      <point r="1.840">  0.53266954E+00 -0.28421358E+00 -0.69704243E-01  0.44825574E+00</point>
      <point r="1.860">  0.52873708E+00 -0.28787812E+00 -0.80747482E-01  0.44262673E+00</point>
      <point r="1.880">  0.52477056E+00 -0.29139547E+00 -0.91588470E-01  0.43700107E+00</point>
      <point r="1.900">  0.52077017E+00 -0.29476480E+00 -0.10222326E+00  0.43138082E+00</point>
      <point r="1.920">  0.51673617E+00 -0.29798548E+00 -0.11264817E+00  0.42576796E+00</point>
      <point r="1.940">  0.51266893E+00 -0.30105704E+00 -0.12285981E+00  0.42016446E+00</point>
      <point r="1.960">  0.50856886E+00 -0.30397923E+00 -0.13285507E+00  0.41457221E+00</point>
      <point r="1.980">  0.50443648E+00 -0.30675193E+00 -0.14263110E+00  0.40899308E+00</point>
      <point r="2.000">  0.50027236E+00 -0.30937521E+00 -0.15218531E+00  0.40342888E+00</point>
      <point r="2.020">  0.49607713E+00 -0.31184930E+00 -0.16151540E+00  0.39788138E+00</point>
      <point r="2.040">  0.49185150E+00 -0.31417459E+00 -0.17061929E+00  0.39235228E+00</point>
      <point r="2.060">  0.48759625E+00 -0.31635160E+00 -0.17949517E+00  0.38684327E+00</point>
      <point r="2.080">  0.48331218E+00 -0.31838102E+00 -0.18814148E+00  0.38135595E+00</point>
      <point r="2.100">  0.47900018E+00 -0.32026366E+00 -0.19655689E+00  0.37589191E+00</point>
      <point r="2.120">  0.47466118E+00 -0.32200048E+00 -0.20474031E+00  0.37045267E+00</point>
      <point r="2.140">  0.47029615E+00 -0.32359256E+00 -0.21269088E+00  0.36503970E+00</point>
      <point r="2.160">  0.46590612E+00 -0.32504109E+00 -0.22040797E+00  0.35965443E+00</point>
      <point r="2.180">  0.46149215E+00 -0.32634739E+00 -0.22789115E+00  0.35429825E+00</point>
      <point r="2.200">  0.45705536E+00 -0.32751290E+00 -0.23514023E+00  0.34897250E+00</point>
      <point r="2.220">  0.45259687E+00 -0.32853914E+00 -0.24215520E+00  0.34367844E+00</point>
      <point r="2.240">  0.44811787E+00 -0.32942776E+00 -0.24893626E+00  0.33841734E+00</point>
      <point r="2.260">  0.44361955E+00 -0.33018048E+00 -0.25548382E+00  0.33319038E+00</point>
      <point r="2.280">  0.43910316E+00 -0.33079913E+00 -0.26179845E+00  0.32799870E+00</point>
      <point r="2.300">  0.43456994E+00 -0.33128560E+00 -0.26788093E+00  0.32284342E+00</point>
      <point r="2.320">  0.43002118E+00 -0.33164189E+00 -0.27373220E+00  0.31772558E+00</point>
      <point r="2.340">  0.42545818E+00 -0.33187005E+00 -0.27935338E+00  0.31264620E+00</point>
      <point r="2.360">  0.42088226E+00 -0.33197222E+00 -0.28474574E+00  0.30760624E+00</point>
      <point r="2.380">  0.41629474E+00 -0.33195058E+00 -0.28991073E+00  0.30260664E+00</point>
      <point r="2.400">  0.41169697E+00 -0.33180740E+00 -0.29484993E+00  0.29764826E+00</point>
      <point r="2.420">  0.40709030E+00 -0.33154498E+00 -0.29956508E+00  0.29273194E+00</point>
      <point r="2.440">  0.40247609E+00 -0.33116569E+00 -0.30405806E+00  0.28785848E+00</point>
      <point r="2.460">  0.39785570E+00 -0.33067194E+00 -0.30833089E+00  0.28302863E+00</point>
      <point r="2.480">  0.39323052E+00 -0.33006618E+00 -0.31238569E+00  0.27824311E+00</point>
      <point r="2.500">  0.38860190E+00 -0.32935091E+00 -0.31622473E+00  0.27350258E+00</point>
      <point r="2.520">  0.38397122E+00 -0.32852865E+00 -0.31985040E+00  0.26880768E+00</point>
      <point r="2.540">  0.37933983E+00 -0.32760195E+00 -0.32326518E+00  0.26415899E+00</point>
      <point r="2.560">  0.37470911E+00 -0.32657341E+00 -0.32647167E+00  0.25955708E+00</point>
      <point r="2.580">  0.37008041E+00 -0.32544562E+00 -0.32947256E+00  0.25500246E+00</point>
      <point r="2.600">  0.36545507E+00 -0.32422123E+00 -0.33227065E+00  0.25049561E+00</point>
      <point r="2.620">  0.36083444E+00 -0.32290287E+00 -0.33486881E+00  0.24603697E+00</point>
      <point r="2.640">  0.35621983E+00 -0.32149319E+00 -0.33727000E+00  0.24162695E+00</point>
      <point r="2.660">  0.35161256E+00 -0.31999487E+00 -0.33947727E+00  0.23726592E+00</point>
      <point r="2.680">  0.34701392E+00 -0.31841058E+00 -0.34149373E+00  0.23295422E+00</point>
      <point r="2.700">  0.34242521E+00 -0.31674298E+00 -0.34332256E+00  0.22869216E+00</point>
      <point r="2.720">  0.33784768E+00 -0.31499476E+00 -0.34496699E+00  0.22448002E+00</point>
      <point r="2.740">  0.33328259E+00 -0.31316860E+00 -0.34643034E+00  0.22031802E+00</point>
      <point r="2.760">  0.32873117E+00 -0.31126715E+00 -0.34771596E+00  0.21620638E+00</point>
      <point r="2.780">  0.32419461E+00 -0.30929308E+00 -0.34882724E+00  0.21214529E+00</point>
      <point r="2.800">  0.31967411E+00 -0.30724905E+00 -0.34976763E+00  0.20813488E+00</point>
      <point r="2.820">  0.31517084E+00 -0.30513767E+00 -0.35054063E+00  0.20417529E+00</point>
      <point r="2.840">  0.31068594E+00 -0.30296159E+00 -0.35114974E+00  0.20026660E+00</point>
      <point r="2.860">  0.30622051E+00 -0.30072340E+00 -0.35159851E+00  0.19640888E+00</point>
      <point r="2.880">  0.30177567E+00 -0.29842570E+00 -0.35189053E+00  0.19260217E+00</point>
      <point r="2.900">  0.29735248E+00 -0.29607103E+00 -0.35202939E+00  0.18884648E+00</point>
      <point r="2.920">  0.29295197E+00 -0.29366195E+00 -0.35201870E+00  0.18514180E+00</point>
      <point r="2.940">  0.28857517E+00 -0.29120097E+00 -0.35186210E+00  0.18148810E+00</point>
      <point r="2.960">  0.28422307E+00 -0.28869059E+00 -0.35156322E+00  0.17788531E+00</point>
      <point r="2.980">  0.27989662E+00 -0.28613325E+00 -0.35112572E+00  0.17433336E+00</point>
      <point r="3.000">  0.27559677E+00 -0.28353140E+00 -0.35055324E+00  0.17083213E+00</point>
      <point r="3.020">  0.27132441E+00 -0.28088744E+00 -0.34984944E+00  0.16738151E+00</point>
      <point r="3.040">  0.26708043E+00 -0.27820372E+00 -0.34901796E+00  0.16398134E+00</point>
      <point r="3.060">  0.26286567E+00 -0.27548259E+00 -0.34806244E+00  0.16063146E+00</point>
      <point r="3.080">  0.25868095E+00 -0.27272634E+00 -0.34698652E+00  0.15733168E+00</point>
      <point r="3.100">  0.25452707E+00 -0.26993724E+00 -0.34579381E+00  0.15408180E+00</point>
      <point r="3.120">  0.25040479E+00 -0.26711751E+00 -0.34448793E+00  0.15088159E+00</point>
      <point r="3.140">  0.24631484E+00 -0.26426933E+00 -0.34307244E+00  0.14773083E+00</point>
      <point r="3.160">  0.24225792E+00 -0.26139486E+00 -0.34155092E+00  0.14462924E+00</point>
      <point r="3.180">  0.23823471E+00 -0.25849621E+00 -0.33992692E+00  0.14157655E+00</point>
      <point r="3.200">  0.23424585E+00 -0.25557545E+00 -0.33820393E+00  0.13857249E+00</point>
      <point r="3.220">  0.23029195E+00 -0.25263459E+00 -0.33638547E+00  0.13561674E+00</point>
      <point r="3.240">  0.22637361E+00 -0.24967563E+00 -0.33447498E+00  0.13270900E+00</point>
      <point r="3.260">  0.22249138E+00 -0.24670051E+00 -0.33247589E+00  0.12984893E+00</point>
      <point r="3.280">  0.21864580E+00 -0.24371114E+00 -0.33039159E+00  0.12703619E+00</point>
      <point r="3.300">  0.21483735E+00 -0.24070936E+00 -0.32822544E+00  0.12427042E+00</point>
      <point r="3.320">  0.21106652E+00 -0.23769700E+00 -0.32598076E+00  0.12155126E+00</point>
      <point r="3.340">  0.20733375E+00 -0.23467583E+00 -0.32366083E+00  0.11887834E+00</point>
      <point r="3.360">  0.20363944E+00 -0.23164758E+00 -0.32126888E+00  0.11625126E+00</point>
      <point r="3.380">  0.19998401E+00 -0.22861392E+00 -0.31880811E+00  0.11366963E+00</point>
      <point r="3.400">  0.19636779E+00 -0.22557651E+00 -0.31628168E+00  0.11113305E+00</point>
      <point r="3.420">  0.19279113E+00 -0.22253694E+00 -0.31369269E+00  0.10864109E+00</point>
      <point r="3.440">  0.18925435E+00 -0.21949676E+00 -0.31104421E+00  0.10619333E+00</point>
      <point r="3.460">  0.18575771E+00 -0.21645748E+00 -0.30833924E+00  0.10378934E+00</point>
      <point r="3.480">  0.18230148E+00 -0.21342057E+00 -0.30558075E+00  0.10142869E+00</point>
      <point r="3.500">  0.17888588E+00 -0.21038745E+00 -0.30277166E+00  0.99110923E-01</point>
      <point r="3.520">  0.17551113E+00 -0.20735950E+00 -0.29991483E+00  0.96835593E-01</point>
      <point r="3.540">  0.17217741E+00 -0.20433805E+00 -0.29701308E+00  0.94602242E-01</point>
      <point r="3.560">  0.16888487E+00 -0.20132440E+00 -0.29406917E+00  0.92410406E-01</point>
      <point r="3.580">  0.16563365E+00 -0.19831979E+00 -0.29108581E+00  0.90259618E-01</point>
      <point r="3.600">  0.16242387E+00 -0.19532544E+00 -0.28806565E+00  0.88149405E-01</point>
      <point r="3.620">  0.15925561E+00 -0.19234251E+00 -0.28501131E+00  0.86079291E-01</point>
      <point r="3.640">  0.15612893E+00 -0.18937212E+00 -0.28192533E+00  0.84048796E-01</point>
      <point r="3.660">  0.15304389E+00 -0.18641536E+00 -0.27881020E+00  0.82057436E-01</point>
      <point r="3.680">  0.15000052E+00 -0.18347328E+00 -0.27566837E+00  0.80104725E-01</point>
      <point r="3.700">  0.14699880E+00 -0.18054686E+00 -0.27250223E+00  0.78190174E-01</point>
      <point r="3.720">  0.14403874E+00 -0.17763708E+00 -0.26931409E+00  0.76313292E-01</point>
      <point r="3.740">  0.14112030E+00 -0.17474486E+00 -0.26610626E+00  0.74473586E-01</point>
      <point r="3.760">  0.13824341E+00 -0.17187109E+00 -0.26288093E+00  0.72670562E-01</point>
      <point r="3.780">  0.13540802E+00 -0.16901660E+00 -0.25964028E+00  0.70903724E-01</point>
      <point r="3.800">  0.13261402E+00 -0.16618222E+00 -0.25638642E+00  0.69172576E-01</point>
      <point r="3.820">  0.12986132E+00 -0.16336871E+00 -0.25312141E+00  0.67476621E-01</point>
      <point r="3.840">  0.12714979E+00 -0.16057681E+00 -0.24984724E+00  0.65815363E-01</point>
      <point r="3.860">  0.12447929E+00 -0.15780721E+00 -0.24656586E+00  0.64188304E-01</point>
      <point r="3.880">  0.12184967E+00 -0.15506059E+00 -0.24327917E+00  0.62594947E-01</point>
      <point r="3.900">  0.11926075E+00 -0.15233757E+00 -0.23998901E+00  0.61034797E-01</point>
      <point r="3.920">  0.11671234E+00 -0.14963876E+00 -0.23669715E+00  0.59507358E-01</point>
      <point r="3.940">  0.11420426E+00 -0.14696470E+00 -0.23340533E+00  0.58012137E-01</point>
      <point r="3.960">  0.11173628E+00 -0.14431595E+00 -0.23011522E+00  0.56548642E-01</point>
      <point r="3.980">  0.10930817E+00 -0.14169298E+00 -0.22682845E+00  0.55116380E-01</point>
      <point r="4.000">  0.10691971E+00 -0.13909628E+00 -0.22354660E+00  0.53714862E-01</point>
      <point r="4.020">  0.10457064E+00 -0.13652628E+00 -0.22027119E+00  0.52343603E-01</point>
      <point r="4.040">  0.10226069E+00 -0.13398339E+00 -0.21700368E+00  0.51002116E-01</point>
      <point r="4.060">  0.99989602E-01 -0.13146798E+00 -0.21374549E+00  0.49689920E-01</point>
      <point r="4.080">  0.97757079E-01 -0.12898041E+00 -0.21049800E+00  0.48406534E-01</point>
      <point r="4.100">  0.95562829E-01 -0.12652100E+00 -0.20726252E+00  0.47151482E-01</point>
      <point r="4.120">  0.93406550E-01 -0.12409004E+00 -0.20404033E+00  0.45924289E-01</point>
      <point r="4.140">  0.91287928E-01 -0.12168780E+00 -0.20083265E+00  0.44724484E-01</point>
      <point r="4.160">  0.89206642E-01 -0.11931453E+00 -0.19764066E+00  0.43551601E-01</point>
      <point r="4.180">  0.87162360E-01 -0.11697044E+00 -0.19446548E+00  0.42405174E-01</point>
      <point r="4.200">  0.85154746E-01 -0.11465573E+00 -0.19130819E+00  0.41284743E-01</point>
      <point r="4.220">  0.83183453E-01 -0.11237057E+00 -0.18816985E+00  0.40189851E-01</point>
      <point r="4.240">  0.81248130E-01 -0.11011510E+00 -0.18505143E+00  0.39120046E-01</point>
      <point r="4.260">  0.79348418E-01 -0.10788946E+00 -0.18195389E+00  0.38074877E-01</point>
      <point r="4.280">  0.77483951E-01 -0.10569374E+00 -0.17887814E+00  0.37053900E-01</point>
      <point r="4.300">  0.75654359E-01 -0.10352802E+00 -0.17582503E+00  0.36056674E-01</point>
      <point r="4.320">  0.73859266E-01 -0.10139237E+00 -0.17279540E+00  0.35082763E-01</point>
      <point r="4.340">  0.72098292E-01 -0.99286836E-01 -0.16979002E+00  0.34131733E-01</point>
      <point r="4.360">  0.70371050E-01 -0.97211434E-01 -0.16680963E+00  0.33203158E-01</point>
      <point r="4.380">  0.68677151E-01 -0.95166172E-01 -0.16385495E+00  0.32296613E-01</point>
      <point r="4.400">  0.67016202E-01 -0.93151038E-01 -0.16092662E+00  0.31411680E-01</point>
      <point r="4.420">  0.65387806E-01 -0.91166001E-01 -0.15802528E+00  0.30547946E-01</point>
      <point r="4.440">  0.63791562E-01 -0.89211015E-01 -0.15515152E+00  0.29705000E-01</point>
      <point r="4.460">  0.62227069E-01 -0.87286019E-01 -0.15230589E+00  0.28882438E-01</point>
      <point r="4.480">  0.60693920E-01 -0.85390935E-01 -0.14948892E+00  0.28079860E-01</point>
      <point r="4.500">  0.59191709E-01 -0.83525671E-01 -0.14670108E+00  0.27296871E-01</point>
      <point r="4.520">  0.57720026E-01 -0.81690121E-01 -0.14394283E+00  0.26533082E-01</point>
      <point r="4.540">  0.56278460E-01 -0.79884165E-01 -0.14121458E+00  0.25788107E-01</point>
      <point r="4.560">  0.54866600E-01 -0.78107670E-01 -0.13851673E+00  0.25061567E-01</point>
      <point r="4.580">  0.53484033E-01 -0.76360492E-01 -0.13584963E+00  0.24353086E-01</point>
      <point r="4.600">  0.52130345E-01 -0.74642473E-01 -0.13321360E+00  0.23662294E-01</point>
      <point r="4.620">  0.50805122E-01 -0.72953445E-01 -0.13060895E+00  0.22988827E-01</point>
      <point r="4.640">  0.49507950E-01 -0.71293228E-01 -0.12803593E+00  0.22332324E-01</point>
      <point r="4.660">  0.48238415E-01 -0.69661633E-01 -0.12549480E+00  0.21692431E-01</point>
      <point r="4.680">  0.46996104E-01 -0.68058460E-01 -0.12298575E+00  0.21068799E-01</point>
      <point r="4.700">  0.45780602E-01 -0.66483500E-01 -0.12050899E+00  0.20461082E-01</point>
      <point r="4.720">  0.44591498E-01 -0.64936536E-01 -0.11806467E+00  0.19868941E-01</point>
      <point r="4.740">  0.43428380E-01 -0.63417341E-01 -0.11565292E+00  0.19292042E-01</point>
      <point r="4.760">  0.42290838E-01 -0.61925680E-01 -0.11327387E+00  0.18730056E-01</point>
      <point r="4.780">  0.41178463E-01 -0.60461313E-01 -0.11092759E+00  0.18182658E-01</point>
      <point r="4.800">  0.40090848E-01 -0.59023990E-01 -0.10861417E+00  0.17649530E-01</point>
      <point r="4.820">  0.39027587E-01 -0.57613456E-01 -0.10633363E+00  0.17130357E-01</point>
      <point r="4.840">  0.37988278E-01 -0.56229449E-01 -0.10408602E+00  0.16624831E-01</point>
      <point r="4.860">  0.36972519E-01 -0.54871701E-01 -0.10187133E+00  0.16132646E-01</point>
      <point r="4.880">  0.35979909E-01 -0.53539939E-01 -0.99689539E-01  0.15653506E-01</point>
      <point r="4.900">  0.35010054E-01 -0.52233886E-01 -0.97540627E-01  0.15187114E-01</point>
      <point r="4.920">  0.34062559E-01 -0.50953258E-01 -0.95424536E-01  0.14733184E-01</point>
      <point r="4.940">  0.33137032E-01 -0.49697769E-01 -0.93341197E-01  0.14291429E-01</point>
      <point r="4.960">  0.32233085E-01 -0.48467126E-01 -0.91290524E-01  0.13861572E-01</point>
      <point r="4.980">  0.31350332E-01 -0.47261035E-01 -0.89272414E-01  0.13443338E-01</point>
      <point r="5.000">  0.30488391E-01 -0.46079199E-01 -0.87286751E-01  0.13036457E-01</point>
      <point r="5.020">  0.29646883E-01 -0.44921315E-01 -0.85333400E-01  0.12640665E-01</point>
      <point r="5.040">  0.28825431E-01 -0.43787080E-01 -0.83412215E-01  0.12255701E-01</point>
      <point r="5.060">  0.28023663E-01 -0.42676187E-01 -0.81523036E-01  0.11881311E-01</point>
      <point r="5.080">  0.27241210E-01 -0.41588329E-01 -0.79665689E-01  0.11517244E-01</point>
      <point r="5.100">  0.26477707E-01 -0.40523195E-01 -0.77839988E-01  0.11163254E-01</point>
      <point r="5.120">  0.25732792E-01 -0.39480474E-01 -0.76045734E-01  0.10819100E-01</point>
      <point r="5.140">  0.25006107E-01 -0.38459852E-01 -0.74282719E-01  0.10484545E-01</point>
      <point r="5.160">  0.24297297E-01 -0.37461016E-01 -0.72550723E-01  0.10159356E-01</point>
      <point r="5.180">  0.23606012E-01 -0.36483650E-01 -0.70849515E-01  0.98433060E-02</point>
      <point r="5.200">  0.22931907E-01 -0.35527440E-01 -0.69178856E-01  0.95361715E-02</point>
      <point r="5.220">  0.22274638E-01 -0.34592069E-01 -0.67538496E-01  0.92377335E-02</point>
      <point r="5.240">  0.21633867E-01 -0.33677223E-01 -0.65928177E-01  0.89477771E-02</point>
      <point r="5.260">  0.21009261E-01 -0.32782585E-01 -0.64347634E-01  0.86660919E-02</point>
      <point r="5.280">  0.20400490E-01 -0.31907840E-01 -0.62796592E-01  0.83924716E-02</point>
      <point r="5.300">  0.19807226E-01 -0.31052673E-01 -0.61274771E-01  0.81267141E-02</point>
      <point r="5.320">  0.19229150E-01 -0.30216771E-01 -0.59781882E-01  0.78686215E-02</point>
      <point r="5.340">  0.18665944E-01 -0.29399819E-01 -0.58317632E-01  0.76179997E-02</point>
      <point r="5.360">  0.18117294E-01 -0.28601506E-01 -0.56881721E-01  0.73746589E-02</point>
      <point r="5.380">  0.17582892E-01 -0.27821521E-01 -0.55473841E-01  0.71384132E-02</point>
      <point r="5.400">  0.17062434E-01 -0.27059554E-01 -0.54093684E-01  0.69090803E-02</point>
      <point r="5.420">  0.16555619E-01 -0.26315298E-01 -0.52740933E-01  0.66864821E-02</point>
      <point r="5.440">  0.16062153E-01 -0.25588445E-01 -0.51415268E-01  0.64704440E-02</point>
      <point r="5.460">  0.15581742E-01 -0.24878691E-01 -0.50116365E-01  0.62607955E-02</point>
      <point r="5.480">  0.15114102E-01 -0.24185733E-01 -0.48843898E-01  0.60573695E-02</point>
      <point r="5.500">  0.14658949E-01 -0.23509272E-01 -0.47597534E-01  0.58600025E-02</point>
      <point r="5.520">  0.14216004E-01 -0.22849007E-01 -0.46376942E-01  0.56685347E-02</point>
      <point r="5.540">  0.13784995E-01 -0.22204644E-01 -0.45181784E-01  0.54828098E-02</point>
      <point r="5.560">  0.13365651E-01 -0.21575887E-01 -0.44011722E-01  0.53026750E-02</point>
      <point r="5.580">  0.12957707E-01 -0.20962446E-01 -0.42866415E-01  0.51279807E-02</point>
      <point r="5.600">  0.12560904E-01 -0.20364032E-01 -0.41745522E-01  0.49585810E-02</point>
      <point r="5.620">  0.12174983E-01 -0.19780359E-01 -0.40648699E-01  0.47943329E-02</point>
      <point r="5.640">  0.11799693E-01 -0.19211143E-01 -0.39575601E-01  0.46350971E-02</point>
      <point r="5.660">  0.11434787E-01 -0.18656104E-01 -0.38525883E-01  0.44807370E-02</point>
      <point r="5.680">  0.11080020E-01 -0.18114963E-01 -0.37499198E-01  0.43311195E-02</point>
      <point r="5.700">  0.10735153E-01 -0.17587447E-01 -0.36495201E-01  0.41861144E-02</point>
      <point r="5.720">  0.10399951E-01 -0.17073283E-01 -0.35513545E-01  0.40455947E-02</point>
      <point r="5.740">  0.10074183E-01 -0.16572202E-01 -0.34553884E-01  0.39094362E-02</point>
      <point r="5.760">  0.97576221E-02 -0.16083939E-01 -0.33615871E-01  0.37775177E-02</point>
      <point r="5.780">  0.94500460E-02 -0.15608232E-01 -0.32699162E-01  0.36497209E-02</point>
      <point r="5.800">  0.91512360E-02 -0.15144821E-01 -0.31803411E-01  0.35259303E-02</point>
      <point r="5.820">  0.88609776E-02 -0.14693449E-01 -0.30928276E-01  0.34060331E-02</point>
      <point r="5.840">  0.85790602E-02 -0.14253865E-01 -0.30073413E-01  0.32899194E-02</point>
      <point r="5.860">  0.83052776E-02 -0.13825820E-01 -0.29238482E-01  0.31774818E-02</point>
      <point r="5.880">  0.80394272E-02 -0.13409066E-01 -0.28423142E-01  0.30686156E-02</point>
      <point r="5.900">  0.77813104E-02 -0.13003361E-01 -0.27627056E-01  0.29632186E-02</point>
      <point r="5.920">  0.75307326E-02 -0.12608465E-01 -0.26849889E-01  0.28611911E-02</point>
      <point r="5.940">  0.72875030E-02 -0.12224144E-01 -0.26091304E-01  0.27624360E-02</point>
      <point r="5.960">  0.70514346E-02 -0.11850164E-01 -0.25350971E-01  0.26668586E-02</point>
      <point r="5.980">  0.68223441E-02 -0.11486296E-01 -0.24628560E-01  0.25743665E-02</point>
      <point r="6.000">  0.66000519E-02 -0.11132315E-01 -0.23923743E-01  0.24848696E-02</point>
      <point r="6.020">  0.63843823E-02 -0.10787998E-01 -0.23236195E-01  0.23982802E-02</point>
      <point r="6.040">  0.61751630E-02 -0.10453127E-01 -0.22565594E-01  0.23145127E-02</point>
      <point r="6.060">  0.59722255E-02 -0.10127485E-01 -0.21911620E-01  0.22334839E-02</point>
      <point r="6.080">  0.57754045E-02 -0.98108627E-02 -0.21273956E-01  0.21551124E-02</point>
      <point r="6.100">  0.55845386E-02 -0.95030500E-02 -0.20652289E-01  0.20793192E-02</point>
      <point r="6.120">  0.53994696E-02 -0.92038421E-02 -0.20046308E-01  0.20060273E-02</point>
      <point r="6.140">  0.52200428E-02 -0.89130377E-02 -0.19455704E-01  0.19351617E-02</point>
      <point r="6.160">  0.50461068E-02 -0.86304384E-02 -0.18880173E-01  0.18666491E-02</point>
      <point r="6.180">  0.48775135E-02 -0.83558496E-02 -0.18319414E-01  0.18004186E-02</point>
      <point r="6.200">  0.47141182E-02 -0.80890798E-02 -0.17773128E-01  0.17364009E-02</point>
      <point r="6.220">  0.45557794E-02 -0.78299411E-02 -0.17241021E-01  0.16745285E-02</point>
      <point r="6.240">  0.44023585E-02 -0.75782488E-02 -0.16722801E-01  0.16147359E-02</point>
      <point r="6.260">  0.42537203E-02 -0.73338215E-02 -0.16218180E-01  0.15569592E-02</point>
      <point r="6.280">  0.41097326E-02 -0.70964813E-02 -0.15726874E-01  0.15011364E-02</point>
      <point r="6.300">  0.39702662E-02 -0.68660534E-02 -0.15248603E-01  0.14472068E-02</point>
      <point r="6.320">  0.38351949E-02 -0.66423664E-02 -0.14783088E-01  0.13951119E-02</point>
      <point r="6.340">  0.37043954E-02 -0.64252520E-02 -0.14330057E-01  0.13447944E-02</point>
      <point r="6.360">  0.35777474E-02 -0.62145453E-02 -0.13889239E-01  0.12961987E-02</point>
      <point r="6.380">  0.34551332E-02 -0.60100843E-02 -0.13460368E-01  0.12492707E-02</point>
      <point r="6.400">  0.33364380E-02 -0.58117104E-02 -0.13043182E-01  0.12039579E-02</point>
      <point r="6.420">  0.32215499E-02 -0.56192680E-02 -0.12637422E-01  0.11602093E-02</point>
      <point r="6.440">  0.31103594E-02 -0.54326047E-02 -0.12242833E-01  0.11179751E-02</point>
      <point r="6.460">  0.30027598E-02 -0.52515710E-02 -0.11859164E-01  0.10772070E-02</point>
      <point r="6.480">  0.28986471E-02 -0.50760206E-02 -0.11486167E-01  0.10378583E-02</point>
      <point r="6.500">  0.27979196E-02 -0.49058101E-02 -0.11123599E-01  0.99988328E-03</point>
      <point r="6.520">  0.27004783E-02 -0.47407990E-02 -0.10771220E-01  0.96323772E-03</point>
      <point r="6.540">  0.26062266E-02 -0.45808499E-02 -0.10428793E-01  0.92787862E-03</point>
      <point r="6.560">  0.25150703E-02 -0.44258281E-02 -0.10096087E-01  0.89376423E-03</point>
      <point r="6.580">  0.24269175E-02 -0.42756020E-02 -0.97728729E-02  0.86085401E-03</point>
      <point r="6.600">  0.23416787E-02 -0.41300425E-02 -0.94589260E-02  0.82910859E-03</point>
      <point r="6.620">  0.22592668E-02 -0.39890235E-02 -0.91540255E-02  0.79848974E-03</point>
      <point r="6.640">  0.21795966E-02 -0.38524215E-02 -0.88579541E-02  0.76896036E-03</point>
      <point r="6.660">  0.21025853E-02 -0.37201160E-02 -0.85704984E-02  0.74048445E-03</point>
      <point r="6.680">  0.20281524E-02 -0.35919887E-02 -0.82914487E-02  0.71302707E-03</point>
      <point r="6.700">  0.19562192E-02 -0.34679245E-02 -0.80205990E-02  0.68655430E-03</point>
      <point r="6.720">  0.18867091E-02 -0.33478103E-02 -0.77577471E-02  0.66103327E-03</point>
      <point r="6.740">  0.18195478E-02 -0.32315359E-02 -0.75026942E-02  0.63643204E-03</point>
      <point r="6.760">  0.17546626E-02 -0.31189937E-02 -0.72552454E-02  0.61271969E-03</point>
      <point r="6.780">  0.16919830E-02 -0.30100783E-02 -0.70152094E-02  0.58986618E-03</point>
      <point r="6.800">  0.16314403E-02 -0.29046869E-02 -0.67823984E-02  0.56784243E-03</point>
      <point r="6.820">  0.15729676E-02 -0.28027191E-02 -0.65566283E-02  0.54662021E-03</point>
      <point r="6.840">  0.15164999E-02 -0.27040768E-02 -0.63377183E-02  0.52617217E-03</point>
      <point r="6.860">  0.14619739E-02 -0.26086642E-02 -0.61254914E-02  0.50647180E-03</point>
      <point r="6.880">  0.14093280E-02 -0.25163878E-02 -0.59197739E-02  0.48749340E-03</point>
      <point r="6.900">  0.13585026E-02 -0.24271565E-02 -0.57203956E-02  0.46921208E-03</point>
      <point r="6.920">  0.13094393E-02 -0.23408813E-02 -0.55271897E-02  0.45160371E-03</point>
      <point r="6.940">  0.12620817E-02 -0.22574752E-02 -0.53399929E-02  0.43464491E-03</point>
      <point r="6.960">  0.12163747E-02 -0.21768535E-02 -0.51586450E-02  0.41831307E-03</point>
      <point r="6.980">  0.11722650E-02 -0.20989336E-02 -0.49829893E-02  0.40258624E-03</point>
      <point r="7.000">  0.11297006E-02 -0.20236350E-02 -0.48128723E-02  0.38744321E-03</point>
      <point r="7.020">  0.10886311E-02 -0.19508790E-02 -0.46481437E-02  0.37286341E-03</point>
      <point r="7.040">  0.10490076E-02 -0.18805891E-02 -0.44886565E-02  0.35882695E-03</point>
      <point r="7.060">  0.10107824E-02 -0.18126906E-02 -0.43342669E-02  0.34531458E-03</point>
      <point r="7.080">  0.97390943E-03 -0.17471107E-02 -0.41848339E-02  0.33230763E-03</point>
      <point r="7.100">  0.93834377E-03 -0.16837786E-02 -0.40402199E-02  0.31978809E-03</point>
      <point r="7.120">  0.90404191E-03 -0.16226253E-02 -0.39002902E-02  0.30773849E-03</point>
      <point r="7.140">  0.87096157E-03 -0.15635834E-02 -0.37649131E-02  0.29614194E-03</point>
      <point r="7.160">  0.83906177E-03 -0.15065876E-02 -0.36339598E-02  0.28498212E-03</point>
      <point r="7.180">  0.80830270E-03 -0.14515739E-02 -0.35073046E-02  0.27424323E-03</point>
      <point r="7.200">  0.77864578E-03 -0.13984803E-02 -0.33848244E-02  0.26391000E-03</point>
      <point r="7.220">  0.75005357E-03 -0.13472465E-02 -0.32663992E-02  0.25396765E-03</point>
      <point r="7.240">  0.72248974E-03 -0.12978135E-02 -0.31519116E-02  0.24440190E-03</point>
      <point r="7.260">  0.69591911E-03 -0.12501241E-02 -0.30412469E-02  0.23519896E-03</point>
      <point r="7.280">  0.67030753E-03 -0.12041228E-02 -0.29342932E-02  0.22634549E-03</point>
      <point r="7.300">  0.64562192E-03 -0.11597552E-02 -0.28309413E-02  0.21782859E-03</point>
      <point r="7.320">  0.62183021E-03 -0.11169688E-02 -0.27310845E-02  0.20963582E-03</point>
      <point r="7.340">  0.59890132E-03 -0.10757124E-02 -0.26346188E-02  0.20175514E-03</point>
      <point r="7.360">  0.57680516E-03 -0.10359360E-02 -0.25414425E-02  0.19417493E-03</point>
      <point r="7.380">  0.55551256E-03 -0.99759147E-03 -0.24514566E-02  0.18688398E-03</point>
      <point r="7.400">  0.53499527E-03 -0.96063160E-03 -0.23645645E-02  0.17987144E-03</point>
      <point r="7.420">  0.51522594E-03 -0.92501071E-03 -0.22806719E-02  0.17312684E-03</point>
      <point r="7.440">  0.49617808E-03 -0.89068440E-03 -0.21996869E-02  0.16664010E-03</point>
      <point r="7.460">  0.47782606E-03 -0.85760953E-03 -0.21215200E-02  0.16040146E-03</point>
      <point r="7.480">  0.46014506E-03 -0.82574417E-03 -0.20460839E-02  0.15440152E-03</point>
      <point r="7.500">  0.44311107E-03 -0.79504766E-03 -0.19732935E-02  0.14863119E-03</point>
      <point r="7.520">  0.42670084E-03 -0.76548048E-03 -0.19030659E-02  0.14308172E-03</point>
      <point r="7.540">  0.41089190E-03 -0.73700427E-03 -0.18353204E-02  0.13774466E-03</point>
      <point r="7.560">  0.39566251E-03 -0.70958181E-03 -0.17699784E-02  0.13261186E-03</point>
      <point r="7.580">  0.38099163E-03 -0.68317699E-03 -0.17069634E-02  0.12767547E-03</point>
      <point r="7.600">  0.36685895E-03 -0.65775474E-03 -0.16462007E-02  0.12292791E-03</point>
      <point r="7.620">  0.35324479E-03 -0.63328106E-03 -0.15876178E-02  0.11836188E-03</point>
      <point r="7.640">  0.34013017E-03 -0.60972299E-03 -0.15311442E-02  0.11397034E-03</point>
      <point r="7.660">  0.32749672E-03 -0.58704853E-03 -0.14767112E-02  0.10974650E-03</point>
      <point r="7.680">  0.31532670E-03 -0.56522669E-03 -0.14242519E-02  0.10568382E-03</point>
      <point r="7.700">  0.30360297E-03 -0.54422739E-03 -0.13737012E-02  0.10177602E-03</point>
      <point r="7.720">  0.29230897E-03 -0.52402150E-03 -0.13249961E-02  0.98017023E-04</point>
      <point r="7.740">  0.28142872E-03 -0.50458080E-03 -0.12780751E-02  0.94400985E-04</point>
      <point r="7.760">  0.27094678E-03 -0.48587791E-03 -0.12328784E-02  0.90922283E-04</point>
      <point r="7.780">  0.26084825E-03 -0.46788635E-03 -0.11893479E-02  0.87575499E-04</point>
      <point r="7.800">  0.25111873E-03 -0.45058043E-03 -0.11474273E-02  0.84355419E-04</point>
      <point r="7.820">  0.24174435E-03 -0.43393532E-03 -0.11070618E-02  0.81257023E-04</point>
      <point r="7.840">  0.23271172E-03 -0.41792693E-03 -0.10681980E-02  0.78275483E-04</point>
      <point r="7.860">  0.22400791E-03 -0.40253197E-03 -0.10307843E-02  0.75406150E-04</point>
      <point r="7.880">  0.21562045E-03 -0.38772790E-03 -0.99477058E-03  0.72644553E-04</point>
      <point r="7.900">  0.20753734E-03 -0.37349289E-03 -0.96010799E-03  0.69986393E-04</point>
      <point r="7.920">  0.19974697E-03 -0.35980583E-03 -0.92674927E-03  0.67427534E-04</point>
      <point r="7.940">  0.19223818E-03 -0.34664631E-03 -0.89464853E-03  0.64963998E-04</point>
      <point r="7.960">  0.18500020E-03 -0.33399458E-03 -0.86376124E-03  0.62591964E-04</point>
      <point r="7.980">  0.17802265E-03 -0.32183153E-03 -0.83404420E-03  0.60307757E-04</point>
      <point r="8.000">  0.17129555E-03 -0.31013871E-03 -0.80545553E-03  0.58107843E-04</point>
      <point r="8.020">  0.16480925E-03 -0.29889828E-03 -0.77795462E-03  0.55988831E-04</point>
      <point r="8.040">  0.15855450E-03 -0.28809298E-03 -0.75150209E-03  0.53947459E-04</point>
      <point r="8.060">  0.15252234E-03 -0.27770617E-03 -0.72605978E-03  0.51980596E-04</point>
      <point r="8.080">  0.14670420E-03 -0.26772173E-03 -0.70159071E-03  0.50085234E-04</point>
      <point r="8.100">  0.14109179E-03 -0.25812412E-03 -0.67805907E-03  0.48258485E-04</point>
      <point r="8.120">  0.13567716E-03 -0.24889834E-03 -0.65543015E-03  0.46497575E-04</point>
      <point r="8.140">  0.13045263E-03 -0.24002990E-03 -0.63367035E-03  0.44799844E-04</point>
      <point r="8.160">  0.12541085E-03 -0.23150480E-03 -0.61274714E-03  0.43162736E-04</point>
      <point r="8.180">  0.12054471E-03 -0.22330954E-03 -0.59262902E-03  0.41583799E-04</point>
      <point r="8.200">  0.11584741E-03 -0.21543110E-03 -0.57328552E-03  0.40060683E-04</point>
      <point r="8.220">  0.11131239E-03 -0.20785692E-03 -0.55468716E-03  0.38591130E-04</point>
      <point r="8.240">  0.10693335E-03 -0.20057488E-03 -0.53680539E-03  0.37172977E-04</point>
      <point r="8.260">  0.10270424E-03 -0.19357329E-03 -0.51961265E-03  0.35804150E-04</point>
      <point r="8.280">  0.98619239E-04 -0.18684090E-03 -0.50308224E-03  0.34482659E-04</point>
      <point r="8.300">  0.94672759E-04 -0.18036685E-03 -0.48718839E-03  0.33206596E-04</point>
      <point r="8.320">  0.90859435E-04 -0.17414068E-03 -0.47190617E-03  0.31974135E-04</point>
      <point r="8.340">  0.87174111E-04 -0.16815231E-03 -0.45721149E-03  0.30783523E-04</point>
      <point r="8.360">  0.83611840E-04 -0.16239203E-03 -0.44308110E-03  0.29633084E-04</point>
      <point r="8.380">  0.80167872E-04 -0.15685051E-03 -0.42949253E-03  0.28521208E-04</point>
      <point r="8.400">  0.76837648E-04 -0.15151874E-03 -0.41642407E-03  0.27446358E-04</point>
      <point r="8.420">  0.73616796E-04 -0.14638806E-03 -0.40385480E-03  0.26407059E-04</point>
      <point r="8.440">  0.70501118E-04 -0.14145014E-03 -0.39176449E-03  0.25401899E-04</point>
      <point r="8.460">  0.67486591E-04 -0.13669697E-03 -0.38013365E-03  0.24429526E-04</point>
      <point r="8.480">  0.64569356E-04 -0.13212084E-03 -0.36894347E-03  0.23488649E-04</point>
      <point r="8.500">  0.61745711E-04 -0.12771432E-03 -0.35817580E-03  0.22578028E-04</point>
      <point r="8.520">  0.59012112E-04 -0.12347030E-03 -0.34781315E-03  0.21696480E-04</point>
      <point r="8.540">  0.56365157E-04 -0.11938192E-03 -0.33783866E-03  0.20842872E-04</point>
      <point r="8.560">  0.53801591E-04 -0.11544261E-03 -0.32823610E-03  0.20016120E-04</point>
      <point r="8.580">  0.51318292E-04 -0.11164603E-03 -0.31898980E-03  0.19215187E-04</point>
      <point r="8.600">  0.48912273E-04 -0.10798613E-03 -0.31008471E-03  0.18439083E-04</point>
      <point r="8.620">  0.46580670E-04 -0.10445706E-03 -0.30150629E-03  0.17686861E-04</point>
      <point r="8.640">  0.44320743E-04 -0.10105325E-03 -0.29324060E-03  0.16957613E-04</point>
      <point r="8.660">  0.42129869E-04 -0.97769307E-04 -0.28527419E-03  0.16250476E-04</point>
      <point r="8.680">  0.40005536E-04 -0.94600094E-04 -0.27759412E-03  0.15564620E-04</point>
      <point r="8.700">  0.37945340E-04 -0.91540671E-04 -0.27018796E-03  0.14899256E-04</point>
      <point r="8.720">  0.35946983E-04 -0.88586301E-04 -0.26304375E-03  0.14253628E-04</point>
      <point r="8.740">  0.34008265E-04 -0.85732447E-04 -0.25615001E-03  0.13627015E-04</point>
      <point r="8.760">  0.32127081E-04 -0.82974759E-04 -0.24949569E-03  0.13018726E-04</point>
      <point r="8.780">  0.30301419E-04 -0.80309072E-04 -0.24307019E-03  0.12428104E-04</point>
      <point r="8.800">  0.28529355E-04 -0.77731397E-04 -0.23686332E-03  0.11854520E-04</point>
      <point r="8.820">  0.26809050E-04 -0.75237914E-04 -0.23086531E-03  0.11297373E-04</point>
      <point r="8.840">  0.25138744E-04 -0.72824969E-04 -0.22506677E-03  0.10756090E-04</point>
      <point r="8.860">  0.23516757E-04 -0.70489065E-04 -0.21945872E-03  0.10230123E-04</point>
      <point r="8.880">  0.21941482E-04 -0.68226857E-04 -0.21403252E-03  0.97189499E-05</point>
      <point r="8.900">  0.20411384E-04 -0.66035149E-04 -0.20877989E-03  0.92220712E-05</point>
      <point r="8.920">  0.18924995E-04 -0.63910885E-04 -0.20369292E-03  0.87390106E-05</point>
      <point r="8.940">  0.17480915E-04 -0.61851143E-04 -0.19876402E-03  0.82693135E-05</point>
      <point r="8.960">  0.16077804E-04 -0.59853136E-04 -0.19398590E-03  0.78125457E-05</point>
      <point r="8.980">  0.14714383E-04 -0.57914199E-04 -0.18935162E-03  0.73682926E-05</point>
      <point r="9.000">  0.13389431E-04 -0.56031792E-04 -0.18485452E-03  0.69361589E-05</point>
      <point r="9.020">  0.12101780E-04 -0.54203489E-04 -0.18048823E-03  0.65157666E-05</point>
      <point r="9.040">  0.10850315E-04 -0.52426977E-04 -0.17624667E-03  0.61067554E-05</point>
      <point r="9.060">  0.96339720E-05 -0.50700051E-04 -0.17212401E-03  0.57087810E-05</point>
      <point r="9.080">  0.84517332E-05 -0.49020609E-04 -0.16811471E-03  0.53215146E-05</point>
      <point r="9.100">  0.73026273E-05 -0.47386648E-04 -0.16421345E-03  0.49446425E-05</point>
      <point r="9.120">  0.61857262E-05 -0.45796262E-04 -0.16041518E-03  0.45778650E-05</point>
      <point r="9.140">  0.51001434E-05 -0.44247635E-04 -0.15671507E-03  0.42208956E-05</point>
      <point r="9.160">  0.40450317E-05 -0.42739041E-04 -0.15310852E-03  0.38734609E-05</point>
      <point r="9.180">  0.30195818E-05 -0.41268835E-04 -0.14959113E-03  0.35352996E-05</point>
      <point r="9.200">  0.20230202E-05 -0.39835456E-04 -0.14615873E-03  0.32061620E-05</point>
      <point r="9.220">  0.10546076E-05 -0.38437419E-04 -0.14280735E-03  0.28858092E-05</point>
      <point r="9.240">  0.11363728E-06 -0.37073314E-04 -0.13953319E-03  0.25740132E-05</point>
      <point r="9.260"> -0.80056659E-06 -0.35741801E-04 -0.13633267E-03  0.22705554E-05</point>
      <point r="9.280"> -0.16886504E-05 -0.34441610E-04 -0.13320235E-03  0.19752270E-05</point>
      <point r="9.300"> -0.25512325E-05 -0.33171535E-04 -0.13013900E-03  0.16878282E-05</point>
      <point r="9.320"> -0.33889047E-05 -0.31930432E-04 -0.12713953E-03  0.14081675E-05</point>
      <point r="9.340"> -0.42022336E-05 -0.30717217E-04 -0.12420101E-03  0.11360616E-05</point>
      <point r="9.360"> -0.49917618E-05 -0.29530863E-04 -0.12132067E-03  0.87133471E-06</point>
      <point r="9.380"> -0.57580092E-05 -0.28370399E-04 -0.11849589E-03  0.61381853E-06</point>
      <point r="9.400"> -0.65014743E-05 -0.27234903E-04 -0.11572417E-03  0.36335145E-06</point>
      <point r="9.420"> -0.72226349E-05 -0.26123506E-04 -0.11300317E-03  0.11977843E-06</point>
      <point r="9.440"> -0.79219497E-05 -0.25035384E-04 -0.11033066E-03 -0.11704942E-06</point>
      <point r="9.460"> -0.85998590E-05 -0.23969759E-04 -0.10770454E-03 -0.34727516E-06</point>
      <point r="9.480"> -0.92567860E-05 -0.22925896E-04 -0.10512283E-03 -0.57103632E-06</point>
      <point r="9.500"> -0.98931373E-05 -0.21903102E-04 -0.10258366E-03 -0.78846524E-06</point>
      <point r="9.520"> -0.10509304E-04 -0.20900721E-04 -0.10008526E-03 -0.99968933E-06</point>
      <point r="9.540"> -0.11105663E-04 -0.19918136E-04 -0.97625974E-04 -0.12048313E-05</point>
      <point r="9.560"> -0.11682577E-04 -0.18954765E-04 -0.95204248E-04 -0.14040096E-05</point>
      <point r="9.580"> -0.12240397E-04 -0.18010059E-04 -0.92818614E-04 -0.15973383E-05</point>
      <point r="9.600"> -0.12779459E-04 -0.17083502E-04 -0.90467697E-04 -0.17849279E-05</point>
      <point r="9.620"> -0.13300090E-04 -0.16174608E-04 -0.88150208E-04 -0.19668848E-05</point>
      <point r="9.640"> -0.13802606E-04 -0.15282918E-04 -0.85864938E-04 -0.21433123E-05</point>
      <point r="9.660"> -0.14287312E-04 -0.14408004E-04 -0.83610760E-04 -0.23143102E-05</point>
      <point r="9.680"> -0.14754503E-04 -0.13549460E-04 -0.81386618E-04 -0.24799754E-05</point>
      <point r="9.700"> -0.15204466E-04 -0.12706907E-04 -0.79191528E-04 -0.26404016E-05</point>
      <point r="9.720"> -0.15637479E-04 -0.11879987E-04 -0.77024574E-04 -0.27956800E-05</point>
      <point r="9.740"> -0.16053813E-04 -0.11068367E-04 -0.74884903E-04 -0.29458994E-05</point>
      <point r="9.760"> -0.16453732E-04 -0.10271732E-04 -0.72771724E-04 -0.30911458E-05</point>
      <point r="9.780"> -0.16837490E-04 -0.94897869E-05 -0.70684302E-04 -0.32315031E-05</point>
      <point r="9.800"> -0.17205337E-04 -0.87222556E-05 -0.68621960E-04 -0.33670533E-05</point>
      <point r="9.820"> -0.17557519E-04 -0.79688792E-05 -0.66584069E-04 -0.34978760E-05</point>
      <point r="9.840"> -0.17894271E-04 -0.72294147E-05 -0.64570052E-04 -0.36240492E-05</point>
      <point r="9.860"> -0.18215828E-04 -0.65036348E-05 -0.62579377E-04 -0.37456491E-05</point>
      <point r="9.880"> -0.18522417E-04 -0.57913262E-05 -0.60611557E-04 -0.38627501E-05</point>
      <point r="9.900"> -0.18814261E-04 -0.50922891E-05 -0.58666145E-04 -0.39754254E-05</point>
      <point r="9.920"> -0.19091580E-04 -0.44063366E-05 -0.56742734E-04 -0.40837464E-05</point>
      <point r="9.940"> -0.19354590E-04 -0.37332931E-05 -0.54840955E-04 -0.41877832E-05</point>
      <point r="9.960"> -0.19603501E-04 -0.30729942E-05 -0.52960471E-04 -0.42876048E-05</point>
      <point r="9.980"> -0.19838522E-04 -0.24252857E-05 -0.51100981E-04 -0.43832789E-05</point>
      <point r="10.000"> -0.20059858E-04 -0.17900230E-05 -0.49262210E-04 -0.44748721E-05</point>
      <point r="10.020"> -0.20267711E-04 -0.11670701E-05 -0.47443916E-04 -0.45624500E-05</point>
      <point r="10.040"> -0.20462282E-04 -0.55629947E-06 -0.45645882E-04 -0.46460770E-05</point>
      <point r="10.060"> -0.20643768E-04  0.42409162E-07 -0.43867916E-04 -0.47258169E-05</point>
      <point r="10.080"> -0.20812363E-04  0.62916891E-06 -0.42109849E-04 -0.48017324E-05</point>
      <point r="10.100"> -0.20968260E-04  0.12040865E-05 -0.40371534E-04 -0.48738855E-05</point>
      <point r="10.120"> -0.21111650E-04  0.17672629E-05 -0.38652844E-04 -0.49423374E-05</point>
      <point r="10.140"> -0.21242724E-04  0.23187935E-05 -0.36953672E-04 -0.50071486E-05</point>
      <point r="10.160"> -0.21361667E-04  0.28587690E-05 -0.35273926E-04 -0.50683790E-05</point>
      <point r="10.180"> -0.21468667E-04  0.33872754E-05 -0.33613531E-04 -0.51260878E-05</point>
      <point r="10.200"> -0.21563908E-04  0.39043950E-05 -0.31972428E-04 -0.51803334E-05</point>
      <point r="10.220"> -0.21647573E-04  0.44102064E-05 -0.30350568E-04 -0.52311740E-05</point>
      <point r="10.240"> -0.21719846E-04  0.49047847E-05 -0.28747917E-04 -0.52786671E-05</point>
      <point r="10.260"> -0.21780908E-04  0.53882026E-05 -0.27164452E-04 -0.53228697E-05</point>
      <point r="10.280"> -0.21830938E-04  0.58605300E-05 -0.25600158E-04 -0.53638383E-05</point>
      <point r="10.300"> -0.21870116E-04  0.63218347E-05 -0.24055031E-04 -0.54016289E-05</point>
      <point r="10.320"> -0.21898622E-04  0.67721825E-05 -0.22529073E-04 -0.54362971E-05</point>
      <point r="10.340"> -0.21916632E-04  0.72116377E-05 -0.21022295E-04 -0.54678983E-05</point>
      <point r="10.360"> -0.21924325E-04  0.76402634E-05 -0.19534714E-04 -0.54964871E-05</point>
      <point r="10.380"> -0.21921876E-04  0.80581214E-05 -0.18066351E-04 -0.55221180E-05</point>
      <point r="10.400">  0.00000000E+00</point>
    </S_spline>
    <Vrep_spline>
      <point r="0.000">  0.37750294E+03</point>
      <point r="0.100">  0.28193386E+03</point>
      <point r="0.200">  0.21033306E+03</point>
      <point r="0.300">  0.15668940E+03</point>
      <point r="0.400">  0.11649933E+03</point>
      <point r="0.500">  0.86388733E+02</point>
      <point r="0.600">  0.63829735E+02</point>
      <point r="0.700">  0.46928429E+02</point>
      <point r="0.800">  0.34265891E+02</point>
      <point r="0.900">  0.24779059E+02</point>
      <point r="1.000">  0.17671480E+02</point>
      <point r="1.058">  0.14395043E+02</point>
      <point r="1.117">  0.11723444E+02</point>
      <point r="1.175">  0.95534825E+01</point>
      <point r="1.234">  0.77972440E+01</point>
      <point r="1.292">  0.63802484E+01</point>
      <point r="1.351">  0.52397649E+01</point>
      <point r="1.409">  0.43232896E+01</point>
      <point r="1.468">  0.35871732E+01</point>
      <point r="1.526">  0.29953890E+01</point>
      <point r="1.585">  0.25184317E+01</point>
      <point r="1.643">  0.21323373E+01</point>
      <point r="1.701">  0.18178159E+01</point>
      <point r="1.760">  0.15594878E+01</point>
      <point r="1.818">  0.13452147E+01</point>
      <point r="1.877">  0.11655188E+01</point>
      <point r="1.935">  0.10130810E+01</point>
      <point r="1.994">  0.88231199E+00</point>
      <point r="2.052">  0.76898840E+00</point>
      <point r="2.111">  0.66994751E+00</point>
      <point r="2.169">  0.58283470E+00</point>
      <point r="2.228">  0.50589712E+00</point>
      <point r="2.286">  0.43781821E+00</point>
      <point r="2.344">  0.37758773E+00</point>
      <point r="2.403">  0.32440228E+00</point>
      <point r="2.461">  0.27759177E+00</point>
      <point r="2.520">  0.23656748E+00</point>
      <point r="2.578">  0.20078787E+00</point>
      <point r="2.637">  0.16973837E+00</point>
      <point r="2.695">  0.14292205E+00</point>
      <point r="2.754">  0.11985815E+00</point>
      <point r="2.812">  0.10008580E+00</point>
      <point r="2.870">  0.83170890E-01</point>
      <point r="2.929">  0.68713934E-01</point>
      <point r="2.987">  0.56357484E-01</point>
      <point r="3.046">  0.45791817E-01</point>
      <point r="3.104">  0.36757980E-01</point>
      <point r="3.163">  0.29047634E-01</point>
      <point r="3.221">  0.22499458E-01</point>
      <point r="3.280">  0.16992229E-01</point>
      <point r="3.338">  0.12435013E-01</point>
      <point r="3.397">  0.87552610E-02</point>
      <point r="3.455">  0.58859192E-02</point>
      <point r="3.513">  0.37530168E-02</point>
      <point r="3.572">  0.22655281E-02</point>
      <point r="3.630">  0.13096400E-02</point>
      <point r="3.689">  0.74989755E-03</point>
      <point r="3.747">  0.44003621E-03</point>
      <point r="3.806">  0.24664847E-03</point>
      <point r="3.864">  0.89168817E-04</point>
      <point r="3.923">  0.00000000E+00</point>
    </Vrep_spline>
  </per_pair_data>
  <per_pair_data type1="2" type2="2" SK_cutoff="10.4000000000000004"
    Vrep_cutoff="4.7999999999999998" SK_npts="521" Vrep_npts="61">
    <H_spline>
      <point r="0.000">  0.48000000E+01  0.00000000E+00 -0.14661800E-01  0.26325200E-02  0.00000000E+00 -0.12866900E-01  0.25452700E-01  0.28086000E+02  0.00000000E+00  0.15585500E-01</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400"> -0.93528652E+00  0.25314733E+00 -0.32181570E-01 -0.84989682E+00</point>
      <point r="0.420"> -0.92514432E+00  0.28097523E+00  0.89243600E-02 -0.80991597E+00</point>
      <point r="0.440"> -0.91492284E+00  0.30456373E+00  0.43825333E-01 -0.77185049E+00</point>
      <point r="0.460"> -0.90427665E+00  0.32403977E+00  0.72919315E-01 -0.73570057E+00</point>
      <point r="0.480"> -0.89297483E+00  0.33959181E+00  0.96634826E-01 -0.70146376E+00</point>
      <point r="0.500"> -0.88086460E+00  0.35144784E+00  0.11541699E+00 -0.66911928E+00</point>
      <point r="0.520"> -0.86786226E+00  0.35986341E+00  0.12969617E+00 -0.63862563E+00</point>
      <point r="0.540"> -0.85394855E+00  0.36511436E+00  0.13989157E+00 -0.60993961E+00</point>
      <point r="0.560"> -0.83913914E+00  0.36748360E+00  0.14641716E+00 -0.58300560E+00</point>
      <point r="0.580"> -0.82348157E+00  0.36725573E+00  0.14967037E+00 -0.55776033E+00</point>
      <point r="0.600"> -0.80705328E+00  0.36471299E+00  0.15002240E+00 -0.53413934E+00</point>
      <point r="0.620"> -0.78994713E+00  0.36012936E+00  0.14782146E+00 -0.51207354E+00</point>
      <point r="0.640"> -0.77226425E+00  0.35376748E+00  0.14339696E+00 -0.49149092E+00</point>
      <point r="0.660"> -0.75411219E+00  0.34587688E+00  0.13705650E+00 -0.47231950E+00</point>
      <point r="0.680"> -0.73560155E+00  0.33669220E+00  0.12908209E+00 -0.45448722E+00</point>
      <point r="0.700"> -0.71684235E+00  0.32643233E+00  0.11973134E+00 -0.43792226E+00</point>
      <point r="0.720"> -0.69794102E+00  0.31529998E+00  0.10924080E+00 -0.42255384E+00</point>
      <point r="0.740"> -0.67899855E+00  0.30348162E+00  0.97827610E-01 -0.40831272E+00</point>
      <point r="0.760"> -0.66010969E+00  0.29114765E+00  0.85689381E-01 -0.39513152E+00</point>
      <point r="0.780"> -0.64136250E+00  0.27845286E+00  0.73004284E-01 -0.38294501E+00</point>
      <point r="0.800"> -0.62283767E+00  0.26553692E+00  0.59932234E-01 -0.37169027E+00</point>
      <point r="0.820"> -0.60460639E+00  0.25252479E+00  0.46620739E-01 -0.36130670E+00</point>
      <point r="0.840"> -0.58673874E+00  0.23952877E+00  0.33185280E-01 -0.35173660E+00</point>
      <point r="0.860"> -0.56928847E+00  0.22664730E+00  0.19748118E-01 -0.34292438E+00</point>
      <point r="0.880"> -0.55230382E+00  0.21396688E+00  0.64130810E-02 -0.33481678E+00</point>
      <point r="0.900"> -0.53582718E+00  0.20156326E+00 -0.67299623E-02 -0.32736367E+00</point>
      <point r="0.920"> -0.51989527E+00  0.18950216E+00 -0.19604635E-01 -0.32051774E+00</point>
      <point r="0.940"> -0.50453853E+00  0.17783971E+00 -0.32146089E-01 -0.31423402E+00</point>
      <point r="0.960"> -0.48978075E+00  0.16662318E+00 -0.44298397E-01 -0.30846967E+00</point>
      <point r="0.980"> -0.47563943E+00  0.15589171E+00 -0.56012493E-01 -0.30318413E+00</point>
      <point r="1.000"> -0.46212642E+00  0.14567711E+00 -0.67245301E-01 -0.29833926E+00</point>
      <point r="1.020"> -0.44924872E+00  0.13600452E+00 -0.77959946E-01 -0.29389927E+00</point>
      <point r="1.040"> -0.43700935E+00  0.12689311E+00 -0.88126341E-01 -0.28983059E+00</point>
      <point r="1.060"> -0.42540796E+00  0.11835670E+00 -0.97721479E-01 -0.28610173E+00</point>
      <point r="1.080"> -0.41444127E+00  0.11040432E+00 -0.10672906E+00 -0.28268320E+00</point>
      <point r="1.100"> -0.40410319E+00  0.10304065E+00 -0.11513850E+00 -0.27954735E+00</point>
      <point r="1.120"> -0.39438480E+00  0.96266445E-01 -0.12294365E+00 -0.27666829E+00</point>
      <point r="1.140"> -0.38527447E+00  0.90078910E-01 -0.13014160E+00 -0.27402179E+00</point>
      <point r="1.160"> -0.37675809E+00  0.84472087E-01 -0.13673194E+00 -0.27158521E+00</point>
      <point r="1.180"> -0.36881946E+00  0.79437242E-01 -0.14271650E+00 -0.26933750E+00</point>
      <point r="1.200"> -0.36144080E+00  0.74963250E-01 -0.14809940E+00 -0.26725908E+00</point>
      <point r="1.220"> -0.35460318E+00  0.71036962E-01 -0.15288740E+00 -0.26533182E+00</point>
      <point r="1.240"> -0.34828688E+00  0.67643520E-01 -0.15708998E+00 -0.26353896E+00</point>
      <point r="1.260"> -0.34247165E+00  0.64766638E-01 -0.16071925E+00 -0.26186496E+00</point>
      <point r="1.280"> -0.33713680E+00  0.62388825E-01 -0.16378962E+00 -0.26029549E+00</point>
      <point r="1.300"> -0.33226122E+00  0.60491583E-01 -0.16631721E+00 -0.25881728E+00</point>
      <point r="1.320"> -0.32782344E+00  0.59055580E-01 -0.16831925E+00 -0.25741808E+00</point>
      <point r="1.340"> -0.32380163E+00  0.58060808E-01 -0.16981363E+00 -0.25608659E+00</point>
      <point r="1.360"> -0.32017378E+00  0.57486747E-01 -0.17081852E+00 -0.25481245E+00</point>
      <point r="1.380"> -0.31691778E+00  0.57312516E-01 -0.17135232E+00 -0.25358615E+00</point>
      <point r="1.400"> -0.31401166E+00  0.57517031E-01 -0.17143367E+00 -0.25239902E+00</point>
      <point r="1.420"> -0.31143373E+00  0.58079141E-01 -0.17108155E+00 -0.25124315E+00</point>
      <point r="1.440"> -0.30916276E+00  0.58977753E-01 -0.17031544E+00 -0.25011139E+00</point>
      <point r="1.460"> -0.30717803E+00  0.60191933E-01 -0.16915525E+00 -0.24899724E+00</point>
      <point r="1.480"> -0.30545943E+00  0.61700994E-01 -0.16762131E+00 -0.24789485E+00</point>
      <point r="1.500"> -0.30398741E+00  0.63484560E-01 -0.16573413E+00 -0.24679893E+00</point>
      <point r="1.520"> -0.30274303E+00  0.65522619E-01 -0.16351426E+00 -0.24570474E+00</point>
      <point r="1.540"> -0.30170792E+00  0.67795574E-01 -0.16098199E+00 -0.24460801E+00</point>
      <point r="1.560"> -0.30086428E+00  0.70284289E-01 -0.15815721E+00 -0.24350497E+00</point>
      <point r="1.580"> -0.30019490E+00  0.72970127E-01 -0.15505932E+00 -0.24239225E+00</point>
      <point r="1.600"> -0.29968320E+00  0.75834997E-01 -0.15170714E+00 -0.24126690E+00</point>
      <point r="1.620"> -0.29931326E+00  0.78861391E-01 -0.14811893E+00 -0.24012634E+00</point>
      <point r="1.640"> -0.29906986E+00  0.82032416E-01 -0.14431246E+00 -0.23896836E+00</point>
      <point r="1.660"> -0.29893323E+00  0.85329176E-01 -0.14029147E+00 -0.23779120E+00</point>
      <point r="1.680"> -0.29890195E+00  0.88742165E-01 -0.13610327E+00 -0.23659325E+00</point>
      <point r="1.700"> -0.29895643E+00  0.92253264E-01 -0.13174909E+00 -0.23537297E+00</point>
      <point r="1.720"> -0.29908395E+00  0.95848076E-01 -0.12724374E+00 -0.23412925E+00</point>
      <point r="1.740"> -0.29927248E+00  0.99512861E-01 -0.12260139E+00 -0.23286123E+00</point>
      <point r="1.760"> -0.29951071E+00  0.10323455E+00 -0.11783576E+00 -0.23156826E+00</point>
      <point r="1.780"> -0.29978811E+00  0.10700077E+00 -0.11296021E+00 -0.23024993E+00</point>
      <point r="1.800"> -0.30009490E+00  0.11079980E+00 -0.10798787E+00 -0.22890598E+00</point>
      <point r="1.820"> -0.30042204E+00  0.11462061E+00 -0.10293164E+00 -0.22753633E+00</point>
      <point r="1.840"> -0.30076117E+00  0.11845282E+00 -0.97804147E-01 -0.22614102E+00</point>
      <point r="1.860"> -0.30110461E+00  0.12228669E+00 -0.92617723E-01 -0.22472023E+00</point>
      <point r="1.880"> -0.30144525E+00  0.12611307E+00 -0.87384257E-01 -0.22327425E+00</point>
      <point r="1.900"> -0.30177652E+00  0.12992341E+00 -0.82115117E-01 -0.22180343E+00</point>
      <point r="1.920"> -0.30209233E+00  0.13370971E+00 -0.76821049E-01 -0.22030824E+00</point>
      <point r="1.940"> -0.30238705E+00  0.13746447E+00 -0.71512100E-01 -0.21878919E+00</point>
      <point r="1.960"> -0.30265544E+00  0.14118072E+00 -0.66197572E-01 -0.21724688E+00</point>
      <point r="1.980"> -0.30289264E+00  0.14485197E+00 -0.60886014E-01 -0.21568196E+00</point>
      <point r="2.000"> -0.30309419E+00  0.14847221E+00 -0.55585235E-01 -0.21409513E+00</point>
      <point r="2.020"> -0.30325596E+00  0.15203585E+00 -0.50302354E-01 -0.21248715E+00</point>
      <point r="2.040"> -0.30337423E+00  0.15553778E+00 -0.45043868E-01 -0.21085881E+00</point>
      <point r="2.060"> -0.30344560E+00  0.15897328E+00 -0.39815729E-01 -0.20921096E+00</point>
      <point r="2.080"> -0.30346707E+00  0.16233808E+00 -0.34623428E-01 -0.20754448E+00</point>
      <point r="2.100"> -0.30343599E+00  0.16562830E+00 -0.29472070E-01 -0.20586027E+00</point>
      <point r="2.120"> -0.30335005E+00  0.16884042E+00 -0.24366439E-01 -0.20415925E+00</point>
      <point r="2.140"> -0.30320730E+00  0.17197134E+00 -0.19311054E-01 -0.20244237E+00</point>
      <point r="2.160"> -0.30300611E+00  0.17501829E+00 -0.14310195E-01 -0.20071059E+00</point>
      <point r="2.180"> -0.30274517E+00  0.17797883E+00 -0.93679265E-02 -0.19896486E+00</point>
      <point r="2.200"> -0.30242341E+00  0.18085085E+00 -0.44880971E-02 -0.19720614E+00</point>
      <point r="2.220"> -0.30204006E+00  0.18363252E+00  0.32566929E-03 -0.19543540E+00</point>
      <point r="2.240"> -0.30159456E+00  0.18632231E+00  0.50699889E-02 -0.19365357E+00</point>
      <point r="2.260"> -0.30108653E+00  0.18891891E+00  0.97417370E-02 -0.19186160E+00</point>
      <point r="2.280"> -0.30051582E+00  0.19142128E+00  0.14338067E-01 -0.19006042E+00</point>
      <point r="2.300"> -0.29988240E+00  0.19382858E+00  0.18856424E-01 -0.18825092E+00</point>
      <point r="2.320"> -0.29918640E+00  0.19614019E+00  0.23294556E-01 -0.18643401E+00</point>
      <point r="2.340"> -0.29842808E+00  0.19835566E+00  0.27650514E-01 -0.18461056E+00</point>
      <point r="2.360"> -0.29760779E+00  0.20047472E+00  0.31922646E-01 -0.18278144E+00</point>
      <point r="2.380"> -0.29672602E+00  0.20249729E+00  0.36109588E-01 -0.18094748E+00</point>
      <point r="2.400"> -0.29578336E+00  0.20442342E+00  0.40210242E-01 -0.17910952E+00</point>
      <point r="2.420"> -0.29478047E+00  0.20625331E+00  0.44223758E-01 -0.17726836E+00</point>
      <point r="2.440"> -0.29371812E+00  0.20798731E+00  0.48149502E-01 -0.17542479E+00</point>
      <point r="2.460"> -0.29259717E+00  0.20962589E+00  0.51987038E-01 -0.17357958E+00</point>
      <point r="2.480"> -0.29141856E+00  0.21116966E+00  0.55736094E-01 -0.17173349E+00</point>
      <point r="2.500"> -0.29018331E+00  0.21261934E+00  0.59396544E-01 -0.16988724E+00</point>
      <point r="2.520"> -0.28889253E+00  0.21397575E+00  0.62968382E-01 -0.16804154E+00</point>
      <point r="2.540"> -0.28754739E+00  0.21523984E+00  0.66451711E-01 -0.16619709E+00</point>
      <point r="2.560"> -0.28614914E+00  0.21641263E+00  0.69846720E-01 -0.16435455E+00</point>
      <point r="2.580"> -0.28469908E+00  0.21749525E+00  0.73153681E-01 -0.16251458E+00</point>
      <point r="2.600"> -0.28319856E+00  0.21848890E+00  0.76372940E-01 -0.16067780E+00</point>
      <point r="2.620"> -0.28164900E+00  0.21939487E+00  0.79504912E-01 -0.15884481E+00</point>
      <point r="2.640"> -0.28005184E+00  0.22021450E+00  0.82550079E-01 -0.15701619E+00</point>
      <point r="2.660"> -0.27840856E+00  0.22094921E+00  0.85508991E-01 -0.15519250E+00</point>
      <point r="2.680"> -0.27672067E+00  0.22160045E+00  0.88382264E-01 -0.15337429E+00</point>
      <point r="2.700"> -0.27498968E+00  0.22216974E+00  0.91170584E-01 -0.15156205E+00</point>
      <point r="2.720"> -0.27321713E+00  0.22265862E+00  0.93874703E-01 -0.14975630E+00</point>
      <point r="2.740"> -0.27140457E+00  0.22306867E+00  0.96495441E-01 -0.14795748E+00</point>
      <point r="2.760"> -0.26955354E+00  0.22340149E+00  0.99033684E-01 -0.14616607E+00</point>
      <point r="2.780"> -0.26766559E+00  0.22365873E+00  0.10149038E+00 -0.14438248E+00</point>
      <point r="2.800"> -0.26574224E+00  0.22384203E+00  0.10386654E+00 -0.14260713E+00</point>
      <point r="2.820"> -0.26378504E+00  0.22395305E+00  0.10616323E+00 -0.14084041E+00</point>
      <point r="2.840"> -0.26179551E+00  0.22399347E+00  0.10838157E+00 -0.13908270E+00</point>
      <point r="2.860"> -0.25977516E+00  0.22396498E+00  0.11052270E+00 -0.13733433E+00</point>
      <point r="2.880"> -0.25772548E+00  0.22386926E+00  0.11258782E+00 -0.13559566E+00</point>
      <point r="2.900"> -0.25564796E+00  0.22370800E+00  0.11457817E+00 -0.13386700E+00</point>
      <point r="2.920"> -0.25354408E+00  0.22348290E+00  0.11649497E+00 -0.13214865E+00</point>
      <point r="2.940"> -0.25141527E+00  0.22319563E+00  0.11833948E+00 -0.13044089E+00</point>
      <point r="2.960"> -0.24926298E+00  0.22284788E+00  0.12011298E+00 -0.12874401E+00</point>
      <point r="2.980"> -0.24708863E+00  0.22244133E+00  0.12181672E+00 -0.12705826E+00</point>
      <point r="3.000"> -0.24489360E+00  0.22197762E+00  0.12345197E+00 -0.12538387E+00</point>
      <point r="3.020"> -0.24267927E+00  0.22145842E+00  0.12501998E+00 -0.12372107E+00</point>
      <point r="3.040"> -0.24044699E+00  0.22088534E+00  0.12652199E+00 -0.12207008E+00</point>
      <point r="3.060"> -0.23819808E+00  0.22026002E+00  0.12795924E+00 -0.12043109E+00</point>
      <point r="3.080"> -0.23593386E+00  0.21958405E+00  0.12933294E+00 -0.11880428E+00</point>
      <point r="3.100"> -0.23365559E+00  0.21885900E+00  0.13064428E+00 -0.11718985E+00</point>
      <point r="3.120"> -0.23136453E+00  0.21808645E+00  0.13189445E+00 -0.11558793E+00</point>
      <point r="3.140"> -0.22906190E+00  0.21726794E+00  0.13308462E+00 -0.11399868E+00</point>
      <point r="3.160"> -0.22674890E+00  0.21640497E+00  0.13421595E+00 -0.11242225E+00</point>
      <point r="3.180"> -0.22442670E+00  0.21549905E+00  0.13528956E+00 -0.11085875E+00</point>
      <point r="3.200"> -0.22209645E+00  0.21455166E+00  0.13630660E+00 -0.10930831E+00</point>
      <point r="3.220"> -0.21975925E+00  0.21356425E+00  0.13726816E+00 -0.10777102E+00</point>
      <point r="3.240"> -0.21741619E+00  0.21253823E+00  0.13817536E+00 -0.10624699E+00</point>
      <point r="3.260"> -0.21506834E+00  0.21147503E+00  0.13902927E+00 -0.10473631E+00</point>
      <point r="3.280"> -0.21271671E+00  0.21037600E+00  0.13983096E+00 -0.10323905E+00</point>
      <point r="3.300"> -0.21036230E+00  0.20924251E+00  0.14058150E+00 -0.10175528E+00</point>
      <point r="3.320"> -0.20800608E+00  0.20807588E+00  0.14128191E+00 -0.10028506E+00</point>
      <point r="3.340"> -0.20564863E+00  0.20687662E+00  0.14193215E+00 -0.98828579E-01</point>
      <point r="3.360"> -0.20329203E+00  0.20564822E+00  0.14253607E+00 -0.97385638E-01</point>
      <point r="3.380"> -0.20093634E+00  0.20439048E+00  0.14309289E+00 -0.95956387E-01</point>
      <point r="3.400"> -0.19858241E+00  0.20310459E+00  0.14360356E+00 -0.94540861E-01</point>
      <point r="3.420"> -0.19623103E+00  0.20179171E+00  0.14406897E+00 -0.93139087E-01</point>
      <point r="3.440"> -0.19388299E+00  0.20045298E+00  0.14449003E+00 -0.91751086E-01</point>
      <point r="3.460"> -0.19153905E+00  0.19908952E+00  0.14486763E+00 -0.90376874E-01</point>
      <point r="3.480"> -0.18919996E+00  0.19770244E+00  0.14520263E+00 -0.89016462E-01</point>
      <point r="3.500"> -0.18686643E+00  0.19629281E+00  0.14549590E+00 -0.87669855E-01</point>
      <point r="3.520"> -0.18453916E+00  0.19486170E+00  0.14574830E+00 -0.86337052E-01</point>
      <point r="3.540"> -0.18221883E+00  0.19341015E+00  0.14596067E+00 -0.85018051E-01</point>
      <point r="3.560"> -0.17990609E+00  0.19193919E+00  0.14613385E+00 -0.83712841E-01</point>
      <point r="3.580"> -0.17760159E+00  0.19044984E+00  0.14626866E+00 -0.82421410E-01</point>
      <point r="3.600"> -0.17530593E+00  0.18894309E+00  0.14636592E+00 -0.81143739E-01</point>
      <point r="3.620"> -0.17301972E+00  0.18741989E+00  0.14642644E+00 -0.79879807E-01</point>
      <point r="3.640"> -0.17074352E+00  0.18588122E+00  0.14645099E+00 -0.78629589E-01</point>
      <point r="3.660"> -0.16847788E+00  0.18432798E+00  0.14644037E+00 -0.77393054E-01</point>
      <point r="3.680"> -0.16622334E+00  0.18276108E+00  0.14639532E+00 -0.76170170E-01</point>
      <point r="3.700"> -0.16398041E+00  0.18118142E+00  0.14631660E+00 -0.74960901E-01</point>
      <point r="3.720"> -0.16174956E+00  0.17958984E+00  0.14620493E+00 -0.73765206E-01</point>
      <point r="3.740"> -0.15953127E+00  0.17798718E+00  0.14606102E+00 -0.72583045E-01</point>
      <point r="3.760"> -0.15732596E+00  0.17637426E+00  0.14588556E+00 -0.71414370E-01</point>
      <point r="3.780"> -0.15513407E+00  0.17475184E+00  0.14567923E+00 -0.70259135E-01</point>
      <point r="3.800"> -0.15295599E+00  0.17312071E+00  0.14544268E+00 -0.69117288E-01</point>
      <point r="3.820"> -0.15079209E+00  0.17148159E+00  0.14517656E+00 -0.67988776E-01</point>
      <point r="3.840"> -0.14864273E+00  0.16983520E+00  0.14488149E+00 -0.66873544E-01</point>
      <point r="3.860"> -0.14650826E+00  0.16818224E+00  0.14455807E+00 -0.65771536E-01</point>
      <point r="3.880"> -0.14438898E+00  0.16652338E+00  0.14420690E+00 -0.64682690E-01</point>
      <point r="3.900"> -0.14228521E+00  0.16485928E+00  0.14382856E+00 -0.63606947E-01</point>
      <point r="3.920"> -0.14019722E+00  0.16319056E+00  0.14342363E+00 -0.62544242E-01</point>
      <point r="3.940"> -0.13812529E+00  0.16151786E+00  0.14299267E+00 -0.61494512E-01</point>
      <point r="3.960"> -0.13606967E+00  0.15984176E+00  0.14253624E+00 -0.60457690E-01</point>
      <point r="3.980"> -0.13403059E+00  0.15816287E+00  0.14205487E+00 -0.59433708E-01</point>
      <point r="4.000"> -0.13200829E+00  0.15648174E+00  0.14154912E+00 -0.58422498E-01</point>
      <point r="4.020"> -0.13000298E+00  0.15479893E+00  0.14101953E+00 -0.57423988E-01</point>
      <point r="4.040"> -0.12801485E+00  0.15311500E+00  0.14046662E+00 -0.56438106E-01</point>
      <point r="4.060"> -0.12604411E+00  0.15143046E+00  0.13989094E+00 -0.55464780E-01</point>
      <point r="4.080"> -0.12409091E+00  0.14974586E+00  0.13929301E+00 -0.54503935E-01</point>
      <point r="4.100"> -0.12215544E+00  0.14806168E+00  0.13867336E+00 -0.53555495E-01</point>
      <point r="4.120"> -0.12023785E+00  0.14637844E+00  0.13803251E+00 -0.52619385E-01</point>
      <point r="4.140"> -0.11833828E+00  0.14469661E+00  0.13737101E+00 -0.51695526E-01</point>
      <point r="4.160"> -0.11645688E+00  0.14301668E+00  0.13668936E+00 -0.50783839E-01</point>
      <point r="4.180"> -0.11459375E+00  0.14133910E+00  0.13598810E+00 -0.49884246E-01</point>
      <point r="4.200"> -0.11274903E+00  0.13966434E+00  0.13526775E+00 -0.48996664E-01</point>
      <point r="4.220"> -0.11092282E+00  0.13799284E+00  0.13452882E+00 -0.48121014E-01</point>
      <point r="4.240"> -0.10911522E+00  0.13632502E+00  0.13377184E+00 -0.47257211E-01</point>
      <point r="4.260"> -0.10732631E+00  0.13466132E+00  0.13299732E+00 -0.46405172E-01</point>
      <point r="4.280"> -0.10555617E+00  0.13300213E+00  0.13220578E+00 -0.45564814E-01</point>
      <point r="4.300"> -0.10380488E+00  0.13134787E+00  0.13139771E+00 -0.44736051E-01</point>
      <point r="4.320"> -0.10207249E+00  0.12969891E+00  0.13057363E+00 -0.43918798E-01</point>
      <point r="4.340"> -0.10035905E+00  0.12805564E+00  0.12973404E+00 -0.43112967E-01</point>
      <point r="4.360"> -0.98664602E-01  0.12641841E+00  0.12887943E+00 -0.42318472E-01</point>
      <point r="4.380"> -0.96989187E-01  0.12478759E+00  0.12801029E+00 -0.41535224E-01</point>
      <point r="4.400"> -0.95332826E-01  0.12316351E+00  0.12712711E+00 -0.40763135E-01</point>
      <point r="4.420"> -0.93695537E-01  0.12154650E+00  0.12623038E+00 -0.40002116E-01</point>
      <point r="4.440"> -0.92077328E-01  0.11993689E+00  0.12532057E+00 -0.39252078E-01</point>
      <point r="4.460"> -0.90478201E-01  0.11833498E+00  0.12439815E+00 -0.38512929E-01</point>
      <point r="4.480"> -0.88898150E-01  0.11674107E+00  0.12346359E+00 -0.37784578E-01</point>
      <point r="4.500"> -0.87337161E-01  0.11515546E+00  0.12251736E+00 -0.37066936E-01</point>
      <point r="4.520"> -0.85795214E-01  0.11357841E+00  0.12155991E+00 -0.36359910E-01</point>
      <point r="4.540"> -0.84272282E-01  0.11201021E+00  0.12059171E+00 -0.35663407E-01</point>
      <point r="4.560"> -0.82768332E-01  0.11045110E+00  0.11961319E+00 -0.34977335E-01</point>
      <point r="4.580"> -0.81283323E-01  0.10890134E+00  0.11862482E+00 -0.34301602E-01</point>
      <point r="4.600"> -0.79817212E-01  0.10736117E+00  0.11762703E+00 -0.33636114E-01</point>
      <point r="4.620"> -0.78369946E-01  0.10583082E+00  0.11662027E+00 -0.32980777E-01</point>
      <point r="4.640"> -0.76941470E-01  0.10431051E+00  0.11560496E+00 -0.32335498E-01</point>
      <point r="4.660"> -0.75531722E-01  0.10280047E+00  0.11458156E+00 -0.31700182E-01</point>
      <point r="4.680"> -0.74140635E-01  0.10130090E+00  0.11355048E+00 -0.31074734E-01</point>
      <point r="4.700"> -0.72768139E-01  0.99811994E-01  0.11251216E+00 -0.30459061E-01</point>
      <point r="4.720"> -0.71414158E-01  0.98333953E-01  0.11146702E+00 -0.29853067E-01</point>
      <point r="4.740"> -0.70078613E-01  0.96866960E-01  0.11041549E+00 -0.29256657E-01</point>
      <point r="4.760"> -0.68761418E-01  0.95411192E-01  0.10935798E+00 -0.28669736E-01</point>
      <point r="4.780"> -0.67462487E-01  0.93966819E-01  0.10829490E+00 -0.28092209E-01</point>
      <point r="4.800"> -0.66181727E-01  0.92534006E-01  0.10722667E+00 -0.27523979E-01</point>
      <point r="4.820"> -0.64919043E-01  0.91112907E-01  0.10615370E+00 -0.26964951E-01</point>
      <point r="4.840"> -0.63674336E-01  0.89703674E-01  0.10507639E+00 -0.26415029E-01</point>
      <point r="4.860"> -0.62447502E-01  0.88306448E-01  0.10399515E+00 -0.25874117E-01</point>
      <point r="4.880"> -0.61238435E-01  0.86921364E-01  0.10291036E+00 -0.25342120E-01</point>
      <point r="4.900"> -0.60047027E-01  0.85548551E-01  0.10182242E+00 -0.24818941E-01</point>
      <point r="4.920"> -0.58873165E-01  0.84188132E-01  0.10073173E+00 -0.24304485E-01</point>
      <point r="4.940"> -0.57716732E-01  0.82840221E-01  0.99638651E-01 -0.23798654E-01</point>
      <point r="4.960"> -0.56577611E-01  0.81504928E-01  0.98543576E-01 -0.23301354E-01</point>
      <point r="4.980"> -0.55455680E-01  0.80182353E-01  0.97446873E-01 -0.22812488E-01</point>
      <point r="5.000"> -0.54350814E-01  0.78872593E-01  0.96348911E-01 -0.22331960E-01</point>
      <point r="5.020"> -0.53262887E-01  0.77575735E-01  0.95250050E-01 -0.21859674E-01</point>
      <point r="5.040"> -0.52191769E-01  0.76291863E-01  0.94150648E-01 -0.21395535E-01</point>
      <point r="5.060"> -0.51137327E-01  0.75021053E-01  0.93051054E-01 -0.20939447E-01</point>
      <point r="5.080"> -0.50099429E-01  0.73763373E-01  0.91951615E-01 -0.20491315E-01</point>
      <point r="5.100"> -0.49077936E-01  0.72518889E-01  0.90852667E-01 -0.20051042E-01</point>
      <point r="5.120"> -0.48072711E-01  0.71287656E-01  0.89754545E-01 -0.19618534E-01</point>
      <point r="5.140"> -0.47083612E-01  0.70069726E-01  0.88657573E-01 -0.19193697E-01</point>
      <point r="5.160"> -0.46110496E-01  0.68865146E-01  0.87562073E-01 -0.18776434E-01</point>
      <point r="5.180"> -0.45153219E-01  0.67673954E-01  0.86468359E-01 -0.18366652E-01</point>
      <point r="5.200"> -0.44211635E-01  0.66496184E-01  0.85376736E-01 -0.17964256E-01</point>
      <point r="5.220"> -0.43285595E-01  0.65331865E-01  0.84287507E-01 -0.17569153E-01</point>
      <point r="5.240"> -0.42374950E-01  0.64181020E-01  0.83200965E-01 -0.17181248E-01</point>
      <point r="5.260"> -0.41479548E-01  0.63043666E-01  0.82117400E-01 -0.16800449E-01</point>
      <point r="5.280"> -0.40599238E-01  0.61919815E-01  0.81037092E-01 -0.16426662E-01</point>
      <point r="5.300"> -0.39733867E-01  0.60809476E-01  0.79960316E-01 -0.16059795E-01</point>
      <point r="5.320"> -0.38883278E-01  0.59712649E-01  0.78887341E-01 -0.15699755E-01</point>
      <point r="5.340"> -0.38047318E-01  0.58629333E-01  0.77818429E-01 -0.15346450E-01</point>
      <point r="5.360"> -0.37225829E-01  0.57559519E-01  0.76753834E-01 -0.14999790E-01</point>
      <point r="5.380"> -0.36418654E-01  0.56503197E-01  0.75693807E-01 -0.14659682E-01</point>
      <point r="5.400"> -0.35625634E-01  0.55460348E-01  0.74638589E-01 -0.14326036E-01</point>
      <point r="5.420"> -0.34846612E-01  0.54430953E-01  0.73588415E-01 -0.13998763E-01</point>
      <point r="5.440"> -0.34081427E-01  0.53414986E-01  0.72543516E-01 -0.13677771E-01</point>
      <point r="5.460"> -0.33329919E-01  0.52412418E-01  0.71504113E-01 -0.13362973E-01</point>
      <point r="5.480"> -0.32591929E-01  0.51423214E-01  0.70470423E-01 -0.13054279E-01</point>
      <point r="5.500"> -0.31867296E-01  0.50447338E-01  0.69442655E-01 -0.12751602E-01</point>
      <point r="5.520"> -0.31155859E-01  0.49484747E-01  0.68421012E-01 -0.12454853E-01</point>
      <point r="5.540"> -0.30457456E-01  0.48535396E-01  0.67405690E-01 -0.12163945E-01</point>
      <point r="5.560"> -0.29771926E-01  0.47599237E-01  0.66396880E-01 -0.11878792E-01</point>
      <point r="5.580"> -0.29099108E-01  0.46676217E-01  0.65394764E-01 -0.11599308E-01</point>
      <point r="5.600"> -0.28438841E-01  0.45766280E-01  0.64399520E-01 -0.11325407E-01</point>
      <point r="5.620"> -0.27790964E-01  0.44869367E-01  0.63411318E-01 -0.11057005E-01</point>
      <point r="5.640"> -0.27155314E-01  0.43985414E-01  0.62430323E-01 -0.10794018E-01</point>
      <point r="5.660"> -0.26531732E-01  0.43114356E-01  0.61456691E-01 -0.10536361E-01</point>
      <point r="5.680"> -0.25920056E-01  0.42256125E-01  0.60490574E-01 -0.10283953E-01</point>
      <point r="5.700"> -0.25320127E-01  0.41410647E-01  0.59532117E-01 -0.10036710E-01</point>
      <point r="5.720"> -0.24731783E-01  0.40577849E-01  0.58581458E-01 -0.97945517E-02</point>
      <point r="5.740"> -0.24154865E-01  0.39757652E-01  0.57638731E-01 -0.95573963E-02</point>
      <point r="5.760"> -0.23589215E-01  0.38949977E-01  0.56704061E-01 -0.93251639E-02</point>
      <point r="5.780"> -0.23034673E-01  0.38154740E-01  0.55777568E-01 -0.90977749E-02</point>
      <point r="5.800"> -0.22491080E-01  0.37371856E-01  0.54859367E-01 -0.88751504E-02</point>
      <point r="5.820"> -0.21958281E-01  0.36601237E-01  0.53949565E-01 -0.86572122E-02</point>
      <point r="5.840"> -0.21436117E-01  0.35842794E-01  0.53048265E-01 -0.84438829E-02</point>
      <point r="5.860"> -0.20924433E-01  0.35096434E-01  0.52155564E-01 -0.82350859E-02</point>
      <point r="5.880"> -0.20423073E-01  0.34362061E-01  0.51271551E-01 -0.80307450E-02</point>
      <point r="5.900"> -0.19931882E-01  0.33639581E-01  0.50396311E-01 -0.78307850E-02</point>
      <point r="5.920"> -0.19450707E-01  0.32928894E-01  0.49529924E-01 -0.76351314E-02</point>
      <point r="5.940"> -0.18979394E-01  0.32229901E-01  0.48672463E-01 -0.74437105E-02</point>
      <point r="5.960"> -0.18517792E-01  0.31542499E-01  0.47823996E-01 -0.72564492E-02</point>
      <point r="5.980"> -0.18065750E-01  0.30866585E-01  0.46984587E-01 -0.70732753E-02</point>
      <point r="6.000"> -0.17623117E-01  0.30202053E-01  0.46154293E-01 -0.68941173E-02</point>
      <point r="6.020"> -0.17189745E-01  0.29548797E-01  0.45333166E-01 -0.67189046E-02</point>
      <point r="6.040"> -0.16765485E-01  0.28906709E-01  0.44521253E-01 -0.65475672E-02</point>
      <point r="6.060"> -0.16350190E-01  0.28275679E-01  0.43718596E-01 -0.63800359E-02</point>
      <point r="6.080"> -0.15943715E-01  0.27655597E-01  0.42925234E-01 -0.62162424E-02</point>
      <point r="6.100"> -0.15545915E-01  0.27046350E-01  0.42141199E-01 -0.60561190E-02</point>
      <point r="6.120"> -0.15156646E-01  0.26447827E-01  0.41366517E-01 -0.58995992E-02</point>
      <point r="6.140"> -0.14775767E-01  0.25859913E-01  0.40601213E-01 -0.57466167E-02</point>
      <point r="6.160"> -0.14403135E-01  0.25282494E-01  0.39845306E-01 -0.55971065E-02</point>
      <point r="6.180"> -0.14038612E-01  0.24715454E-01  0.39098808E-01 -0.54510040E-02</point>
      <point r="6.200"> -0.13682059E-01  0.24158676E-01  0.38361732E-01 -0.53082458E-02</point>
      <point r="6.220"> -0.13333338E-01  0.23612044E-01  0.37634081E-01 -0.51687690E-02</point>
      <point r="6.240"> -0.12992313E-01  0.23075440E-01  0.36915858E-01 -0.50325116E-02</point>
      <point r="6.260"> -0.12658850E-01  0.22548746E-01  0.36207060E-01 -0.48994123E-02</point>
      <point r="6.280"> -0.12332816E-01  0.22031843E-01  0.35507681E-01 -0.47694108E-02</point>
      <point r="6.300"> -0.12014079E-01  0.21524613E-01  0.34817710E-01 -0.46424474E-02</point>
      <point r="6.320"> -0.11702507E-01  0.21026935E-01  0.34137135E-01 -0.45184633E-02</point>
      <point r="6.340"> -0.11397973E-01  0.20538691E-01  0.33465937E-01 -0.43974005E-02</point>
      <point r="6.360"> -0.11100348E-01  0.20059760E-01  0.32804095E-01 -0.42792018E-02</point>
      <point r="6.380"> -0.10809506E-01  0.19590022E-01  0.32151585E-01 -0.41638107E-02</point>
      <point r="6.400"> -0.10525322E-01  0.19129357E-01  0.31508379E-01 -0.40511715E-02</point>
      <point r="6.420"> -0.10247673E-01  0.18677646E-01  0.30874446E-01 -0.39412295E-02</point>
      <point r="6.440"> -0.99764366E-02  0.18234767E-01  0.30249753E-01 -0.38339306E-02</point>
      <point r="6.460"> -0.97114924E-02  0.17800601E-01  0.29634261E-01 -0.37292214E-02</point>
      <point r="6.480"> -0.94527213E-02  0.17375028E-01  0.29027932E-01 -0.36270494E-02</point>
      <point r="6.500"> -0.92000058E-02  0.16957928E-01  0.28430722E-01 -0.35273630E-02</point>
      <point r="6.520"> -0.89532296E-02  0.16549181E-01  0.27842586E-01 -0.34301112E-02</point>
      <point r="6.540"> -0.87122783E-02  0.16148669E-01  0.27263475E-01 -0.33352439E-02</point>
      <point r="6.560"> -0.84770386E-02  0.15756272E-01  0.26693340E-01 -0.32427115E-02</point>
      <point r="6.580"> -0.82473989E-02  0.15371872E-01  0.26132126E-01 -0.31524654E-02</point>
      <point r="6.600"> -0.80232490E-02  0.14995350E-01  0.25579779E-01 -0.30644578E-02</point>
      <point r="6.620"> -0.78044802E-02  0.14626590E-01  0.25036240E-01 -0.29786416E-02</point>
      <point r="6.640"> -0.75909854E-02  0.14265473E-01  0.24501451E-01 -0.28949702E-02</point>
      <point r="6.660"> -0.73826587E-02  0.13911883E-01  0.23975349E-01 -0.28133982E-02</point>
      <point r="6.680"> -0.71795461E-02  0.13566036E-01  0.23458599E-01 -0.27339227E-02</point>
      <point r="6.700"> -0.69812502E-02  0.13227165E-01  0.22949708E-01 -0.26564159E-02</point>
      <point r="6.720"> -0.67878140E-02  0.12895475E-01  0.22449307E-01 -0.25808757E-02</point>
      <point r="6.740"> -0.65991378E-02  0.12570851E-01  0.21957326E-01 -0.25072595E-02</point>
      <point r="6.760"> -0.64151231E-02  0.12253180E-01  0.21473693E-01 -0.24355252E-02</point>
      <point r="6.780"> -0.62356733E-02  0.11942348E-01  0.20998337E-01 -0.23656315E-02</point>
      <point r="6.800"> -0.60606929E-02  0.11638245E-01  0.20531183E-01 -0.22975379E-02</point>
      <point r="6.820"> -0.58900879E-02  0.11340758E-01  0.20072156E-01 -0.22312045E-02</point>
      <point r="6.840"> -0.57237660E-02  0.11049777E-01  0.19621179E-01 -0.21665921E-02</point>
      <point r="6.860"> -0.55616361E-02  0.10765194E-01  0.19178175E-01 -0.21036621E-02</point>
      <point r="6.880"> -0.54036086E-02  0.10486899E-01  0.18743065E-01 -0.20423767E-02</point>
      <point r="6.900"> -0.52495954E-02  0.10214785E-01  0.18315768E-01 -0.19826987E-02</point>
      <point r="6.920"> -0.50995097E-02  0.99487455E-02  0.17896203E-01 -0.19245918E-02</point>
      <point r="6.940"> -0.49532663E-02  0.96886740E-02  0.17484289E-01 -0.18680199E-02</point>
      <point r="6.960"> -0.48107811E-02  0.94344660E-02  0.17079942E-01 -0.18129480E-02</point>
      <point r="6.980"> -0.46719718E-02  0.91860177E-02  0.16683079E-01 -0.17593414E-02</point>
      <point r="7.000"> -0.45367571E-02  0.89432263E-02  0.16293615E-01 -0.17071664E-02</point>
      <point r="7.020"> -0.44050573E-02  0.87059899E-02  0.15911466E-01 -0.16563896E-02</point>
      <point r="7.040"> -0.42767941E-02  0.84742079E-02  0.15536545E-01 -0.16069783E-02</point>
      <point r="7.060"> -0.41518903E-02  0.82477806E-02  0.15168766E-01 -0.15589006E-02</point>
      <point r="7.080"> -0.40302704E-02  0.80266092E-02  0.14808043E-01 -0.15121251E-02</point>
      <point r="7.100"> -0.39118598E-02  0.78105962E-02  0.14454289E-01 -0.14666208E-02</point>
      <point r="7.120"> -0.37965856E-02  0.75996453E-02  0.14107416E-01 -0.14223576E-02</point>
      <point r="7.140"> -0.36843761E-02  0.73936609E-02  0.13767336E-01 -0.13793059E-02</point>
      <point r="7.160"> -0.35751607E-02  0.71925489E-02  0.13433962E-01 -0.13374365E-02</point>
      <point r="7.180"> -0.34688704E-02  0.69962161E-02  0.13107205E-01 -0.12967210E-02</point>
      <point r="7.200"> -0.33654372E-02  0.68045706E-02  0.12786977E-01 -0.12571316E-02</point>
      <point r="7.220"> -0.32647945E-02  0.66175214E-02  0.12473189E-01 -0.12186407E-02</point>
      <point r="7.240"> -0.31668768E-02  0.64349788E-02  0.12165754E-01 -0.11812215E-02</point>
      <point r="7.260"> -0.30716201E-02  0.62568543E-02  0.11864583E-01 -0.11448479E-02</point>
      <point r="7.280"> -0.29789612E-02  0.60830603E-02  0.11569586E-01 -0.11094941E-02</point>
      <point r="7.300"> -0.28888386E-02  0.59135107E-02  0.11280677E-01 -0.10751348E-02</point>
      <point r="7.320"> -0.28011916E-02  0.57481203E-02  0.10997767E-01 -0.10417453E-02</point>
      <point r="7.340"> -0.27159608E-02  0.55868050E-02  0.10720767E-01 -0.10093016E-02</point>
      <point r="7.360"> -0.26330879E-02  0.54294822E-02  0.10449591E-01 -0.97777975E-03</point>
      <point r="7.380"> -0.25525158E-02  0.52760701E-02  0.10184150E-01 -0.94715672E-03</point>
      <point r="7.400"> -0.24741887E-02  0.51264882E-02  0.99243574E-02 -0.91740979E-03</point>
      <point r="7.420"> -0.23980515E-02  0.49806572E-02  0.96701262E-02 -0.88851673E-03</point>
      <point r="7.440"> -0.23240506E-02  0.48384989E-02  0.94213698E-02 -0.86045580E-03</point>
      <point r="7.460"> -0.22521332E-02  0.46999363E-02  0.91780020E-02 -0.83320572E-03</point>
      <point r="7.480"> -0.21822477E-02  0.45648933E-02  0.89399371E-02 -0.80674567E-03</point>
      <point r="7.500"> -0.21143436E-02  0.44332954E-02  0.87070897E-02 -0.78105528E-03</point>
      <point r="7.520"> -0.20483714E-02  0.43050689E-02  0.84793749E-02 -0.75611462E-03</point>
      <point r="7.540"> -0.19842825E-02  0.41801413E-02  0.82567084E-02 -0.73190422E-03</point>
      <point r="7.560"> -0.19220294E-02  0.40584413E-02  0.80390065E-02 -0.70840500E-03</point>
      <point r="7.580"> -0.18615657E-02  0.39398986E-02  0.78261857E-02 -0.68559834E-03</point>
      <point r="7.600"> -0.18028459E-02  0.38244443E-02  0.76181633E-02 -0.66346602E-03</point>
      <point r="7.620"> -0.17458253E-02  0.37120102E-02  0.74148573E-02 -0.64199023E-03</point>
      <point r="7.640"> -0.16904604E-02  0.36025296E-02  0.72161862E-02 -0.62115355E-03</point>
      <point r="7.660"> -0.16367085E-02  0.34959367E-02  0.70220691E-02 -0.60093896E-03</point>
      <point r="7.680"> -0.15845279E-02  0.33921669E-02  0.68324257E-02 -0.58132985E-03</point>
      <point r="7.700"> -0.15338777E-02  0.32911564E-02  0.66471766E-02 -0.56230996E-03</point>
      <point r="7.720"> -0.14847180E-02  0.31928429E-02  0.64662429E-02 -0.54386342E-03</point>
      <point r="7.740"> -0.14370096E-02  0.30971649E-02  0.62895465E-02 -0.52597472E-03</point>
      <point r="7.760"> -0.13907144E-02  0.30040621E-02  0.61170099E-02 -0.50862872E-03</point>
      <point r="7.780"> -0.13457951E-02  0.29134750E-02  0.59485565E-02 -0.49181062E-03</point>
      <point r="7.800"> -0.13022150E-02  0.28253456E-02  0.57841104E-02 -0.47550598E-03</point>
      <point r="7.820"> -0.12599385E-02  0.27396165E-02  0.56235964E-02 -0.45970070E-03</point>
      <point r="7.840"> -0.12189306E-02  0.26562315E-02  0.54669402E-02 -0.44438100E-03</point>
      <point r="7.860"> -0.11791573E-02  0.25751355E-02  0.53140682E-02 -0.42953345E-03</point>
      <point r="7.880"> -0.11405852E-02  0.24962742E-02  0.51649076E-02 -0.41514493E-03</point>
      <point r="7.900"> -0.11031818E-02  0.24195946E-02  0.50193864E-02 -0.40120263E-03</point>
      <point r="7.920"> -0.10669153E-02  0.23450444E-02  0.48774335E-02 -0.38769407E-03</point>
      <point r="7.940"> -0.10317544E-02  0.22725725E-02  0.47389787E-02 -0.37460705E-03</point>
      <point r="7.960"> -0.99766903E-03  0.22021285E-02  0.46039524E-02 -0.36192969E-03</point>
      <point r="7.980"> -0.96462940E-03  0.21336633E-02  0.44722860E-02 -0.34965039E-03</point>
      <point r="8.000"> -0.93260660E-03  0.20671285E-02  0.43439118E-02 -0.33775782E-03</point>
      <point r="8.020"> -0.90157239E-03  0.20024767E-02  0.42187629E-02 -0.32624098E-03</point>
      <point r="8.040"> -0.87149918E-03  0.19396614E-02  0.40967732E-02 -0.31508908E-03</point>
      <point r="8.060"> -0.84236006E-03  0.18786372E-02  0.39778777E-02 -0.30429166E-03</point>
      <point r="8.080"> -0.81412877E-03  0.18193594E-02  0.38620120E-02 -0.29383849E-03</point>
      <point r="8.100"> -0.78677967E-03  0.17617843E-02  0.37491127E-02 -0.28371959E-03</point>
      <point r="8.120"> -0.76028776E-03  0.17058689E-02  0.36391173E-02 -0.27392527E-03</point>
      <point r="8.140"> -0.73462866E-03  0.16515715E-02  0.35319642E-02 -0.26444606E-03</point>
      <point r="8.160"> -0.70977857E-03  0.15988507E-02  0.34275928E-02 -0.25527272E-03</point>
      <point r="8.180"> -0.68571430E-03  0.15476664E-02  0.33259431E-02 -0.24639629E-03</point>
      <point r="8.200"> -0.66241322E-03  0.14979793E-02  0.32269562E-02 -0.23780802E-03</point>
      <point r="8.220"> -0.63985329E-03  0.14497505E-02  0.31305741E-02 -0.22949936E-03</point>
      <point r="8.240"> -0.61801300E-03  0.14029426E-02  0.30367396E-02 -0.22146204E-03</point>
      <point r="8.260"> -0.59687140E-03  0.13575183E-02  0.29453965E-02 -0.21368795E-03</point>
      <point r="8.280"> -0.57640809E-03  0.13134417E-02  0.28564893E-02 -0.20616924E-03</point>
      <point r="8.300"> -0.55660314E-03  0.12706773E-02  0.27699638E-02 -0.19889823E-03</point>
      <point r="8.320"> -0.53743720E-03  0.12291905E-02  0.26857661E-02 -0.19186748E-03</point>
      <point r="8.340"> -0.51889136E-03  0.11889476E-02  0.26038438E-02 -0.18506972E-03</point>
      <point r="8.360"> -0.50094724E-03  0.11499153E-02  0.25241449E-02 -0.17849788E-03</point>
      <point r="8.380"> -0.48358694E-03  0.11120613E-02  0.24466185E-02 -0.17214509E-03</point>
      <point r="8.400"> -0.46679299E-03  0.10753542E-02  0.23712146E-02 -0.16600467E-03</point>
      <point r="8.420"> -0.45054844E-03  0.10397628E-02  0.22978840E-02 -0.16007010E-03</point>
      <point r="8.440"> -0.43483674E-03  0.10052570E-02  0.22265783E-02 -0.15433505E-03</point>
      <point r="8.460"> -0.41964181E-03  0.97180743E-03  0.21572502E-02 -0.14879337E-03</point>
      <point r="8.480"> -0.40494798E-03  0.93938512E-03  0.20898530E-02 -0.14343907E-03</point>
      <point r="8.500"> -0.39074003E-03  0.90796195E-03  0.20243410E-02 -0.13826633E-03</point>
      <point r="8.520"> -0.37700312E-03  0.87751042E-03  0.19606692E-02 -0.13326949E-03</point>
      <point r="8.540"> -0.36372283E-03  0.84800370E-03  0.18987935E-02 -0.12844304E-03</point>
      <point r="8.560"> -0.35088514E-03  0.81941556E-03  0.18386709E-02 -0.12378164E-03</point>
      <point r="8.580"> -0.33847639E-03  0.79172042E-03  0.17802587E-02 -0.11928009E-03</point>
      <point r="8.600"> -0.32648332E-03  0.76489331E-03  0.17235153E-02 -0.11493333E-03</point>
      <point r="8.620"> -0.31489303E-03  0.73890983E-03  0.16684001E-02 -0.11073645E-03</point>
      <point r="8.640"> -0.30369296E-03  0.71374622E-03  0.16148728E-02 -0.10668468E-03</point>
      <point r="8.660"> -0.29287093E-03  0.68937926E-03  0.15628944E-02 -0.10277340E-03</point>
      <point r="8.680"> -0.28241508E-03  0.66578633E-03  0.15124263E-02 -0.98998092E-04</point>
      <point r="8.700"> -0.27231389E-03  0.64294532E-03  0.14634308E-02 -0.95354389E-04</point>
      <point r="8.720"> -0.26255617E-03  0.62083472E-03  0.14158711E-02 -0.91838046E-04</point>
      <point r="8.740"> -0.25313104E-03  0.59943351E-03  0.13697108E-02 -0.88444938E-04</point>
      <point r="8.760"> -0.24402792E-03  0.57872123E-03  0.13249146E-02 -0.85171062E-04</point>
      <point r="8.780"> -0.23523656E-03  0.55867790E-03  0.12814477E-02 -0.82012533E-04</point>
      <point r="8.800"> -0.22674698E-03  0.53928407E-03  0.12392760E-02 -0.78965576E-04</point>
      <point r="8.820"> -0.21854949E-03  0.52052077E-03  0.11983664E-02 -0.76026531E-04</point>
      <point r="8.840"> -0.21063469E-03  0.50236951E-03  0.11586862E-02 -0.73191844E-04</point>
      <point r="8.860"> -0.20299344E-03  0.48481227E-03  0.11202034E-02 -0.70458066E-04</point>
      <point r="8.880"> -0.19561689E-03  0.46783151E-03  0.10828869E-02 -0.67821852E-04</point>
      <point r="8.900"> -0.18849640E-03  0.45141012E-03  0.10467060E-02 -0.65279956E-04</point>
      <point r="8.920"> -0.18162365E-03  0.43553145E-03  0.10116307E-02 -0.62829231E-04</point>
      <point r="8.940"> -0.17499050E-03  0.42017925E-03  0.97763191E-03 -0.60466621E-04</point>
      <point r="8.960"> -0.16858909E-03  0.40533774E-03  0.94468085E-03 -0.58189165E-04</point>
      <point r="8.980"> -0.16241179E-03  0.39099152E-03  0.91274953E-03 -0.55993993E-04</point>
      <point r="9.000"> -0.15645118E-03  0.37712560E-03  0.88181052E-03 -0.53878319E-04</point>
      <point r="9.020"> -0.15070007E-03  0.36372540E-03  0.85183703E-03 -0.51839445E-04</point>
      <point r="9.040"> -0.14515148E-03  0.35077671E-03  0.82280283E-03 -0.49874755E-04</point>
      <point r="9.060"> -0.13979866E-03  0.33826569E-03  0.79468229E-03 -0.47981712E-04</point>
      <point r="9.080"> -0.13463504E-03  0.32617891E-03  0.76745035E-03 -0.46157861E-04</point>
      <point r="9.100"> -0.12965425E-03  0.31450325E-03  0.74108251E-03 -0.44400818E-04</point>
      <point r="9.120"> -0.12485014E-03  0.30322599E-03  0.71555484E-03 -0.42708279E-04</point>
      <point r="9.140"> -0.12021671E-03  0.29233470E-03  0.69084393E-03 -0.41078008E-04</point>
      <point r="9.160"> -0.11574818E-03  0.28181735E-03  0.66692692E-03 -0.39507840E-04</point>
      <point r="9.180"> -0.11143892E-03  0.27166218E-03  0.64378148E-03 -0.37995681E-04</point>
      <point r="9.200"> -0.10728348E-03  0.26185779E-03  0.62138578E-03 -0.36539499E-04</point>
      <point r="9.220"> -0.10327658E-03  0.25239308E-03  0.59971851E-03 -0.35137331E-04</point>
      <point r="9.240"> -0.99413121E-04  0.24325726E-03  0.57875886E-03 -0.33787272E-04</point>
      <point r="9.260"> -0.95688137E-04  0.23443983E-03  0.55848649E-03 -0.32487483E-04</point>
      <point r="9.280"> -0.92096832E-04  0.22593060E-03  0.53888158E-03 -0.31236179E-04</point>
      <point r="9.300"> -0.88634556E-04  0.21771964E-03  0.51992474E-03 -0.30031638E-04</point>
      <point r="9.320"> -0.85296809E-04  0.20979733E-03  0.50159706E-03 -0.28872189E-04</point>
      <point r="9.340"> -0.82079231E-04  0.20215429E-03  0.48388008E-03 -0.27756219E-04</point>
      <point r="9.360"> -0.78977601E-04  0.19478144E-03  0.46675581E-03 -0.26682166E-04</point>
      <point r="9.380"> -0.75987835E-04  0.18766992E-03  0.45020668E-03 -0.25648520E-04</point>
      <point r="9.400"> -0.73105977E-04  0.18081115E-03  0.43421553E-03 -0.24653819E-04</point>
      <point r="9.420"> -0.70328201E-04  0.17419679E-03  0.41876566E-03 -0.23696652E-04</point>
      <point r="9.440"> -0.67650803E-04  0.16781875E-03  0.40384077E-03 -0.22775653E-04</point>
      <point r="9.460"> -0.65070201E-04  0.16166915E-03  0.38942496E-03 -0.21889503E-04</point>
      <point r="9.480"> -0.62582928E-04  0.15574037E-03  0.37550274E-03 -0.21036926E-04</point>
      <point r="9.500"> -0.60185633E-04  0.15002500E-03  0.36205899E-03 -0.20216690E-04</point>
      <point r="9.520"> -0.57875074E-04  0.14451583E-03  0.34907901E-03 -0.19427604E-04</point>
      <point r="9.540"> -0.55648115E-04  0.13920590E-03  0.33654845E-03 -0.18668517E-04</point>
      <point r="9.560"> -0.53501727E-04  0.13408842E-03  0.32445333E-03 -0.17938318E-04</point>
      <point r="9.580"> -0.51432980E-04  0.12915684E-03  0.31278005E-03 -0.17235933E-04</point>
      <point r="9.600"> -0.49439043E-04  0.12440478E-03  0.30151535E-03 -0.16560327E-04</point>
      <point r="9.620"> -0.47517179E-04  0.11982606E-03  0.29064631E-03 -0.15910498E-04</point>
      <point r="9.640"> -0.45664747E-04  0.11541468E-03  0.28016039E-03 -0.15285482E-04</point>
      <point r="9.660"> -0.43879192E-04  0.11116484E-03  0.27004533E-03 -0.14684344E-04</point>
      <point r="9.680"> -0.42158048E-04  0.10707092E-03  0.26028924E-03 -0.14106187E-04</point>
      <point r="9.700"> -0.40498934E-04  0.10312744E-03  0.25088053E-03 -0.13550141E-04</point>
      <point r="9.720"> -0.38899551E-04  0.99329129E-04  0.24180793E-03 -0.13015368E-04</point>
      <point r="9.740"> -0.37357678E-04  0.95670858E-04  0.23306049E-03 -0.12501062E-04</point>
      <point r="9.760"> -0.35871173E-04  0.92147665E-04  0.22462754E-03 -0.12006442E-04</point>
      <point r="9.780"> -0.34437969E-04  0.88754744E-04  0.21649872E-03 -0.11530758E-04</point>
      <point r="9.800"> -0.33056071E-04  0.85487441E-04  0.20866394E-03 -0.11073284E-04</point>
      <point r="9.820"> -0.31723554E-04  0.82341250E-04  0.20111342E-03 -0.10633322E-04</point>
      <point r="9.840"> -0.30438562E-04  0.79311808E-04  0.19383764E-03 -0.10210199E-04</point>
      <point r="9.860"> -0.29199305E-04  0.76394895E-04  0.18682735E-03 -0.98032669E-05</point>
      <point r="9.880"> -0.28004058E-04  0.73586425E-04  0.18007356E-03 -0.94118995E-05</point>
      <point r="9.900"> -0.26851156E-04  0.70882446E-04  0.17356756E-03 -0.90354950E-05</point>
      <point r="9.920"> -0.25738997E-04  0.68279135E-04  0.16730087E-03 -0.86734730E-05</point>
      <point r="9.940"> -0.24666034E-04  0.65772795E-04  0.16126526E-03 -0.83252750E-05</point>
      <point r="9.960"> -0.23630780E-04  0.63359849E-04  0.15545276E-03 -0.79903627E-05</point>
      <point r="9.980"> -0.22631799E-04  0.61036842E-04  0.14985561E-03 -0.76682183E-05</point>
      <point r="10.000"> -0.21667712E-04  0.58800432E-04  0.14446630E-03 -0.73583429E-05</point>
      <point r="10.020"> -0.20737188E-04  0.56647390E-04  0.13927753E-03 -0.70602567E-05</point>
      <point r="10.040"> -0.19838946E-04  0.54574595E-04  0.13428224E-03 -0.67734977E-05</point>
      <point r="10.060"> -0.18971755E-04  0.52579034E-04  0.12947356E-03 -0.64976214E-05</point>
      <point r="10.080"> -0.18134427E-04  0.50657796E-04  0.12484484E-03 -0.62322003E-05</point>
      <point r="10.100"> -0.17325822E-04  0.48808068E-04  0.12038964E-03 -0.59768231E-05</point>
      <point r="10.120"> -0.16544842E-04  0.47027135E-04  0.11610171E-03 -0.57310942E-05</point>
      <point r="10.140"> -0.15790431E-04  0.45312378E-04  0.11197501E-03 -0.54946334E-05</point>
      <point r="10.160"> -0.15061574E-04  0.43661267E-04  0.10800366E-03 -0.52670750E-05</point>
      <point r="10.180"> -0.14357293E-04  0.42071360E-04  0.10418199E-03 -0.50480674E-05</point>
      <point r="10.200"> -0.13676651E-04  0.40540305E-04  0.10050450E-03 -0.48372729E-05</point>
      <point r="10.220"> -0.13018746E-04  0.39065828E-04  0.96965864E-04 -0.46343668E-05</point>
      <point r="10.240"> -0.12382711E-04  0.37645740E-04  0.93560938E-04 -0.44390372E-05</point>
      <point r="10.260"> -0.11767713E-04  0.36277929E-04  0.90284732E-04 -0.42509844E-05</point>
      <point r="10.280"> -0.11172953E-04  0.34960358E-04  0.87132422E-04 -0.40699207E-05</point>
      <point r="10.300"> -0.10597663E-04  0.33691066E-04  0.84099344E-04 -0.38955696E-05</point>
      <point r="10.320"> -0.10041105E-04  0.32468161E-04  0.81180984E-04 -0.37276658E-05</point>
      <point r="10.340"> -0.95025726E-05  0.31289823E-04  0.78372981E-04 -0.35659543E-05</point>
      <point r="10.360"> -0.89813852E-05  0.30154295E-04  0.75671118E-04 -0.34101908E-05</point>
      <point r="10.380"> -0.84768916E-05  0.29059888E-04  0.73071318E-04 -0.32601404E-05</point>
      <point r="10.400">  0.00000000E+00</point>
    </H_spline>
    <S_spline>
      <point r="0.000">  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.48000000E+01</point>
      <point r="0.020">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.040">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.060">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.080">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.100">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.120">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.140">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.160">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.180">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.200">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.220">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.240">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.260">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.280">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.300">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.320">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.340">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.360">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.380">  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01  0.10000000E+01</point>
      <point r="0.400">  0.93581115E+00 -0.15301420E+00  0.85914858E+00  0.94682025E+00</point>
      <point r="0.420">  0.93094071E+00 -0.16005805E+00  0.84869532E+00  0.94239667E+00</point>
      <point r="0.440">  0.92593997E+00 -0.16690118E+00  0.83821654E+00  0.93789942E+00</point>
      <point r="0.460">  0.92080925E+00 -0.17353582E+00  0.82773044E+00  0.93333742E+00</point>
      <point r="0.480">  0.91554992E+00 -0.17995667E+00  0.81725176E+00  0.92871878E+00</point>
      <point r="0.500">  0.91016443E+00 -0.18616076E+00  0.80679208E+00  0.92405082E+00</point>
      <point r="0.520">  0.90465630E+00 -0.19214726E+00  0.79636008E+00  0.91934016E+00</point>
      <point r="0.540">  0.89903007E+00 -0.19791727E+00  0.78596185E+00  0.91459267E+00</point>
      <point r="0.560">  0.89329120E+00 -0.20347362E+00  0.77560116E+00  0.90981358E+00</point>
      <point r="0.580">  0.88744595E+00 -0.20882066E+00  0.76527974E+00  0.90500750E+00</point>
      <point r="0.600">  0.88150136E+00 -0.21396408E+00  0.75499755E+00  0.90017843E+00</point>
      <point r="0.620">  0.87546505E+00 -0.21891068E+00  0.74475301E+00  0.89532984E+00</point>
      <point r="0.640">  0.86934517E+00 -0.22366821E+00  0.73454321E+00  0.89046470E+00</point>
      <point r="0.660">  0.86315024E+00 -0.22824519E+00  0.72436418E+00  0.88558551E+00</point>
      <point r="0.680">  0.85688911E+00 -0.23265076E+00  0.71421105E+00  0.88069434E+00</point>
      <point r="0.700">  0.85057079E+00 -0.23689449E+00  0.70407824E+00  0.87579288E+00</point>
      <point r="0.720">  0.84420441E+00 -0.24098629E+00  0.69395961E+00  0.87088246E+00</point>
      <point r="0.740">  0.83779909E+00 -0.24493625E+00  0.68384865E+00  0.86596412E+00</point>
      <point r="0.760">  0.83136390E+00 -0.24875456E+00  0.67373858E+00  0.86103859E+00</point>
      <point r="0.780">  0.82490776E+00 -0.25245137E+00  0.66362249E+00  0.85610635E+00</point>
      <point r="0.800">  0.81843939E+00 -0.25603673E+00  0.65349343E+00  0.85116767E+00</point>
      <point r="0.820">  0.81196727E+00 -0.25952052E+00  0.64334452E+00  0.84622263E+00</point>
      <point r="0.840">  0.80549954E+00 -0.26291235E+00  0.63316902E+00  0.84127112E+00</point>
      <point r="0.860">  0.79904403E+00 -0.26622153E+00  0.62296041E+00  0.83631290E+00</point>
      <point r="0.880">  0.79260818E+00 -0.26945701E+00  0.61271244E+00  0.83134761E+00</point>
      <point r="0.900">  0.78619902E+00 -0.27262734E+00  0.60241919E+00  0.82637477E+00</point>
      <point r="0.920">  0.77982317E+00 -0.27574065E+00  0.59207512E+00  0.82139383E+00</point>
      <point r="0.940">  0.77348680E+00 -0.27880461E+00  0.58167507E+00  0.81640417E+00</point>
      <point r="0.960">  0.76719564E+00 -0.28182643E+00  0.57121433E+00  0.81140513E+00</point>
      <point r="0.980">  0.76095495E+00 -0.28481282E+00  0.56068863E+00  0.80639598E+00</point>
      <point r="1.000">  0.75476954E+00 -0.28777001E+00  0.55009417E+00  0.80137601E+00</point>
      <point r="1.020">  0.74864378E+00 -0.29070373E+00  0.53942763E+00  0.79634446E+00</point>
      <point r="1.040">  0.74258157E+00 -0.29361921E+00  0.52868617E+00  0.79130059E+00</point>
      <point r="1.060">  0.73658638E+00 -0.29652119E+00  0.51786742E+00  0.78624366E+00</point>
      <point r="1.080">  0.73066125E+00 -0.29941393E+00  0.50696953E+00  0.78117295E+00</point>
      <point r="1.100">  0.72480880E+00 -0.30230122E+00  0.49599108E+00  0.77608775E+00</point>
      <point r="1.120">  0.71903124E+00 -0.30518637E+00  0.48493115E+00  0.77098740E+00</point>
      <point r="1.140">  0.71333040E+00 -0.30807225E+00  0.47378926E+00  0.76587125E+00</point>
      <point r="1.160">  0.70770772E+00 -0.31096129E+00  0.46256539E+00  0.76073871E+00</point>
      <point r="1.180">  0.70216429E+00 -0.31385552E+00  0.45125997E+00  0.75558923E+00</point>
      <point r="1.200">  0.69670088E+00 -0.31675655E+00  0.43987381E+00  0.75042231E+00</point>
      <point r="1.220">  0.69131791E+00 -0.31966561E+00  0.42840815E+00  0.74523747E+00</point>
      <point r="1.240">  0.68601552E+00 -0.32258356E+00  0.41686462E+00  0.74003433E+00</point>
      <point r="1.260">  0.68079356E+00 -0.32551093E+00  0.40524522E+00  0.73481254E+00</point>
      <point r="1.280">  0.67565161E+00 -0.32844791E+00  0.39355228E+00  0.72957179E+00</point>
      <point r="1.300">  0.67058900E+00 -0.33139438E+00  0.38178850E+00  0.72431186E+00</point>
      <point r="1.320">  0.66560486E+00 -0.33434994E+00  0.36995686E+00  0.71903256E+00</point>
      <point r="1.340">  0.66069807E+00 -0.33731391E+00  0.35806068E+00  0.71373378E+00</point>
      <point r="1.360">  0.65586734E+00 -0.34028536E+00  0.34610352E+00  0.70841544E+00</point>
      <point r="1.380">  0.65111120E+00 -0.34326314E+00  0.33408922E+00  0.70307755E+00</point>
      <point r="1.400">  0.64642801E+00 -0.34624586E+00  0.32202188E+00  0.69772015E+00</point>
      <point r="1.420">  0.64181600E+00 -0.34923195E+00  0.30990578E+00  0.69234336E+00</point>
      <point r="1.440">  0.63727326E+00 -0.35221965E+00  0.29774545E+00  0.68694733E+00</point>
      <point r="1.460">  0.63279777E+00 -0.35520704E+00  0.28554560E+00  0.68153228E+00</point>
      <point r="1.480">  0.62838740E+00 -0.35819205E+00  0.27331111E+00  0.67609848E+00</point>
      <point r="1.500">  0.62403995E+00 -0.36117247E+00  0.26104702E+00  0.67064625E+00</point>
      <point r="1.520">  0.61975314E+00 -0.36414598E+00  0.24875850E+00  0.66517596E+00</point>
      <point r="1.540">  0.61552462E+00 -0.36711013E+00  0.23645086E+00  0.65968804E+00</point>
      <point r="1.560">  0.61135198E+00 -0.37006242E+00  0.22412952E+00  0.65418295E+00</point>
      <point r="1.580">  0.60723280E+00 -0.37300024E+00  0.21179999E+00  0.64866121E+00</point>
      <point r="1.600">  0.60316460E+00 -0.37592092E+00  0.19946786E+00  0.64312337E+00</point>
      <point r="1.620">  0.59914490E+00 -0.37882175E+00  0.18713879E+00  0.63757003E+00</point>
      <point r="1.640">  0.59517119E+00 -0.38169995E+00  0.17481851E+00  0.63200183E+00</point>
      <point r="1.660">  0.59124098E+00 -0.38455274E+00  0.16251276E+00  0.62641945E+00</point>
      <point r="1.680">  0.58735176E+00 -0.38737728E+00  0.15022734E+00  0.62082360E+00</point>
      <point r="1.700">  0.58350106E+00 -0.39017076E+00  0.13796805E+00  0.61521503E+00</point>
      <point r="1.720">  0.57968640E+00 -0.39293034E+00  0.12574070E+00  0.60959452E+00</point>
      <point r="1.740">  0.57590535E+00 -0.39565319E+00  0.11355110E+00  0.60396289E+00</point>
      <point r="1.760">  0.57215550E+00 -0.39833650E+00  0.10140503E+00  0.59832097E+00</point>
      <point r="1.780">  0.56843448E+00 -0.40097746E+00  0.89308263E-01  0.59266964E+00</point>
      <point r="1.800">  0.56473996E+00 -0.40357333E+00  0.77266529E-01  0.58700978E+00</point>
      <point r="1.820">  0.56106966E+00 -0.40612137E+00  0.65285513E-01  0.58134232E+00</point>
      <point r="1.840">  0.55742134E+00 -0.40861888E+00  0.53370850E-01  0.57566819E+00</point>
      <point r="1.860">  0.55379284E+00 -0.41106324E+00  0.41528110E-01  0.56998835E+00</point>
      <point r="1.880">  0.55018203E+00 -0.41345186E+00  0.29762799E-01  0.56430378E+00</point>
      <point r="1.900">  0.54658687E+00 -0.41578220E+00  0.18080344E-01  0.55861548E+00</point>
      <point r="1.920">  0.54300534E+00 -0.41805180E+00  0.64860877E-02  0.55292444E+00</point>
      <point r="1.940">  0.53943554E+00 -0.42025826E+00 -0.50147147E-02  0.54723169E+00</point>
      <point r="1.960">  0.53587559E+00 -0.42239926E+00 -0.16416906E-01  0.54153826E+00</point>
      <point r="1.980">  0.53232373E+00 -0.42447255E+00 -0.27715431E-01  0.53584519E+00</point>
      <point r="2.000">  0.52877822E+00 -0.42647594E+00 -0.38905342E-01  0.53015355E+00</point>
      <point r="2.020">  0.52523742E+00 -0.42840735E+00 -0.49981805E-01  0.52446437E+00</point>
      <point r="2.040">  0.52169978E+00 -0.43026478E+00 -0.60940103E-01  0.51877873E+00</point>
      <point r="2.060">  0.51816379E+00 -0.43204629E+00 -0.71775642E-01  0.51309769E+00</point>
      <point r="2.080">  0.51462804E+00 -0.43375005E+00 -0.82483956E-01  0.50742232E+00</point>
      <point r="2.100">  0.51109119E+00 -0.43537433E+00 -0.93060706E-01  0.50175370E+00</point>
      <point r="2.120">  0.50755197E+00 -0.43691748E+00 -0.10350169E+00  0.49609288E+00</point>
      <point r="2.140">  0.50400920E+00 -0.43837793E+00 -0.11380284E+00  0.49044095E+00</point>
      <point r="2.160">  0.50046176E+00 -0.43975423E+00 -0.12396024E+00  0.48479896E+00</point>
      <point r="2.180">  0.49690861E+00 -0.44104501E+00 -0.13397009E+00  0.47916798E+00</point>
      <point r="2.200">  0.49334880E+00 -0.44224899E+00 -0.14382877E+00  0.47354907E+00</point>
      <point r="2.220">  0.48978144E+00 -0.44336502E+00 -0.15353278E+00  0.46794328E+00</point>
      <point r="2.240">  0.48620571E+00 -0.44439201E+00 -0.16307878E+00  0.46235165E+00</point>
      <point r="2.260">  0.48262087E+00 -0.44532898E+00 -0.17246358E+00  0.45677522E+00</point>
      <point r="2.280">  0.47902627E+00 -0.44617506E+00 -0.18168415E+00  0.45121502E+00</point>
      <point r="2.300">  0.47542129E+00 -0.44692945E+00 -0.19073759E+00  0.44567207E+00</point>
      <point r="2.320">  0.47180543E+00 -0.44759147E+00 -0.19962118E+00  0.44014736E+00</point>
      <point r="2.340">  0.46817820E+00 -0.44816053E+00 -0.20833234E+00  0.43464190E+00</point>
      <point r="2.360">  0.46453923E+00 -0.44863614E+00 -0.21686864E+00  0.42915667E+00</point>
      <point r="2.380">  0.46088819E+00 -0.44901787E+00 -0.22522783E+00  0.42369264E+00</point>
      <point r="2.400">  0.45722481E+00 -0.44930544E+00 -0.23340778E+00  0.41825076E+00</point>
      <point r="2.420">  0.45354891E+00 -0.44949862E+00 -0.24140654E+00  0.41283198E+00</point>
      <point r="2.440">  0.44986033E+00 -0.44959727E+00 -0.24922231E+00  0.40743721E+00</point>
      <point r="2.460">  0.44615900E+00 -0.44960138E+00 -0.25685343E+00  0.40206737E+00</point>
      <point r="2.480">  0.44244490E+00 -0.44951097E+00 -0.26429841E+00  0.39672334E+00</point>
      <point r="2.500">  0.43871805E+00 -0.44932620E+00 -0.27155590E+00  0.39140602E+00</point>
      <point r="2.520">  0.43497856E+00 -0.44904728E+00 -0.27862470E+00  0.38611624E+00</point>
      <point r="2.540">  0.43122655E+00 -0.44867452E+00 -0.28550377E+00  0.38085485E+00</point>
      <point r="2.560">  0.42746222E+00 -0.44820829E+00 -0.29219222E+00  0.37562267E+00</point>
      <point r="2.580">  0.42368580E+00 -0.44764907E+00 -0.29868927E+00  0.37042050E+00</point>
      <point r="2.600">  0.41989758E+00 -0.44699740E+00 -0.30499433E+00  0.36524913E+00</point>
      <point r="2.620">  0.41609788E+00 -0.44625388E+00 -0.31110693E+00  0.36010932E+00</point>
      <point r="2.640">  0.41228708E+00 -0.44541922E+00 -0.31702673E+00  0.35500180E+00</point>
      <point r="2.660">  0.40846558E+00 -0.44449416E+00 -0.32275356E+00  0.34992730E+00</point>
      <point r="2.680">  0.40463383E+00 -0.44347954E+00 -0.32828735E+00  0.34488653E+00</point>
      <point r="2.700">  0.40079233E+00 -0.44237625E+00 -0.33362818E+00  0.33988016E+00</point>
      <point r="2.720">  0.39694158E+00 -0.44118525E+00 -0.33877627E+00  0.33490884E+00</point>
      <point r="2.740">  0.39308214E+00 -0.43990756E+00 -0.34373194E+00  0.32997323E+00</point>
      <point r="2.760">  0.38921459E+00 -0.43854426E+00 -0.34849567E+00  0.32507393E+00</point>
      <point r="2.780">  0.38533955E+00 -0.43709649E+00 -0.35306804E+00  0.32021155E+00</point>
      <point r="2.800">  0.38145765E+00 -0.43556544E+00 -0.35744975E+00  0.31538665E+00</point>
      <point r="2.820">  0.37756957E+00 -0.43395234E+00 -0.36164162E+00  0.31059978E+00</point>
      <point r="2.840">  0.37367598E+00 -0.43225850E+00 -0.36564460E+00  0.30585148E+00</point>
      <point r="2.860">  0.36977759E+00 -0.43048526E+00 -0.36945972E+00  0.30114225E+00</point>
      <point r="2.880">  0.36587515E+00 -0.42863400E+00 -0.37308814E+00  0.29647259E+00</point>
      <point r="2.900">  0.36196938E+00 -0.42670616E+00 -0.37653111E+00  0.29184295E+00</point>
      <point r="2.920">  0.35806107E+00 -0.42470319E+00 -0.37979001E+00  0.28725378E+00</point>
      <point r="2.940">  0.35415097E+00 -0.42262661E+00 -0.38286627E+00  0.28270550E+00</point>
      <point r="2.960">  0.35023990E+00 -0.42047797E+00 -0.38576147E+00  0.27819851E+00</point>
      <point r="2.980">  0.34632864E+00 -0.41825883E+00 -0.38847723E+00  0.27373319E+00</point>
      <point r="3.000">  0.34241801E+00 -0.41597081E+00 -0.39101530E+00  0.26930990E+00</point>
      <point r="3.020">  0.33850884E+00 -0.41361554E+00 -0.39337748E+00  0.26492896E+00</point>
      <point r="3.040">  0.33460195E+00 -0.41119469E+00 -0.39556568E+00  0.26059071E+00</point>
      <point r="3.060">  0.33069817E+00 -0.40870994E+00 -0.39758188E+00  0.25629542E+00</point>
      <point r="3.080">  0.32679835E+00 -0.40616300E+00 -0.39942813E+00  0.25204338E+00</point>
      <point r="3.100">  0.32290332E+00 -0.40355561E+00 -0.40110654E+00  0.24783483E+00</point>
      <point r="3.120">  0.31901392E+00 -0.40088950E+00 -0.40261932E+00  0.24367002E+00</point>
      <point r="3.140">  0.31513101E+00 -0.39816644E+00 -0.40396871E+00  0.23954914E+00</point>
      <point r="3.160">  0.31125541E+00 -0.39538821E+00 -0.40515704E+00  0.23547240E+00</point>
      <point r="3.180">  0.30738797E+00 -0.39255659E+00 -0.40618669E+00  0.23143996E+00</point>
      <point r="3.200">  0.30352952E+00 -0.38967338E+00 -0.40706008E+00  0.22745199E+00</point>
      <point r="3.220">  0.29968090E+00 -0.38674038E+00 -0.40777970E+00  0.22350860E+00</point>
      <point r="3.240">  0.29584293E+00 -0.38375941E+00 -0.40834808E+00  0.21960993E+00</point>
      <point r="3.260">  0.29201642E+00 -0.38073228E+00 -0.40876781E+00  0.21575606E+00</point>
      <point r="3.280">  0.28820220E+00 -0.37766080E+00 -0.40904149E+00  0.21194708E+00</point>
      <point r="3.300">  0.28440106E+00 -0.37454680E+00 -0.40917179E+00  0.20818305E+00</point>
      <point r="3.320">  0.28061380E+00 -0.37139209E+00 -0.40916141E+00  0.20446401E+00</point>
      <point r="3.340">  0.27684120E+00 -0.36819847E+00 -0.40901308E+00  0.20078998E+00</point>
      <point r="3.360">  0.27308404E+00 -0.36496777E+00 -0.40872955E+00  0.19716098E+00</point>
      <point r="3.380">  0.26934307E+00 -0.36170177E+00 -0.40831361E+00  0.19357699E+00</point>
      <point r="3.400">  0.26561905E+00 -0.35840227E+00 -0.40776808E+00  0.19003800E+00</point>
      <point r="3.420">  0.26191271E+00 -0.35507107E+00 -0.40709578E+00  0.18654396E+00</point>
      <point r="3.440">  0.25822477E+00 -0.35170992E+00 -0.40629957E+00  0.18309482E+00</point>
      <point r="3.460">  0.25455595E+00 -0.34832060E+00 -0.40538232E+00  0.17969051E+00</point>
      <point r="3.480">  0.25090693E+00 -0.34490485E+00 -0.40434690E+00  0.17633095E+00</point>
      <point r="3.500">  0.24727839E+00 -0.34146441E+00 -0.40319621E+00  0.17301602E+00</point>
      <point r="3.520">  0.24367100E+00 -0.33800100E+00 -0.40193314E+00  0.16974562E+00</point>
      <point r="3.540">  0.24008541E+00 -0.33451631E+00 -0.40056062E+00  0.16651962E+00</point>
      <point r="3.560">  0.23652223E+00 -0.33101204E+00 -0.39908153E+00  0.16333789E+00</point>
      <point r="3.580">  0.23298209E+00 -0.32748985E+00 -0.39749880E+00  0.16020025E+00</point>
      <point r="3.600">  0.22946559E+00 -0.32395138E+00 -0.39581534E+00  0.15710655E+00</point>
      <point r="3.620">  0.22597329E+00 -0.32039826E+00 -0.39403404E+00  0.15405661E+00</point>
      <point r="3.640">  0.22250576E+00 -0.31683209E+00 -0.39215782E+00  0.15105024E+00</point>
      <point r="3.660">  0.21906354E+00 -0.31325445E+00 -0.39018955E+00  0.14808722E+00</point>
      <point r="3.680">  0.21564716E+00 -0.30966690E+00 -0.38813213E+00  0.14516736E+00</point>
      <point r="3.700">  0.21225712E+00 -0.30607098E+00 -0.38598842E+00  0.14229041E+00</point>
      <point r="3.720">  0.20889391E+00 -0.30246818E+00 -0.38376129E+00  0.13945616E+00</point>
      <point r="3.740">  0.20555800E+00 -0.29886000E+00 -0.38145356E+00  0.13666434E+00</point>
      <point r="3.760">  0.20224983E+00 -0.29524788E+00 -0.37906808E+00  0.13391470E+00</point>
      <point r="3.780">  0.19896984E+00 -0.29163327E+00 -0.37660763E+00  0.13120698E+00</point>
      <point r="3.800">  0.19571844E+00 -0.28801755E+00 -0.37407501E+00  0.12854091E+00</point>
      <point r="3.820">  0.19249603E+00 -0.28440212E+00 -0.37147298E+00  0.12591619E+00</point>
      <point r="3.840">  0.18930297E+00 -0.28078830E+00 -0.36880427E+00  0.12333255E+00</point>
      <point r="3.860">  0.18613963E+00 -0.27717743E+00 -0.36607160E+00  0.12078967E+00</point>
      <point r="3.880">  0.18300634E+00 -0.27357079E+00 -0.36327766E+00  0.11828725E+00</point>
      <point r="3.900">  0.17990342E+00 -0.26996965E+00 -0.36042510E+00  0.11582497E+00</point>
      <point r="3.920">  0.17683117E+00 -0.26637523E+00 -0.35751654E+00  0.11340252E+00</point>
      <point r="3.940">  0.17378987E+00 -0.26278873E+00 -0.35455460E+00  0.11101957E+00</point>
      <point r="3.960">  0.17077978E+00 -0.25921133E+00 -0.35154184E+00  0.10867577E+00</point>
      <point r="3.980">  0.16780115E+00 -0.25564417E+00 -0.34848079E+00  0.10637079E+00</point>
      <point r="4.000">  0.16485421E+00 -0.25208837E+00 -0.34537394E+00  0.10410428E+00</point>
      <point r="4.020">  0.16193916E+00 -0.24854499E+00 -0.34222378E+00  0.10187589E+00</point>
      <point r="4.040">  0.15905620E+00 -0.24501510E+00 -0.33903273E+00  0.99685257E-01</point>
      <point r="4.060">  0.15620551E+00 -0.24149971E+00 -0.33580318E+00  0.97532024E-01</point>
      <point r="4.080">  0.15338723E+00 -0.23799982E+00 -0.33253748E+00  0.95415822E-01</point>
      <point r="4.100">  0.15060151E+00 -0.23451638E+00 -0.32923797E+00  0.93336279E-01</point>
      <point r="4.120">  0.14784848E+00 -0.23105033E+00 -0.32590692E+00  0.91293020E-01</point>
      <point r="4.140">  0.14512823E+00 -0.22760256E+00 -0.32254657E+00  0.89285667E-01</point>
      <point r="4.160">  0.14244087E+00 -0.22417395E+00 -0.31915912E+00  0.87313837E-01</point>
      <point r="4.180">  0.13978646E+00 -0.22076533E+00 -0.31574674E+00  0.85377144E-01</point>
      <point r="4.200">  0.13716507E+00 -0.21737752E+00 -0.31231155E+00  0.83475201E-01</point>
      <point r="4.220">  0.13457675E+00 -0.21401129E+00 -0.30885562E+00  0.81607616E-01</point>
      <point r="4.240">  0.13202151E+00 -0.21066739E+00 -0.30538100E+00  0.79773996E-01</point>
      <point r="4.260">  0.12949938E+00 -0.20734656E+00 -0.30188968E+00  0.77973946E-01</point>
      <point r="4.280">  0.12701035E+00 -0.20404948E+00 -0.29838362E+00  0.76207069E-01</point>
      <point r="4.300">  0.12455442E+00 -0.20077681E+00 -0.29486473E+00  0.74472965E-01</point>
      <point r="4.320">  0.12213155E+00 -0.19752919E+00 -0.29133488E+00  0.72771236E-01</point>
      <point r="4.340">  0.11974171E+00 -0.19430723E+00 -0.28779589E+00  0.71101479E-01</point>
      <point r="4.360">  0.11738484E+00 -0.19111152E+00 -0.28424956E+00  0.69463292E-01</point>
      <point r="4.380">  0.11506087E+00 -0.18794259E+00 -0.28069762E+00  0.67856272E-01</point>
      <point r="4.400">  0.11276974E+00 -0.18480098E+00 -0.27714178E+00  0.66280017E-01</point>
      <point r="4.420">  0.11051133E+00 -0.18168719E+00 -0.27358368E+00  0.64734121E-01</point>
      <point r="4.440">  0.10828556E+00 -0.17860169E+00 -0.27002495E+00  0.63218182E-01</point>
      <point r="4.460">  0.10609232E+00 -0.17554493E+00 -0.26646715E+00  0.61731796E-01</point>
      <point r="4.480">  0.10393146E+00 -0.17251732E+00 -0.26291182E+00  0.60274560E-01</point>
      <point r="4.500">  0.10180287E+00 -0.16951927E+00 -0.25936044E+00  0.58846069E-01</point>
      <point r="4.520">  0.99706387E-01 -0.16655114E+00 -0.25581445E+00  0.57445923E-01</point>
      <point r="4.540">  0.97641865E-01 -0.16361329E+00 -0.25227526E+00  0.56073720E-01</point>
      <point r="4.560">  0.95609137E-01 -0.16070603E+00 -0.24874422E+00  0.54729059E-01</point>
      <point r="4.580">  0.93608026E-01 -0.15782966E+00 -0.24522266E+00  0.53411541E-01</point>
      <point r="4.600">  0.91638351E-01 -0.15498446E+00 -0.24171185E+00  0.52120768E-01</point>
      <point r="4.620">  0.89699918E-01 -0.15217069E+00 -0.23821303E+00  0.50856343E-01</point>
      <point r="4.640">  0.87792527E-01 -0.14938857E+00 -0.23472740E+00  0.49617871E-01</point>
      <point r="4.660">  0.85915970E-01 -0.14663832E+00 -0.23125611E+00  0.48404959E-01</point>
      <point r="4.680">  0.84070031E-01 -0.14392013E+00 -0.22780028E+00  0.47217217E-01</point>
      <point r="4.700">  0.82254485E-01 -0.14123416E+00 -0.22436099E+00  0.46054253E-01</point>
      <point r="4.720">  0.80469102E-01 -0.13858056E+00 -0.22093926E+00  0.44915682E-01</point>
      <point r="4.740">  0.78713646E-01 -0.13595946E+00 -0.21753612E+00  0.43801117E-01</point>
      <point r="4.760">  0.76987873E-01 -0.13337098E+00 -0.21415250E+00  0.42710177E-01</point>
      <point r="4.780">  0.75291533E-01 -0.13081519E+00 -0.21078935E+00  0.41642481E-01</point>
      <point r="4.800">  0.73624371E-01 -0.12829217E+00 -0.20744754E+00  0.40597651E-01</point>
      <point r="4.820">  0.71986125E-01 -0.12580198E+00 -0.20412792E+00  0.39575313E-01</point>
      <point r="4.840">  0.70376531E-01 -0.12334465E+00 -0.20083132E+00  0.38575094E-01</point>
      <point r="4.860">  0.68795317E-01 -0.12092020E+00 -0.19755850E+00  0.37596625E-01</point>
      <point r="4.880">  0.67242207E-01 -0.11852864E+00 -0.19431021E+00  0.36639539E-01</point>
      <point r="4.900">  0.65716922E-01 -0.11616995E+00 -0.19108717E+00  0.35703474E-01</point>
      <point r="4.920">  0.64219179E-01 -0.11384410E+00 -0.18789004E+00  0.34788068E-01</point>
      <point r="4.940">  0.62748689E-01 -0.11155105E+00 -0.18471947E+00  0.33892964E-01</point>
      <point r="4.960">  0.61305161E-01 -0.10929075E+00 -0.18157607E+00  0.33017809E-01</point>
      <point r="4.980">  0.59888302E-01 -0.10706311E+00 -0.17846041E+00  0.32162252E-01</point>
      <point r="5.000">  0.58497813E-01 -0.10486807E+00 -0.17537304E+00  0.31325946E-01</point>
      <point r="5.020">  0.57133395E-01 -0.10270550E+00 -0.17231448E+00  0.30508546E-01</point>
      <point r="5.040">  0.55794746E-01 -0.10057532E+00 -0.16928521E+00  0.29709713E-01</point>
      <point r="5.060">  0.54481561E-01 -0.98477383E-01 -0.16628568E+00  0.28929110E-01</point>
      <point r="5.080">  0.53193532E-01 -0.96411563E-01 -0.16331632E+00  0.28166403E-01</point>
      <point r="5.100">  0.51930351E-01 -0.94377712E-01 -0.16037753E+00  0.27421262E-01</point>
      <point r="5.120">  0.50691709E-01 -0.92375671E-01 -0.15746967E+00  0.26693363E-01</point>
      <point r="5.140">  0.49477293E-01 -0.90405272E-01 -0.15459309E+00  0.25982382E-01</point>
      <point r="5.160">  0.48286791E-01 -0.88466336E-01 -0.15174810E+00  0.25288001E-01</point>
      <point r="5.180">  0.47119890E-01 -0.86558674E-01 -0.14893500E+00  0.24609905E-01</point>
      <point r="5.200">  0.45976274E-01 -0.84682088E-01 -0.14615404E+00  0.23947783E-01</point>
      <point r="5.220">  0.44855630E-01 -0.82836372E-01 -0.14340546E+00  0.23301327E-01</point>
      <point r="5.240">  0.43757641E-01 -0.81021310E-01 -0.14068949E+00  0.22670235E-01</point>
      <point r="5.260">  0.42681992E-01 -0.79236679E-01 -0.13800630E+00  0.22054207E-01</point>
      <point r="5.280">  0.41628369E-01 -0.77482247E-01 -0.13535608E+00  0.21452947E-01</point>
      <point r="5.300">  0.40596455E-01 -0.75757776E-01 -0.13273896E+00  0.20866164E-01</point>
      <point r="5.320">  0.39585935E-01 -0.74063021E-01 -0.13015506E+00  0.20293568E-01</point>
      <point r="5.340">  0.38596495E-01 -0.72397728E-01 -0.12760450E+00  0.19734877E-01</point>
      <point r="5.360">  0.37627821E-01 -0.70761639E-01 -0.12508736E+00  0.19189810E-01</point>
      <point r="5.380">  0.36679600E-01 -0.69154490E-01 -0.12260369E+00  0.18658090E-01</point>
      <point r="5.400">  0.35751520E-01 -0.67576011E-01 -0.12015355E+00  0.18139446E-01</point>
      <point r="5.420">  0.34843269E-01 -0.66025925E-01 -0.11773695E+00  0.17633609E-01</point>
      <point r="5.440">  0.33954538E-01 -0.64503953E-01 -0.11535390E+00  0.17140314E-01</point>
      <point r="5.460">  0.33085018E-01 -0.63009809E-01 -0.11300439E+00  0.16659300E-01</point>
      <point r="5.480">  0.32234402E-01 -0.61543205E-01 -0.11068839E+00  0.16190312E-01</point>
      <point r="5.500">  0.31402384E-01 -0.60103846E-01 -0.10840587E+00  0.15733095E-01</point>
      <point r="5.520">  0.30588661E-01 -0.58691436E-01 -0.10615675E+00  0.15287400E-01</point>
      <point r="5.540">  0.29792929E-01 -0.57305673E-01 -0.10394096E+00  0.14852983E-01</point>
      <point r="5.560">  0.29014890E-01 -0.55946255E-01 -0.10175841E+00  0.14429603E-01</point>
      <point r="5.580">  0.28254244E-01 -0.54612875E-01 -0.99608991E-01  0.14017020E-01</point>
      <point r="5.600">  0.27510695E-01 -0.53305223E-01 -0.97492593E-01  0.13615002E-01</point>
      <point r="5.620">  0.26783950E-01 -0.52022989E-01 -0.95409078E-01  0.13223318E-01</point>
      <point r="5.640">  0.26073717E-01 -0.50765857E-01 -0.93358302E-01  0.12841743E-01</point>
      <point r="5.660">  0.25379705E-01 -0.49533514E-01 -0.91340109E-01  0.12470053E-01</point>
      <point r="5.680">  0.24701629E-01 -0.48325641E-01 -0.89354328E-01  0.12108030E-01</point>
      <point r="5.700">  0.24039203E-01 -0.47141921E-01 -0.87400779E-01  0.11755458E-01</point>
      <point r="5.720">  0.23392146E-01 -0.45982034E-01 -0.85479273E-01  0.11412127E-01</point>
      <point r="5.740">  0.22760178E-01 -0.44845659E-01 -0.83589607E-01  0.11077827E-01</point>
      <point r="5.760">  0.22143022E-01 -0.43732476E-01 -0.81731571E-01  0.10752356E-01</point>
      <point r="5.780">  0.21540404E-01 -0.42642162E-01 -0.79904943E-01  0.10435511E-01</point>
      <point r="5.800">  0.20952054E-01 -0.41574395E-01 -0.78109495E-01  0.10127097E-01</point>
      <point r="5.820">  0.20377702E-01 -0.40528854E-01 -0.76344989E-01  0.98269192E-02</point>
      <point r="5.840">  0.19817083E-01 -0.39505217E-01 -0.74611178E-01  0.95347875E-02</point>
      <point r="5.860">  0.19269936E-01 -0.38503161E-01 -0.72907808E-01  0.92505153E-02</point>
      <point r="5.880">  0.18735999E-01 -0.37522365E-01 -0.71234619E-01  0.89739193E-02</point>
      <point r="5.900">  0.18215016E-01 -0.36562508E-01 -0.69591342E-01  0.87048196E-02</point>
      <point r="5.920">  0.17706734E-01 -0.35623271E-01 -0.67977704E-01  0.84430395E-02</point>
      <point r="5.940">  0.17210903E-01 -0.34704335E-01 -0.66393423E-01  0.81884060E-02</point>
      <point r="5.960">  0.16727274E-01 -0.33805380E-01 -0.64838214E-01  0.79407489E-02</point>
      <point r="5.980">  0.16255604E-01 -0.32926091E-01 -0.63311784E-01  0.76999015E-02</point>
      <point r="6.000">  0.15795651E-01 -0.32066153E-01 -0.61813838E-01  0.74657004E-02</point>
      <point r="6.020">  0.15347176E-01 -0.31225250E-01 -0.60344073E-01  0.72379850E-02</point>
      <point r="6.040">  0.14909947E-01 -0.30403070E-01 -0.58902185E-01  0.70165982E-02</point>
      <point r="6.060">  0.14483729E-01 -0.29599304E-01 -0.57487864E-01  0.68013859E-02</point>
      <point r="6.080">  0.14068296E-01 -0.28813641E-01 -0.56100796E-01  0.65921969E-02</point>
      <point r="6.100">  0.13663421E-01 -0.28045776E-01 -0.54740665E-01  0.63888831E-02</point>
      <point r="6.120">  0.13268882E-01 -0.27295402E-01 -0.53407150E-01  0.61912995E-02</point>
      <point r="6.140">  0.12884461E-01 -0.26562218E-01 -0.52099931E-01  0.59993039E-02</point>
      <point r="6.160">  0.12509942E-01 -0.25845923E-01 -0.50818681E-01  0.58127569E-02</point>
      <point r="6.180">  0.12145111E-01 -0.25146218E-01 -0.49563073E-01  0.56315223E-02</point>
      <point r="6.200">  0.11789761E-01 -0.24462809E-01 -0.48332779E-01  0.54554663E-02</point>
      <point r="6.220">  0.11443685E-01 -0.23795401E-01 -0.47127467E-01  0.52844582E-02</point>
      <point r="6.240">  0.11106679E-01 -0.23143704E-01 -0.45946805E-01  0.51183699E-02</point>
      <point r="6.260">  0.10778545E-01 -0.22507430E-01 -0.44790460E-01  0.49570760E-02</point>
      <point r="6.280">  0.10459085E-01 -0.21886294E-01 -0.43658098E-01  0.48004538E-02</point>
      <point r="6.300">  0.10148106E-01 -0.21280013E-01 -0.42549383E-01  0.46483832E-02</point>
      <point r="6.320">  0.98454183E-02 -0.20688307E-01 -0.41463981E-01  0.45007467E-02</point>
      <point r="6.340">  0.95508344E-02 -0.20110900E-01 -0.40401554E-01  0.43574293E-02</point>
      <point r="6.360">  0.92641705E-02 -0.19547518E-01 -0.39361768E-01  0.42183185E-02</point>
      <point r="6.380">  0.89852456E-02 -0.18997890E-01 -0.38344286E-01  0.40833044E-02</point>
      <point r="6.400">  0.87138820E-02 -0.18461747E-01 -0.37348773E-01  0.39522794E-02</point>
      <point r="6.420">  0.84499050E-02 -0.17938826E-01 -0.36374894E-01  0.38251382E-02</point>
      <point r="6.440">  0.81931430E-02 -0.17428865E-01 -0.35422316E-01  0.37017780E-02</point>
      <point r="6.460">  0.79434273E-02 -0.16931604E-01 -0.34490704E-01  0.35820983E-02</point>
      <point r="6.480">  0.77005924E-02 -0.16446790E-01 -0.33579726E-01  0.34660007E-02</point>
      <point r="6.500">  0.74644756E-02 -0.15974169E-01 -0.32689051E-01  0.33533893E-02</point>
      <point r="6.520">  0.72349172E-02 -0.15513492E-01 -0.31818350E-01  0.32441701E-02</point>
      <point r="6.540">  0.70117606E-02 -0.15064515E-01 -0.30967294E-01  0.31382516E-02</point>
      <point r="6.560">  0.67948519E-02 -0.14626995E-01 -0.30135556E-01  0.30355440E-02</point>
      <point r="6.580">  0.65840402E-02 -0.14200692E-01 -0.29322811E-01  0.29359600E-02</point>
      <point r="6.600">  0.63791773E-02 -0.13785371E-01 -0.28528736E-01  0.28394139E-02</point>
      <point r="6.620">  0.61801181E-02 -0.13380799E-01 -0.27753010E-01  0.27458225E-02</point>
      <point r="6.640">  0.59867201E-02 -0.12986747E-01 -0.26995315E-01  0.26551041E-02</point>
      <point r="6.660">  0.57988435E-02 -0.12602989E-01 -0.26255333E-01  0.25671793E-02</point>
      <point r="6.680">  0.56163513E-02 -0.12229304E-01 -0.25532750E-01  0.24819704E-02</point>
      <point r="6.700">  0.54391094E-02 -0.11865470E-01 -0.24827253E-01  0.23994015E-02</point>
      <point r="6.720">  0.52669861E-02 -0.11511274E-01 -0.24138535E-01  0.23193988E-02</point>
      <point r="6.740">  0.50998525E-02 -0.11166501E-01 -0.23466287E-01  0.22418901E-02</point>
      <point r="6.760">  0.49375822E-02 -0.10830943E-01 -0.22810205E-01  0.21668048E-02</point>
      <point r="6.780">  0.47800514E-02 -0.10504394E-01 -0.22169989E-01  0.20940745E-02</point>
      <point r="6.800">  0.46271391E-02 -0.10186651E-01 -0.21545339E-01  0.20236320E-02</point>
      <point r="6.820">  0.44787263E-02 -0.98775154E-02 -0.20935960E-01  0.19554120E-02</point>
      <point r="6.840">  0.43346970E-02 -0.95767903E-02 -0.20341560E-01  0.18893507E-02</point>
      <point r="6.860">  0.41949373E-02 -0.92842833E-02 -0.19761848E-01  0.18253862E-02</point>
      <point r="6.880">  0.40593359E-02 -0.89998046E-02 -0.19196538E-01  0.17634577E-02</point>
      <point r="6.900">  0.39277837E-02 -0.87231680E-02 -0.18645348E-01  0.17035063E-02</point>
      <point r="6.920">  0.38001741E-02 -0.84541904E-02 -0.18107996E-01  0.16454744E-02</point>
      <point r="6.940">  0.36764028E-02 -0.81926919E-02 -0.17584206E-01  0.15893060E-02</point>
      <point r="6.960">  0.35563677E-02 -0.79384956E-02 -0.17073704E-01  0.15349463E-02</point>
      <point r="6.980">  0.34399690E-02 -0.76914281E-02 -0.16576221E-01  0.14823422E-02</point>
      <point r="7.000">  0.33271090E-02 -0.74513188E-02 -0.16091489E-01  0.14314418E-02</point>
      <point r="7.020">  0.32176923E-02 -0.72180003E-02 -0.15619244E-01  0.13821945E-02</point>
      <point r="7.040">  0.31116257E-02 -0.69913084E-02 -0.15159226E-01  0.13345513E-02</point>
      <point r="7.060">  0.30088178E-02 -0.67710819E-02 -0.14711180E-01  0.12884641E-02</point>
      <point r="7.080">  0.29091797E-02 -0.65571624E-02 -0.14274851E-01  0.12438864E-02</point>
      <point r="7.100">  0.28126241E-02 -0.63493948E-02 -0.13849990E-01  0.12007727E-02</point>
      <point r="7.120">  0.27190661E-02 -0.61476268E-02 -0.13436351E-01  0.11590788E-02</point>
      <point r="7.140">  0.26284224E-02 -0.59517092E-02 -0.13033690E-01  0.11187618E-02</point>
      <point r="7.160">  0.25406120E-02 -0.57614954E-02 -0.12641769E-01  0.10797796E-02</point>
      <point r="7.180">  0.24555555E-02 -0.55768421E-02 -0.12260352E-01  0.10420917E-02</point>
      <point r="7.200">  0.23731755E-02 -0.53976084E-02 -0.11889206E-01  0.10056583E-02</point>
      <point r="7.220">  0.22933965E-02 -0.52236566E-02 -0.11528104E-01  0.97044087E-03</point>
      <point r="7.240">  0.22161447E-02 -0.50548515E-02 -0.11176819E-01  0.93640194E-03</point>
      <point r="7.260">  0.21413481E-02 -0.48910608E-02 -0.10835131E-01  0.90350500E-03</point>
      <point r="7.280">  0.20689365E-02 -0.47321549E-02 -0.10502821E-01  0.87171460E-03</point>
      <point r="7.300">  0.19988413E-02 -0.45780068E-02 -0.10179674E-01  0.84099622E-03</point>
      <point r="7.320">  0.19309957E-02 -0.44284923E-02 -0.98654798E-02  0.81131634E-03</point>
      <point r="7.340">  0.18653344E-02 -0.42834896E-02 -0.95600301E-02  0.78264236E-03</point>
      <point r="7.360">  0.18017938E-02 -0.41428797E-02 -0.92631211E-02  0.75494259E-03</point>
      <point r="7.380">  0.17403120E-02 -0.40065460E-02 -0.89745518E-02  0.72818624E-03</point>
      <point r="7.400">  0.16808284E-02 -0.38743744E-02 -0.86941253E-02  0.70234339E-03</point>
      <point r="7.420">  0.16232840E-02 -0.37462534E-02 -0.84216474E-02  0.67738497E-03</point>
      <point r="7.440">  0.15676215E-02 -0.36220739E-02 -0.81569279E-02  0.65328273E-03</point>
      <point r="7.460">  0.15137849E-02 -0.35017290E-02 -0.78997798E-02  0.63000923E-03</point>
      <point r="7.480">  0.14617195E-02 -0.33851145E-02 -0.76500192E-02  0.60753784E-03</point>
      <point r="7.500">  0.14113723E-02 -0.32721282E-02 -0.74074658E-02  0.58584268E-03</point>
      <point r="7.520">  0.13626914E-02 -0.31626705E-02 -0.71719427E-02  0.56489861E-03</point>
      <point r="7.540">  0.13156264E-02 -0.30566439E-02 -0.69432761E-02  0.54468125E-03</point>
      <point r="7.560">  0.12701282E-02 -0.29539530E-02 -0.67212954E-02  0.52516692E-03</point>
      <point r="7.580">  0.12261489E-02 -0.28545048E-02 -0.65058334E-02  0.50633264E-03</point>
      <point r="7.600">  0.11836421E-02 -0.27582085E-02 -0.62967261E-02  0.48815608E-03</point>
      <point r="7.620">  0.11425623E-02 -0.26649751E-02 -0.60938126E-02  0.47061562E-03</point>
      <point r="7.640">  0.11028655E-02 -0.25747180E-02 -0.58969351E-02  0.45369025E-03</point>
      <point r="7.660">  0.10645088E-02 -0.24873525E-02 -0.57059390E-02  0.43735959E-03</point>
      <point r="7.680">  0.10274504E-02 -0.24027959E-02 -0.55206728E-02  0.42160387E-03</point>
      <point r="7.700">  0.99164957E-03 -0.23209676E-02 -0.53409879E-02  0.40640394E-03</point>
      <point r="7.720">  0.95706689E-03 -0.22417887E-02 -0.51667389E-02  0.39174119E-03</point>
      <point r="7.740">  0.92366391E-03 -0.21651825E-02 -0.49977832E-02  0.37759761E-03</point>
      <point r="7.760">  0.89140322E-03 -0.20910740E-02 -0.48339812E-02  0.36395571E-03</point>
      <point r="7.780">  0.86024850E-03 -0.20193900E-02 -0.46751963E-02  0.35079856E-03</point>
      <point r="7.800">  0.83016439E-03 -0.19500593E-02 -0.45212947E-02  0.33810973E-03</point>
      <point r="7.820">  0.80111657E-03 -0.18830122E-02 -0.43721454E-02  0.32587331E-03</point>
      <point r="7.840">  0.77307166E-03 -0.18181811E-02 -0.42276203E-02  0.31407388E-03</point>
      <point r="7.860">  0.74599721E-03 -0.17554996E-02 -0.40875939E-02  0.30269650E-03</point>
      <point r="7.880">  0.71986172E-03 -0.16949036E-02 -0.39519435E-02  0.29172669E-03</point>
      <point r="7.900">  0.69463457E-03 -0.16363300E-02 -0.38205493E-02  0.28115043E-03</point>
      <point r="7.920">  0.67028604E-03 -0.15797178E-02 -0.36932937E-02  0.27095414E-03</point>
      <point r="7.940">  0.64678723E-03 -0.15250073E-02 -0.35700621E-02  0.26112466E-03</point>
      <point r="7.960">  0.62411013E-03 -0.14721405E-02 -0.34507423E-02  0.25164927E-03</point>
      <point r="7.980">  0.60222749E-03 -0.14210607E-02 -0.33352247E-02  0.24251564E-03</point>
      <point r="8.000">  0.58111291E-03 -0.13717131E-02 -0.32234022E-02  0.23371182E-03</point>
      <point r="8.020">  0.56074073E-03 -0.13240438E-02 -0.31151700E-02  0.22522628E-03</point>
      <point r="8.040">  0.54108605E-03 -0.12780008E-02 -0.30104259E-02  0.21704783E-03</point>
      <point r="8.060">  0.52212475E-03 -0.12335332E-02 -0.29090702E-02  0.20916565E-03</point>
      <point r="8.080">  0.50383338E-03 -0.11905915E-02 -0.28110051E-02  0.20156928E-03</point>
      <point r="8.100">  0.48618921E-03 -0.11491277E-02 -0.27161355E-02  0.19424859E-03</point>
      <point r="8.120">  0.46917022E-03 -0.11090950E-02 -0.26243685E-02  0.18719380E-03</point>
      <point r="8.140">  0.45275503E-03 -0.10704479E-02 -0.25356133E-02  0.18039543E-03</point>
      <point r="8.160">  0.43692292E-03 -0.10331419E-02 -0.24497813E-02  0.17384432E-03</point>
      <point r="8.180">  0.42165380E-03 -0.99713425E-03 -0.23667863E-02  0.16753162E-03</point>
      <point r="8.200">  0.40692820E-03 -0.96238291E-03 -0.22865438E-02  0.16144876E-03</point>
      <point r="8.220">  0.39272726E-03 -0.92884723E-03 -0.22089718E-02  0.15558748E-03</point>
      <point r="8.240">  0.37903269E-03 -0.89648765E-03 -0.21339901E-02  0.14993977E-03</point>
      <point r="8.260">  0.36582679E-03 -0.86526576E-03 -0.20615205E-02  0.14449791E-03</point>
      <point r="8.280">  0.35309240E-03 -0.83514421E-03 -0.19914868E-02  0.13925443E-03</point>
      <point r="8.300">  0.34081291E-03 -0.80608671E-03 -0.19238149E-02  0.13420211E-03</point>
      <point r="8.320">  0.32897222E-03 -0.77805804E-03 -0.18584323E-02  0.12933399E-03</point>
      <point r="8.340">  0.31755477E-03 -0.75102397E-03 -0.17952685E-02  0.12464334E-03</point>
      <point r="8.360">  0.30654548E-03 -0.72495127E-03 -0.17342549E-02  0.12012365E-03</point>
      <point r="8.380">  0.29592974E-03 -0.69980769E-03 -0.16753246E-02  0.11576866E-03</point>
      <point r="8.400">  0.28569345E-03 -0.67556192E-03 -0.16184124E-02  0.11157230E-03</point>
      <point r="8.420">  0.27582294E-03 -0.65218356E-03 -0.15634549E-02  0.10752872E-03</point>
      <point r="8.440">  0.26630498E-03 -0.62964312E-03 -0.15103905E-02  0.10363228E-03</point>
      <point r="8.460">  0.25712681E-03 -0.60791199E-03 -0.14591589E-02  0.99877539E-04</point>
      <point r="8.480">  0.24827604E-03 -0.58696243E-03 -0.14097017E-02  0.96259231E-04</point>
      <point r="8.500">  0.23974074E-03 -0.56676750E-03 -0.13619621E-02  0.92772289E-04</point>
      <point r="8.520">  0.23150935E-03 -0.54730112E-03 -0.13158847E-02  0.89411821E-04</point>
      <point r="8.540">  0.22357069E-03 -0.52853797E-03 -0.12714156E-02  0.86173109E-04</point>
      <point r="8.560">  0.21591397E-03 -0.51045352E-03 -0.12285025E-02  0.83051603E-04</point>
      <point r="8.580">  0.20852877E-03 -0.49302399E-03 -0.11870945E-02  0.80042916E-04</point>
      <point r="8.600">  0.20140500E-03 -0.47622633E-03 -0.11471422E-02  0.77142818E-04</point>
      <point r="8.620">  0.19453294E-03 -0.46003823E-03 -0.11085975E-02  0.74347233E-04</point>
      <point r="8.640">  0.18790319E-03 -0.44443806E-03 -0.10714136E-02  0.71652230E-04</point>
      <point r="8.660">  0.18150668E-03 -0.42940486E-03 -0.10355451E-02  0.69054022E-04</point>
      <point r="8.680">  0.17533464E-03 -0.41491836E-03 -0.10009480E-02  0.66548962E-04</point>
      <point r="8.700">  0.16937862E-03 -0.40095891E-03 -0.96757937E-03  0.64133532E-04</point>
      <point r="8.720">  0.16363046E-03 -0.38750750E-03 -0.93539770E-03  0.61804346E-04</point>
      <point r="8.740">  0.15808229E-03 -0.37454573E-03 -0.90436259E-03  0.59558143E-04</point>
      <point r="8.760">  0.15272652E-03 -0.36205579E-03 -0.87443485E-03  0.57391782E-04</point>
      <point r="8.780">  0.14755582E-03 -0.35002044E-03 -0.84557644E-03  0.55302238E-04</point>
      <point r="8.800">  0.14256314E-03 -0.33842304E-03 -0.81775047E-03  0.53286599E-04</point>
      <point r="8.820">  0.13774165E-03 -0.32724744E-03 -0.79092112E-03  0.51342062E-04</point>
      <point r="8.840">  0.13308480E-03 -0.31647808E-03 -0.76505367E-03  0.49465929E-04</point>
      <point r="8.860">  0.12858627E-03 -0.30609986E-03 -0.74011444E-03  0.47655604E-04</point>
      <point r="8.880">  0.12423995E-03 -0.29609824E-03 -0.71607079E-03  0.45908589E-04</point>
      <point r="8.900">  0.12003998E-03 -0.28645912E-03 -0.69289106E-03  0.44222479E-04</point>
      <point r="8.920">  0.11598071E-03 -0.27716891E-03 -0.67054456E-03  0.42594963E-04</point>
      <point r="8.940">  0.11205669E-03 -0.26821447E-03 -0.64900158E-03  0.41023817E-04</point>
      <point r="8.960">  0.10826267E-03 -0.25958308E-03 -0.62823331E-03  0.39506902E-04</point>
      <point r="8.980">  0.10459361E-03 -0.25126251E-03 -0.60821185E-03  0.38042162E-04</point>
      <point r="9.000">  0.10104465E-03 -0.24324091E-03 -0.58891019E-03  0.36627619E-04</point>
      <point r="9.020">  0.97611123E-04 -0.23550685E-03 -0.57030216E-03  0.35261373E-04</point>
      <point r="9.040">  0.94288533E-04 -0.22804931E-03 -0.55236245E-03  0.33941596E-04</point>
      <point r="9.060">  0.91072556E-04 -0.22085764E-03 -0.53506654E-03  0.32666533E-04</point>
      <point r="9.080">  0.87959036E-04 -0.21392159E-03 -0.51839071E-03  0.31434497E-04</point>
      <point r="9.100">  0.84943975E-04 -0.20723123E-03 -0.50231201E-03  0.30243866E-04</point>
      <point r="9.120">  0.82023535E-04 -0.20077704E-03 -0.48680826E-03  0.29093083E-04</point>
      <point r="9.140">  0.79194024E-04 -0.19454980E-03 -0.47185799E-03  0.27980652E-04</point>
      <point r="9.160">  0.76451896E-04 -0.18854062E-03 -0.45744044E-03  0.26905136E-04</point>
      <point r="9.180">  0.73793747E-04 -0.18274097E-03 -0.44353556E-03  0.25865155E-04</point>
      <point r="9.200">  0.71216306E-04 -0.17714260E-03 -0.43012396E-03  0.24859385E-04</point>
      <point r="9.220">  0.68716433E-04 -0.17173757E-03 -0.41718690E-03  0.23886552E-04</point>
      <point r="9.240">  0.66291116E-04 -0.16651823E-03 -0.40470629E-03  0.22945435E-04</point>
      <point r="9.260">  0.63937463E-04 -0.16147722E-03 -0.39266466E-03  0.22034863E-04</point>
      <point r="9.280">  0.61652700E-04 -0.15660745E-03 -0.38104511E-03  0.21153709E-04</point>
      <point r="9.300">  0.59434165E-04 -0.15190210E-03 -0.36983137E-03  0.20300894E-04</point>
      <point r="9.320">  0.57279308E-04 -0.14735461E-03 -0.35900772E-03  0.19475382E-04</point>
      <point r="9.340">  0.55185682E-04 -0.14295865E-03 -0.34855896E-03  0.18676179E-04</point>
      <point r="9.360">  0.53150944E-04 -0.13870816E-03 -0.33847048E-03  0.17902329E-04</point>
      <point r="9.380">  0.51172846E-04 -0.13459729E-03 -0.32872815E-03  0.17152919E-04</point>
      <point r="9.400">  0.49249236E-04 -0.13062042E-03 -0.31931835E-03  0.16427070E-04</point>
      <point r="9.420">  0.47378054E-04 -0.12677217E-03 -0.31022798E-03  0.15723940E-04</point>
      <point r="9.440">  0.45557326E-04 -0.12304733E-03 -0.30144437E-03  0.15042722E-04</point>
      <point r="9.460">  0.43785162E-04 -0.11944094E-03 -0.29295534E-03  0.14382641E-04</point>
      <point r="9.480">  0.42059755E-04 -0.11594821E-03 -0.28474915E-03  0.13742953E-04</point>
      <point r="9.500">  0.40379375E-04 -0.11256453E-03 -0.27681450E-03  0.13122948E-04</point>
      <point r="9.520">  0.38742367E-04 -0.10928551E-03 -0.26914049E-03  0.12521941E-04</point>
      <point r="9.540">  0.37147149E-04 -0.10610691E-03 -0.26171666E-03  0.11939279E-04</point>
      <point r="9.560">  0.35592208E-04 -0.10302466E-03 -0.25453291E-03  0.11374332E-04</point>
      <point r="9.580">  0.34076098E-04 -0.10003488E-03 -0.24757954E-03  0.10826500E-04</point>
      <point r="9.600">  0.32597438E-04 -0.97133822E-04 -0.24084723E-03  0.10295204E-04</point>
      <point r="9.620">  0.31154908E-04 -0.94317910E-04 -0.23432700E-03  0.97798921E-05</point>
      <point r="9.640">  0.29747247E-04 -0.91583712E-04 -0.22801023E-03  0.92800334E-05</point>
      <point r="9.660">  0.28373252E-04 -0.88927937E-04 -0.22188862E-03  0.87951193E-05</point>
      <point r="9.680">  0.27031773E-04 -0.86347438E-04 -0.21595422E-03  0.83246621E-05</point>
      <point r="9.700">  0.25721715E-04 -0.83839195E-04 -0.21019937E-03  0.78681945E-05</point>
      <point r="9.720">  0.24442032E-04 -0.81400322E-04 -0.20461674E-03  0.74252682E-05</point>
      <point r="9.740">  0.23191725E-04 -0.79028055E-04 -0.19919926E-03  0.69954534E-05</point>
      <point r="9.760">  0.21969843E-04 -0.76719747E-04 -0.19394018E-03  0.65783381E-05</point>
      <point r="9.780">  0.20775479E-04 -0.74472870E-04 -0.18883300E-03  0.61735272E-05</point>
      <point r="9.800">  0.19607769E-04 -0.72285004E-04 -0.18387151E-03  0.57806417E-05</point>
      <point r="9.820">  0.18465889E-04 -0.70153836E-04 -0.17904972E-03  0.53993182E-05</point>
      <point r="9.840">  0.17349054E-04 -0.68077156E-04 -0.17436193E-03  0.50292080E-05</point>
      <point r="9.860">  0.16256516E-04 -0.66052852E-04 -0.16980265E-03  0.46699767E-05</point>
      <point r="9.880">  0.15187564E-04 -0.64078907E-04 -0.16536663E-03  0.43213035E-05</point>
      <point r="9.900">  0.14141520E-04 -0.62153396E-04 -0.16104885E-03  0.39828805E-05</point>
      <point r="9.920">  0.13117738E-04 -0.60274480E-04 -0.15684449E-03  0.36544122E-05</point>
      <point r="9.940">  0.12115606E-04 -0.58440405E-04 -0.15274896E-03  0.33356149E-05</point>
      <point r="9.960">  0.11134538E-04 -0.56649497E-04 -0.14875784E-03  0.30262163E-05</point>
      <point r="9.980">  0.10173980E-04 -0.54900162E-04 -0.14486692E-03  0.27259549E-05</point>
      <point r="10.000">  0.92334026E-05 -0.53190877E-04 -0.14107219E-03  0.24345796E-05</point>
      <point r="10.020">  0.83123033E-05 -0.51520193E-04 -0.13736979E-03  0.21518490E-05</point>
      <point r="10.040">  0.74102041E-05 -0.49886729E-04 -0.13375606E-03  0.18775312E-05</point>
      <point r="10.060">  0.65266505E-05 -0.48289169E-04 -0.13022748E-03  0.16114033E-05</point>
      <point r="10.080">  0.56612107E-05 -0.46726260E-04 -0.12678070E-03  0.13532509E-05</point>
      <point r="10.100">  0.48134737E-05 -0.45196811E-04 -0.12341254E-03  0.11028678E-05</point>
      <point r="10.120">  0.39830492E-05 -0.43699687E-04 -0.12011994E-03  0.86005572E-06</point>
      <point r="10.140">  0.31695660E-05 -0.42233809E-04 -0.11690000E-03  0.62462358E-06</point>
      <point r="10.160">  0.23726713E-05 -0.40798151E-04 -0.11374996E-03  0.39638757E-06</point>
      <point r="10.180">  0.15920296E-05 -0.39391737E-04 -0.11066718E-03  0.17517056E-06</point>
      <point r="10.200">  0.82732204E-06 -0.38013641E-04 -0.10764914E-03 -0.39198142E-07</point>
      <point r="10.220">  0.78245384E-07 -0.36662983E-04 -0.10469347E-03 -0.24688308E-06</point>
      <point r="10.240"> -0.65548875E-06 -0.35338925E-04 -0.10179787E-03 -0.44804299E-06</point>
      <point r="10.260"> -0.13741547E-05 -0.34040675E-04 -0.98960206E-04 -0.64283101E-06</point>
      <point r="10.280"> -0.20780133E-05 -0.32767477E-04 -0.96178409E-04 -0.83139503E-06</point>
      <point r="10.300"> -0.27673131E-05 -0.31518618E-04 -0.93450530E-04 -0.10138779E-05</point>
      <point r="10.320"> -0.34422905E-05 -0.30293418E-04 -0.90774714E-04 -0.11904176E-05</point>
      <point r="10.340"> -0.41031708E-05 -0.29091235E-04 -0.88149203E-04 -0.13611476E-05</point>
      <point r="10.360"> -0.47501685E-05 -0.27911458E-04 -0.85572327E-04 -0.15261971E-05</point>
      <point r="10.380"> -0.53834883E-05 -0.26753509E-04 -0.83042502E-04 -0.16856910E-05</point>
      <point r="10.400">  0.00000000E+00</point>
    </S_spline>
    <Vrep_spline>
      <point r="0.000">  0.38236698E+02</point>
      <point r="0.200">  0.23847193E+02</point>
      <point r="0.400">  0.14884652E+02</point>
      <point r="0.600">  0.93023087E+01</point>
      <point r="0.800">  0.58253317E+01</point>
      <point r="1.000">  0.36596872E+01</point>
      <point r="1.200">  0.23108098E+01</point>
      <point r="1.400">  0.14706579E+01</point>
      <point r="1.600">  0.94736706E+00</point>
      <point r="1.800">  0.62143396E+00</point>
      <point r="2.000">  0.41842566E+00</point>
      <point r="2.056">  0.37697051E+00</point>
      <point r="2.112">  0.34140709E+00</point>
      <point r="2.168">  0.31079383E+00</point>
      <point r="2.224">  0.28431428E+00</point>
      <point r="2.280">  0.26126506E+00</point>
      <point r="2.336">  0.24104454E+00</point>
      <point r="2.392">  0.22314225E+00</point>
      <point r="2.448">  0.20712897E+00</point>
      <point r="2.504">  0.19264757E+00</point>
      <point r="2.560">  0.17940440E+00</point>
      <point r="2.616">  0.16716145E+00</point>
      <point r="2.672">  0.15572900E+00</point>
      <point r="2.728">  0.14495896E+00</point>
      <point r="2.784">  0.13473869E+00</point>
      <point r="2.840">  0.12498547E+00</point>
      <point r="2.896">  0.11564138E+00</point>
      <point r="2.952">  0.10666873E+00</point>
      <point r="3.008">  0.98046016E-01</point>
      <point r="3.064">  0.89764245E-01</point>
      <point r="3.120">  0.81823733E-01</point>
      <point r="3.176">  0.74231306E-01</point>
      <point r="3.232">  0.66997891E-01</point>
      <point r="3.288">  0.60136453E-01</point>
      <point r="3.344">  0.53660288E-01</point>
      <point r="3.400">  0.47581623E-01</point>
      <point r="3.456">  0.41910508E-01</point>
      <point r="3.512">  0.36653983E-01</point>
      <point r="3.568">  0.31815491E-01</point>
      <point r="3.624">  0.27394515E-01</point>
      <point r="3.680">  0.23386413E-01</point>
      <point r="3.736">  0.19782437E-01</point>
      <point r="3.792">  0.16569907E-01</point>
      <point r="3.848">  0.13732512E-01</point>
      <point r="3.904">  0.11250728E-01</point>
      <point r="3.960">  0.91023201E-02</point>
      <point r="4.016">  0.72629048E-02</point>
      <point r="4.072">  0.57065605E-02</point>
      <point r="4.128">  0.44064510E-02</point>
      <point r="4.184">  0.33354454E-02</point>
      <point r="4.240">  0.24667113E-02</point>
      <point r="4.296">  0.17742565E-02</point>
      <point r="4.352">  0.12333973E-02</point>
      <point r="4.408">  0.82113067E-03</point>
      <point r="4.464">  0.51638707E-03</point>
      <point r="4.520">  0.30014058E-03</point>
      <point r="4.576">  0.15535432E-03</point>
      <point r="4.632">  0.66737588E-04</point>
      <point r="4.688">  0.20292079E-04</point>
      <point r="4.744">  0.26240937E-05</point>
      <point r="4.800">  0.00000000E+00</point>
    </Vrep_spline>
  </per_pair_data>
</DFTB_params>
</eval_test_params>"""

   got_tight_binding = True
   try:
      p = Potential('TB NRL-TB', param_str=NRL_TB_tight_binding_xml)
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

            self.pot = Potential('TB NRL-TB', param_str=NRL_TB_tight_binding_xml)
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
                                                [ -1.44992225,    1.05479624,    0.78616438]]).T)

   got_tight_binding = True
   try:
      p = Potential('TB DFTB', param_str=DFTB_tight_binding_xml)
   except RuntimeError:
      got_tight_binding = False

   if got_tight_binding:
      class TestPotential_DFTB(QuippyTestCase):

         def setUp(self):
            xyz = """8
   Lattice="5.42883523318981 0 0 0 5.42883523318981 0 0 0 5.42883523318981" Properties=Z:I:1:pos:R:3
   6 0.0210809907025043 -0.082432438809103 -0.0525403165481939
   14 -0.0194362141345624 2.79575394428175 2.65958482185915
   14 2.6322497362192 0.0428914772326525 2.75195107653885
   14 2.7348145379646 2.61668595724592 -0.00569985609748024
   14 1.33742072200028 1.26492764509989 1.32415402710619
   14 1.27192747372937 4.12374105659993 4.07496695423226
   14 3.97244199589577 1.36902339889138 4.0668447417454
   14 4.09570476049115 4.02286216155155 1.27329051246382"""

            self.pot = Potential('TB DFTB', param_str=DFTB_tight_binding_xml)
            self.at = Atoms(xyz, format='string')

            verbosity_push(PRINT_SILENT)

         def tearDown(self):
            verbosity_pop()

         def test_energy(self):
            self.pot.calc(self.at, args_str="energy SCF_NONLOCAL_U_DFTB=T")
            self.assertAlmostEqual(self.at.energy, -290.32759387387028482)

         def test_virial(self):
            self.pot.calc(self.at, SCF_NONLOCAL_U_DFTB=True, args_str="virial")
            self.assertArrayAlmostEqual(self.at.virial,
				        farray([[ -0.84692881758198624E+01,  0.50924629468887850E+00,  0.14608485544479166E+01], 
					        [ 0.50924629468888027E+00, -0.85692585595856414E+01,  0.20723041137007878E+00], 
					        [ 0.14608485544479148E+01,  0.20723041137007966E+00, -0.77367535676847403E+01]]))

         def test_force(self):
            self.pot.calc(self.at, SCF_NONLOCAL_U_DFTB=True, args_str="force")
            self.assertArrayAlmostEqual(self.at.force,
				        farray([[  -0.65034095, 0.54765845, 0.23744842],
					        [  0.02870700, -1.58253909, 0.44340254],
					        [  0.97285025, -1.23003376, -0.84897610],
					        [  -1.38141608, 1.29919950, 0.13578729],
					        [  -2.16221795, -0.57459449, -1.50507678],
					        [  -0.64930114, 0.54914774, 1.99718808],
					        [  3.31926807, -1.67338278, 0.49712349],
					        [  0.52245078, 2.66454443, -0.95689694]]).T)


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
         self.assertArrayAlmostEqual(self.a.force, 1.0*np.ones((3,self.a.n)))

      def test_virial(self):
         self.p.calc(self.a, args_str="virial")
         self.assertArrayAlmostEqual(self.a.virial, 1.0*np.ones((3,3)))

      def test_all(self):
         e = farray(0.0)
         f = fzeros((3,self.a.n))
         v = fzeros((3,3))
         self.p.calc(self.a, energy=e, force=f, virial=v)
         self.assertAlmostEqual(e, 1.0)
         self.assertArrayAlmostEqual(v, 1.0*np.ones((3,3)))
         self.assertArrayAlmostEqual(f, 1.0*np.ones((3,self.a.n)))

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
         randomise(self.at.pos, 0.1)
         self.at.set_cutoff(self.pot.cutoff())
         self.at.calc_connect()

         self.nsteps_ref = 6

         self.lat_ref = farray([[  5.43118500e+00,   6.52383282e-04,  -7.71644552e-05],
                                [  6.52353394e-04,   5.43113895e+00,  -1.05226108e-04],
                                [ -7.71047017e-05,  -1.05266395e-04,   5.43111965e+00]]).T

         self.pos_ref = farray([[ -1.46951886e-02,  -4.31724853e-04,  -2.98041952e-03],
                                [  1.34330077e+00,   1.35748083e+00,   1.35458260e+00],
                                [  2.70117370e+00,   2.71545008e+00,  -3.04193919e-03],
                                [ -1.37266267e+00,  -1.35840607e+00,   1.35462470e+00],
                                [  2.70085164e+00,  -1.50958027e-04,   2.71252244e+00],
                                [ -1.37227962e+00,   1.35724737e+00,  -1.36097640e+00],
                                [ -1.43610167e-02,   2.71511882e+00,   2.71250787e+00],
                                [  1.34296689e+00,  -1.35806961e+00,  -1.36098177e+00]]).T

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

      if (hasattr(quippy.Potential,'calc_elastic_constants')):
	 def testelastic_constants(self):
	    self.pot.calc_elastic_constants(self.at, fd=1e-3, c=self.c, c0=self.c0, relax_initial=True)
	    self.assertArrayAlmostEqual(self.c, self.c_ref)
	    self.assertArrayAlmostEqual(self.c0, self.c0_ref)


if __name__ == '__main__':
   unittest.main()
