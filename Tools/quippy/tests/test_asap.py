from quippy import *
import unittest, quippy
from quippytest import *

try:
   p1 = Potential('IP ASAP', '')
   got_asap1 = True
except RuntimeError:
   got_asap1 = False

try:
   p2 = Potential('IP ASAP2', '')
   got_asap2 = True
except RuntimeError:
   got_asap2 = False

# If True, validate ASAP2 against original ASAP potential.
# Otherwise, we compare to reference data in this file.
do_compare_p1_p2 = False

if do_compare_p1_p2 and not got_asap1:
   # Original ASAP is not available
   do_compare_p1_p2 = False
   
if got_asap2:
   class PotTestMixin:
      def compare_p1_p2(self, at, debug=True, df=False):
         e1 = farray(0.0)
         e2 = farray(0.0)
         f1 = fzeros((3, at.n))
         f2 = fzeros((3, at.n))
         v1 = fzeros((3,3))
         v2 = fzeros((3,3))
         df1 = fzeros((3, at.n))
         df2 = fzeros((3, at.n))
         local_e2 = fzeros((at.n,))

         at.add_property('efield', 0.0, n_cols=3)
         at.add_property('dipoles', 0.0, n_cols=3)
         at.add_property('efield_old1', 0.0, n_cols=3)
         at.add_property('efield_old2', 0.0, n_cols=3)
         at.add_property('efield_old3', 0.0, n_cols=3)

         at1 = at.copy()
         at1.set_cutoff(self.cutoff)
         at1.calc_connect()

         at2 = at.copy()
         at2.set_cutoff(self.cutoff)
         at2.calc_connect()

         if df:
            self.p1.calc(at1, e=e1, f=f1, virial=v1, df=df1, calc_dipoles=True)
            self.p2.calc(at2, e=e2, f=f2, virial=v2, df=df2, local_e=local_e2)
         else:
            self.p1.calc(at1, e=e1, f=f1, virial=v1, calc_dipoles=True)
            self.p2.calc(at2, e=e2, f=f2, virial=v2, local_e=local_e2)


         if debug:
            print 'e1 = ', e1
            print 'e2 = ', e2
            print 'f1 = ', f1
            print 'f2 = ', f2
            print 'v1 = ', v1
            print 'v2 = ', v2
            if hasattr(at1, 'dipoles'):
               print 'dip1 = ', at1.dipoles
            if hasattr(at2, 'dipoles'):
               print 'dip2 = ', at2.dipoles

         self.assertAlmostEqual(e1, e2)
         self.assertAlmostEqual(e2, sum(local_e2))
         self.assertArrayAlmostEqual(f1, f2)
         self.assertArrayAlmostEqual(v1, v2)
         if df:
            self.assertArrayAlmostEqual(df1, df2)
            self.assertArrayAlmostEqual(f1, df1, tol=1e-5)
            self.assertArrayAlmostEqual(f2, df2, tol=1e-5)
         if hasattr(at1, 'dipoles') and hasattr(at2, 'dipoles'):
            self.assertArrayAlmostEqual(at1.dipoles, at2.dipoles)

         return e2, f2.T, v2, local_e2, at2.dipoles


      def compare_ref(self, at, ref):
         e = farray(0.0)
         f = fzeros((3, at.n))
         v = fzeros((3,3))
         df = fzeros((3, at.n))
         local_e = fzeros((at.n,))

         at.add_property('efield', 0.0, n_cols=3)
         at.add_property('dipoles', 0.0, n_cols=3)
         at.add_property('efield_old1', 0.0, n_cols=3)
         at.add_property('efield_old2', 0.0, n_cols=3)
         at.add_property('efield_old3', 0.0, n_cols=3)

         at.calc_connect()
         self.p2.calc(at, e=e, f=f, virial=v, local_e=local_e)

         e_ref, f_ref, v_ref, local_e_ref, dip_ref = ref

         self.assertAlmostEqual(e, e_ref)
         self.assertAlmostEqual(e, sum(local_e))
         self.assertArrayAlmostEqual(local_e, local_e_ref)
         self.assertArrayAlmostEqual(f, f_ref)
         self.assertArrayAlmostEqual(v, v_ref)
         self.assertArrayAlmostEqual(at.dipoles, dip_ref)

      def test_dimer(self):
         dimer = Atoms(n=2, lattice=100.0*fidentity(3))
         dimer.pos[1] = [0.0,0.0,0.0]
         dimer.pos[2] = [3.042*BOHR, 0.0, 0.0]
         dimer.set_atoms([14, 8])
         dimer.set_cutoff(self.cutoff)
         if do_compare_p1_p2:
            self.compare_p1_p2(dimer, debug=self.debug, df=False)
         else:
            self.compare_ref(dimer, self.dimer_ref)

      def test_trimer(self):
         trimer = Atoms(n=3, lattice=100.0*fidentity(3))
         trimer.pos[1] = [0.0,0.0,0.0]
         trimer.pos[2] = [3.042*BOHR, 0.0, 0.0]
         trimer.pos[3] = [2.0*3.042*BOHR, 0.0, 0.0]
         trimer.set_atoms([8, 14, 8])
         trimer.set_cutoff(self.cutoff)
         if do_compare_p1_p2:
            self.compare_p1_p2(trimer, debug=self.debug, df=False)
         else:
            self.compare_ref(trimer, self.trimer_ref)

      def test_quartz(self):
         quartz = alpha_quartz(**sio2.quartz_params['CASTEP_LDA'])
         quartz.set_cutoff(self.cutoff)
         if do_compare_p1_p2:
            self.compare_p1_p2(quartz, debug=self.debug)
         else:
            self.compare_ref(quartz, self.quartz_ref)

      def __test_bigquartz(self):
         """This test fails due to round-off errors in original ASAP implementation
            If two atoms' positions differ by exactly half a unit cell, i.e. |s_i - s_j| = 0.5,
            then conversion from scaled to absolute coordinates is unstable.
            (short_range.f:229, nnlist.f:48)"""
         quartz = alpha_quartz(**sio2.quartz_params['CASTEP_LDA'])
         bigquartz = supercell(quartz, 2, 1, 1)
         bigquartz.set_cutoff(self.cutoff)

         if do_compare_p1_p2:
            self.compare_p1_p2(bigquartz, debug=self.debug)
         else:
            self.compare_ref(bigquartz, self.bigquartz_ref)

      def test_bigquartz_randomise(self):
         quartz = alpha_quartz(**sio2.quartz_params['CASTEP_LDA'])
         bigquartz = supercell(quartz, 2, 1, 1)
         numpy.random.seed(1)
         bigquartz.pos += numpy.random.uniform(-0.1,0.1,size=3*bigquartz.n).reshape(3,bigquartz.n)
         bigquartz.set_cutoff(self.cutoff)
         if do_compare_p1_p2:
            self.compare_p1_p2(bigquartz, debug=self.debug)
         else:
            self.compare_ref(bigquartz, self.bigquartz_ref)


   class TestMorseStretch(QuippyTestCase, PotTestMixin):
      """ Turn off charge and dipole terms, to test short-range part of potential only"""

      dimer_ref =  (FortranArray(1.3232765968383859),
                    FortranArray([[-9.05935698,  0.        ,  0.        ],
                                  [ 9.05935698,  0.        ,  0.        ]]),
                    FortranArray([[ 14.58336505,   0.        ,   0.        ],
                                  [  0.        ,   0.        ,   0.        ],
                                  [  0.        ,   0.        ,   0.        ]]),
                    FortranArray([ 0.6616383,  0.6616383]),
                    FortranArray([[ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.]]))

      trimer_ref = (FortranArray(2.6914040732260349),
                    FortranArray([[-9.23936758,  0.        ,  0.        ],
                                  [ 0.        ,  0.        ,  0.        ],
                                  [ 9.23936758,  0.        ,  0.        ]]),
                    FortranArray([[ 29.7462768,   0.       ,   0.       ],
                                  [  0.       ,   0.       ,   0.       ],
                                  [  0.       ,   0.       ,   0.       ]]),
                    FortranArray([ 0.68406374,  1.3232766 ,  0.68406374]),
                    FortranArray([[ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.]]))

      quartz_ref = (FortranArray(20.41542250348801),
                    FortranArray([[  1.81666108e-01,  -3.14654930e-01,   1.35239216e-14],
                                  [  1.81666108e-01,   3.14654930e-01,   1.23107109e-14],
                                  [ -3.63332217e-01,  -3.73615451e-15,  -2.09875030e-16],
                                  [ -2.81150866e+00,   5.89289282e+00,  -3.78417660e+00],
                                  [ -3.69764055e+00,  -5.38128433e+00,  -3.78417660e+00],
                                  [  6.50914921e+00,  -5.11608487e-01,  -3.78417660e+00],
                                  [ -2.81150866e+00,  -5.89289282e+00,   3.78417660e+00],
                                  [  6.50914921e+00,   5.11608487e-01,   3.78417660e+00],
                                  [ -3.69764055e+00,   5.38128433e+00,   3.78417660e+00]]),
                    FortranArray([[  7.34602052e+01,  -1.54217752e-14,  -1.11188130e-14],
                                  [ -1.54217752e-14,   7.34602052e+01,  -2.65643286e-14],
                                  [ -1.11188130e-14,  -2.65643286e-14,   7.44956117e+01]]),
                    FortranArray([ 1.9532191 ,  1.9532191 ,  1.9532191 ,  2.42596087,  2.42596087,
                                   2.42596087,  2.42596087,  2.42596087,  2.42596087]),
                    FortranArray([[ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.],
                                  [ 0.,  0.,  0.]]))

      bigquartz_ref = (FortranArray(48.876005374180011),
                       FortranArray([[  2.93277763,   9.88802741,  -9.49992649],
                                     [-12.29984225,   3.03696733,  -4.61268954],
                                     [  0.32474038,  -1.28242997,  -0.50325171],
                                     [ -6.83023167,   0.39617117,  -7.36381674],
                                     [  7.6798915 ,  -8.27684431, -13.76672841],
                                     [  3.94615918,  -0.08369327,  -2.24527527],
                                     [ -0.33796743, -12.20849207,  -2.44265921],
                                     [  5.93572874,  -1.42901596,   8.40091225],
                                     [  6.62918157,   9.24560364,  12.51037759],
                                     [ -7.01698404,   2.63741145,   3.81348624],
                                     [  7.1530025 ,   5.51186507,  -9.07658237],
                                     [ -5.06562611,   5.24298752,   1.23723701],
                                     [ -2.08984681,   7.55793921,  -2.67194417],
                                     [-17.9623249 ,  -4.44680332,   2.97420223],
                                     [ 10.54089403, -11.89494979,   4.63150636],
                                     [ -1.08234574,  -7.14925412,   2.66404441],
                                     [  7.87211358,  -4.18856031,  10.10962134],
                                     [ -0.32932015,   7.4430703 ,   5.84148648]]),
                       FortranArray([[ 165.26199497,    3.69267784,    2.42134071],
                                     [   3.69267784,  158.95382262,    6.49024662],
                                     [   2.42134071,    6.49024662,  167.64378429]]),
                       FortranArray([ 1.84093094,  2.94860751,  1.58428019,  2.08952622,  3.18160937,
                                      1.82653314,  2.95589861,  2.80423499,  3.3284506 ,  3.27660168,
                                      2.87698788,  2.80079856,  2.63617854,  3.03799766,  2.51515247,
                                      2.5992419 ,  3.74776611,  2.82520902]),
                       FortranArray([[  0.00000000e+000,   0.00000000e+000,   7.95602526e-316],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  7.78282165e-316,   4.94065646e-324,   2.33639533e-310],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  1.05637316e-312,   7.95602447e-316,   0.00000000e+000],
                                     [  0.00000000e+000,   8.07576622e-316,   0.00000000e+000],
                                     [  0.00000000e+000,   4.94065646e-324,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
                                     [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000]]))

      def setUp(self):
         self.xml = """<ASAP_params 
     betapol="0.75" 
     cutoff="20.0 18.0 18.0 18.0"
     cutoff_ms="18.0"
     cutoff_coulomb="20.0"
     tolpol="0.0005" 
     yuksmoothlength="20.0" 
     iesr="2 2 2" 
     a_ew="1e-06" 
     n_types="2" 
     gcut="0.0" 
     pred_order="2" 
     maxipol="60" 
     raggio="0.0" 
     tewald="F" 
     yukalpha="1e-06">

     <per_type_data atomic_num="8" pol="0.0" z="0.0" type="1" />
     <per_type_data atomic_num="14" pol="0.0" z="0.0" type="2" />

     <per_pair_data atnum_i="8"  atnum_j="8"  C_pol="0.46009932"  D_ms="0.00018650185" gamma_ms="11.642637" B_pol="0.87357114" R_ms="8.0465068" />
     <per_pair_data atnum_i="8"  atnum_j="14" C_pol="-1.5091142"  D_ms="0.0053600978"  gamma_ms="10.405794" B_pol="1.977039"   R_ms="4.193651"  />
     <per_pair_data atnum_i="14" atnum_j="14" C_pol="0.0"         D_ms="-0.0021645401" gamma_ms="4.5784138" B_pol="0.0"        R_ms="13.113727" />
   </ASAP_params>
   """
         if do_compare_p1_p2:
            self.p1 = Potential('IP ASAP', self.xml)
         self.p2 = Potential('IP ASAP2', self.xml)
         self.cutoff = 18.0*BOHR
         self.debug = False



   class TestCharge(QuippyTestCase, PotTestMixin):
      """ Turn off short-range and dipole terms, to test charge part of potential only"""

      dimer_ref = (FortranArray(-24.070897337877238), FortranArray([[ 21.07413139,   0.        ,   0.        ],
              [-21.07413139,   0.        ,   0.        ]]), FortranArray([[-33.92423456,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ]]), FortranArray([-12.03544867, -12.03544867]), FortranArray([[ 0.,  0.,  0.],
              [ 0.,  0.,  0.]]))

      trimer_ref = (FortranArray(-44.403166958656648), FortranArray([[ 18.63466235,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [-18.63466235,   0.        ,   0.        ]]), FortranArray([[-59.99456346,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ]]), FortranArray([-10.16613481, -24.07089734, -10.16613481]), FortranArray([[ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.]]))

      quartz_ref = (FortranArray(-181.6944694918098), FortranArray([[ -5.40952808e-01,   9.36957749e-01,   7.09756958e-15],
              [ -5.40952808e-01,  -9.36957749e-01,   8.09518074e-15],
              [  1.08190562e+00,   6.52367311e-15,  -1.95299023e-14],
              [  4.07411523e+00,  -6.85592161e+00,   4.46551787e+00],
              [  3.90034466e+00,   6.95624809e+00,   4.46551787e+00],
              [ -7.97445989e+00,  -1.00326485e-01,   4.46551787e+00],
              [  4.07411523e+00,   6.85592161e+00,  -4.46551787e+00],
              [ -7.97445989e+00,   1.00326485e-01,  -4.46551787e+00],
              [  3.90034466e+00,  -6.95624809e+00,  -4.46551787e+00]]), FortranArray([[ -7.64106751e+01,   5.42848848e-16,  -9.31988864e-15],
              [  5.42848848e-16,  -7.64106751e+01,   2.78505061e-14],
              [ -9.31988864e-15,   2.78505061e-14,  -7.58877866e+01]]), FortranArray([-36.12479391, -36.12479391, -36.12479391, -12.22001463,
              -12.22001463, -12.22001463, -12.22001463, -12.22001463, -12.22001463]), FortranArray([[ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.],
              [ 0.,  0.,  0.]]))

      bigquartz_ref = (FortranArray(-367.16126747929781), FortranArray([[-2.35632642, -2.0209376 ,  3.87904176],
              [ 2.84952011, -0.45094289, -0.53113981],
              [ 2.420267  , -1.71737752,  4.63039802],
              [ 6.93931868, -4.53731247,  4.68397283],
              [-1.1283966 ,  7.02435327,  9.35725181],
              [-5.29444952, -0.69398915,  1.94212592],
              [ 4.23143492,  9.56380949,  0.31116211],
              [-8.9420606 , -0.94686359, -8.73981904],
              [ 0.02200083, -8.73150269, -6.55068154],
              [ 3.67019499,  1.31856074, -1.61204147],
              [-2.22233694, -3.71938351,  3.10773226],
              [ 0.72156838, -4.56370697,  0.0819926 ],
              [ 4.22696576, -4.80580269,  4.27477823],
              [ 9.99852044,  5.66166145,  2.73743777],
              [-8.56627422,  6.52755751, -0.38682621],
              [ 0.88500014,  8.13213328, -5.01001697],
              [-7.85653112,  2.24543026, -7.62948936],
              [ 0.40158417, -8.28568694, -4.54587894]]), FortranArray([[-153.94991269,   -2.45634663,   -5.5806971 ],
              [  -2.45634663, -152.08622231,   -4.68870731],
              [  -5.5806971 ,   -4.68870731, -154.96812059]]), FortranArray([-33.85355213, -37.66036671, -35.88787751, -11.60745322,
              -12.0585835 , -10.22794598, -12.31926618, -12.35573703,
              -13.36160802, -38.11538713, -36.36248279, -36.88969989,
              -12.32753866, -13.23357389, -11.7024293 , -12.1274268 ,
              -14.4787731 , -12.59156565]), FortranArray([[  5.43472210e-323,   2.18550278e-315,   2.05297892e-315],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  2.04652552e-315,   2.17363545e-315,   0.00000000e+000],
              [  4.01707838e-057,   5.64697363e-038,   6.76931563e-043],
              [  1.58687490e-047,   1.08756584e-071,   4.04676264e-086],
              [  7.06264092e-077,   4.91414628e-062,   2.19926270e-313],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  4.05003039e-116,   2.33761622e-310,   2.14844351e-315],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [ -9.91682126e-270,   2.05286699e-315,   2.14843521e-315],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  4.94065646e-324,   3.88558520e-317,   1.53160350e-322],
              [  1.45003458e-047,   3.82048470e-086,   1.58687490e-047],
              [  1.05637317e-312,   7.06264092e-077,   4.91414628e-062],
              [  2.33761623e-310,   0.00000000e+000,   2.04652694e-315],
              [  1.45189716e-047,   3.45374415e-086,   9.05449679e-043],
              [  4.91414628e-062,   2.19926270e-313,   2.16781531e-315]]))

      def setUp(self):
         self.xml = """<ASAP_params 
     betapol="0.75" 
     cutoff="20.0 18.0 18.0 18.0"
     cutoff_ms="18.0"
     cutoff_coulomb="20.0"
     tolpol="0.0005" 
     yuksmoothlength="20.0" 
     iesr="2 2 2" 
     a_ew="1e-06" 
     n_types="2" 
     gcut="0.0" 
     pred_order="2" 
     maxipol="60" 
     raggio="0.0" 
     tewald="F" 
     yukalpha="1e-06">

     <per_type_data atomic_num="8" pol="0.0" z="-1.95853" type="1" />
     <per_type_data atomic_num="14" pol="0.0" z="3.91706" type="2" />

     <per_pair_data atnum_i="8"  atnum_j="8"  C_pol="0.46009932"  D_ms="0.0" gamma_ms="11.642637" B_pol="0.87357114" R_ms="8.0465068" />
     <per_pair_data atnum_i="8"  atnum_j="14" C_pol="-1.5091142"  D_ms="0.0" gamma_ms="10.405794" B_pol="1.977039"   R_ms="4.193651"  />
     <per_pair_data atnum_i="14" atnum_j="14" C_pol="0.0"         D_ms="0.0" gamma_ms="4.5784138" B_pol="0.0"        R_ms="13.113727" />
   </ASAP_params>
   """
         if do_compare_p1_p2:
            self.p1 = Potential('IP ASAP', self.xml)
         self.p2 = Potential('IP ASAP2', self.xml)
         self.cutoff = 20.0*BOHR
         self.debug = False


   class TestDipoleLongRange(QuippyTestCase, PotTestMixin):
      """ Turn off short-range terms, to test charge and dipole parts of potential only"""

      dimer_ref =  (FortranArray(-17.920706383643051), FortranArray([[ 19.45266941,   0.        ,   0.        ],
              [-19.45266941,   0.        ,   0.        ]]), FortranArray([[-31.31407448,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ]]), FortranArray([-9.04196871, -8.87873767]), FortranArray([[ 0.        ,  0.        ,  0.        ],
              [ 0.10953198,  0.        ,  0.        ]]))

      trimer_ref = (FortranArray(-33.77012116015446), FortranArray([[ 17.66745928,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [-17.66745928,   0.        ,   0.        ]]), FortranArray([[-56.88063925,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ]]), FortranArray([ -7.85260528, -18.0649106 ,  -7.85260528]), FortranArray([[-0.10314826,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.10314826,  0.        ,  0.        ]]))
      
      quartz_ref = (FortranArray(-145.00951258119775), FortranArray([[ -4.56597548e-01,   7.90850151e-01,   1.69993886e-15],
              [ -4.56597548e-01,  -7.90850151e-01,  -1.58513935e-17],
              [  9.13195095e-01,   2.46191320e-15,  -7.56265619e-15],
              [  3.73701627e+00,  -6.41941407e+00,   4.19227411e+00],
              [  3.69086753e+00,   6.44605806e+00,   4.19227411e+00],
              [ -7.42788380e+00,  -2.66439852e-02,   4.19227411e+00],
              [  3.73701627e+00,   6.41941407e+00,  -4.19227411e+00],
              [ -7.42788380e+00,   2.66439852e-02,  -4.19227411e+00],
              [  3.69086753e+00,  -6.44605806e+00,  -4.19227411e+00]]), FortranArray([[ -7.33343215e+01,  -5.19794384e-15,   1.59466046e-14],
              [ -5.21772326e-15,  -7.33343215e+01,   3.55057139e-14],
              [  1.60222203e-14,   3.55414169e-14,  -7.32106130e+01]]), FortranArray([-28.2876303 , -28.2876303 , -28.2876303 , -10.02443695,
              -10.02443695, -10.02443695, -10.02443695, -10.02443695, -10.02443695]), FortranArray([[ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [-0.0232381 ,  0.04181159, -0.02749311],
              [-0.02459085, -0.04103058, -0.02749311],
              [ 0.04782895, -0.00078101, -0.02749311],
              [-0.0232381 , -0.04181159,  0.02749311],
              [ 0.04782895,  0.00078101,  0.02749311],
              [-0.02459085,  0.04103058,  0.02749311]]))

      bigquartz_ref =  (FortranArray(-293.89596056205966), FortranArray([[ -2.08556343,  -2.53157461,   4.07065982],
              [  3.67263959,  -0.65153177,  -0.25212019],
              [  1.99745414,  -1.38116121,   4.08378903],
              [  6.43626849,  -4.01107178,   4.61390629],
              [ -1.57184014,   6.63777767,   9.11923249],
              [ -4.98131305,  -0.551726  ,   1.85144973],
              [  3.51302957,   9.30028   ,   0.51941844],
              [ -8.23922897,  -0.88572301,  -8.20950357],
              [ -0.53908634,  -8.09690393,  -6.56122849],
              [  3.63724393,   0.95543137,  -1.65454499],
              [ -2.55779223,  -3.62022476,   3.48268632],
              [  0.96176958,  -4.7668353 ,  -0.09161431],
              [  3.98097028,  -4.57664084,   4.04336359],
              [ 10.0179143 ,   5.25617998,   2.16209757],
              [ -8.21974883,   6.70521117,  -0.7001308 ],
              [  0.88801158,   7.58080766,  -4.64732635],
              [ -7.26778637,   2.39395518,  -7.36287336],
              [  0.3570579 ,  -7.7562498 ,  -4.46726124]]), FortranArray([[-148.19601394,   -2.36596991,   -5.22831181],
              [  -2.36596991, -146.31561484,   -4.64900502],
              [  -5.22831181,   -4.64900502, -149.86176414]]), FortranArray([-26.2926636 , -29.6335351 , -27.89610351,  -9.39674925,
               -9.99340869,  -8.27479928, -10.14468882, -10.23484308,
              -11.11317274, -30.14839652, -28.72089926, -29.19998377,
              -10.15690765, -10.87159232,  -9.46575002,  -9.96148855,
              -12.05372921, -10.33724919]), FortranArray([[ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [-0.04004351,  0.02478763, -0.03143275],
              [ 0.01007097, -0.04307067, -0.05899246],
              [ 0.03190406,  0.00266745, -0.01270899],
              [-0.02042961, -0.06111135, -0.00274793],
              [ 0.05186195,  0.00611388,  0.05287925],
              [ 0.00476878,  0.05182657,  0.04460869],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [-0.02529707,  0.03199534, -0.02651793],
              [-0.06687741, -0.03356669, -0.01308603],
              [ 0.05273026, -0.04459676,  0.00453838],
              [-0.00629478, -0.04849379,  0.03014353],
              [ 0.04890792, -0.01452389,  0.04788236],
              [-0.00444397,  0.05008518,  0.03009146]]))

      def setUp(self):
         self.xml = """<ASAP_params 
     betapol="0.75" 
     cutoff="20.0 18.0 18.0 18.0"
     cutoff_ms="18.0"
     cutoff_coulomb="20.0"
     tolpol="1e-15" 
     yuksmoothlength="20.0" 
     iesr="2 2 2" 
     a_ew="1e-06" 
     n_types="2" 
     gcut="0.0" 
     pred_order="2" 
     maxipol="60" 
     raggio="0.0" 
     tewald="F" 
     yukalpha="0.1"
     tdip_sr="F">

     <per_type_data atomic_num="8" pol="1.0" z="-1.95853" type="1" />
     <per_type_data atomic_num="14" pol="0.0" z="3.91706" type="2" />

     <per_pair_data atnum_i="8"  atnum_j="8"  C_pol="0.46009932"  D_ms="0.0" gamma_ms="11.642637" B_pol="0.87357114" R_ms="8.0465068" />
     <per_pair_data atnum_i="8"  atnum_j="14" C_pol="-1.5091142"  D_ms="0.0" gamma_ms="10.405794" B_pol="1.977039"   R_ms="4.193651"  />
     <per_pair_data atnum_i="14" atnum_j="14" C_pol="0.0"         D_ms="0.0" gamma_ms="4.5784138" B_pol="0.0"        R_ms="13.113727" />
   </ASAP_params>
   """
         if do_compare_p1_p2:
            self.p1 = Potential('IP ASAP', self.xml)
         self.p2 = Potential('IP ASAP2', self.xml)
         self.cutoff = 20.0*BOHR
         self.debug = False


   class TestDipoleShortRange(QuippyTestCase, PotTestMixin):

      dimer_ref = (FortranArray(-17.811006977457581), FortranArray([[ 18.94240662,   0.        ,   0.        ],
              [-18.94240662,   0.        ,   0.        ]]), FortranArray([[-30.49267528,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ]]), FortranArray([-8.93226931, -8.87873767]), FortranArray([[ 0.       ,  0.       ,  0.       ],
              [ 0.0627256,  0.       ,  0.       ]]))

      trimer_ref = (FortranArray(-33.563525193963763), FortranArray([[ 17.19744833,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [-17.19744833,   0.        ,   0.        ]]), FortranArray([[-55.36743225,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ],
              [  0.        ,   0.        ,   0.        ]]), FortranArray([ -7.85585209, -17.85182102,  -7.85585209]), FortranArray([[-0.05527469,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.05527469,  0.        ,  0.        ]]))

      quartz_ref = (FortranArray(-144.79740179687568), FortranArray([[ -4.92225246e-01,   8.52559136e-01,   2.42526623e-15],
              [ -4.92225246e-01,  -8.52559136e-01,   1.04930082e-15],
              [  9.84450493e-01,   2.37044520e-15,  -7.23613974e-15],
              [  3.78732762e+00,  -6.47630354e+00,   4.23003219e+00],
              [  3.71497958e+00,   6.51807370e+00,   4.23003219e+00],
              [ -7.50230719e+00,  -4.17701600e-02,   4.23003219e+00],
              [  3.78732762e+00,   6.47630354e+00,  -4.23003219e+00],
              [ -7.50230719e+00,   4.17701600e-02,  -4.23003219e+00],
              [  3.71497958e+00,  -6.51807370e+00,  -4.23003219e+00]]), FortranArray([[ -7.29045857e+01,  -5.77736150e-15,   1.57770946e-14],
              [ -5.80585293e-15,  -7.29045857e+01,   3.50407202e-14],
              [  1.57152070e-14,   3.50435723e-14,  -7.26577486e+01]]), FortranArray([-28.21066338, -28.21066338, -28.21066338, -10.02756861,
              -10.02756861, -10.02756861, -10.02756861, -10.02756861, -10.02756861]), FortranArray([[ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [ 0.        ,  0.        ,  0.        ],
              [-0.01017154,  0.01531192, -0.01028628],
              [-0.00817474, -0.01646477, -0.01028628],
              [ 0.01834628,  0.00115285, -0.01028628],
              [-0.01017154, -0.01531192,  0.01028628],
              [ 0.01834628, -0.00115285,  0.01028628],
              [-0.00817474,  0.01646477,  0.01028628]]))

      bigquartz_ref = (FortranArray(-293.30383805940176), FortranArray([[-2.08415672, -2.23007886,  3.74869236],
              [ 3.16451334, -0.59185079, -0.31440948],
              [ 2.04159262, -1.41858379,  4.18632826],
              [ 6.41344975, -4.14899642,  4.54377578],
              [-1.16409577,  6.64240981,  8.88166689],
              [-5.01788994, -0.59216369,  1.88547306],
              [ 3.73907778,  9.16476719,  0.30075182],
              [-8.31315417, -0.92273074, -8.22990847],
              [-0.21573686, -8.15142253, -6.35533748],
              [ 3.3386908 ,  1.24920556, -1.50772306],
              [-2.38418468, -3.48960376,  3.15056973],
              [ 0.91358002, -4.37357124, -0.04218135],
              [ 4.0269014 , -4.67971601,  4.05462354],
              [ 9.64370758,  5.31059828,  2.45489693],
              [-8.08774173,  6.32318443, -0.42988191],
              [ 0.93668411,  7.60451508, -4.71403328],
              [-7.48810435,  2.14860634, -7.22439801],
              [ 0.53686683, -7.84456887, -4.38890531]]), FortranArray([[-147.13389319,   -2.39672914,   -5.20789094],
              [  -2.39672914, -145.18823505,   -4.52055773],
              [  -5.20789094,   -4.52055773, -148.43587042]]), FortranArray([-26.29431962, -29.49983881, -27.84408267,  -9.39524437,
              -10.00607828,  -8.27638007, -10.15057588, -10.24336127,
              -11.11087648, -29.96689858, -28.56601455, -29.07165033,
              -10.16627649, -10.86878858,  -9.47024297,  -9.96695783,
              -12.05319564, -10.35305564]), FortranArray([[  0.00000000e+000,   2.61854792e-322,   0.00000000e+000],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [ -1.29260253e-002,   1.63677234e-002,  -3.84011947e-003],
              [ -8.14834924e-003,  -1.16193389e-002,  -1.26546262e-002],
              [  1.24785221e-002,   2.56647966e-003,  -2.70009262e-003],
              [ -1.39423455e-002,  -1.59807080e-002,   1.77527876e-003],
              [  2.33497524e-002,   7.38173497e-003,   1.88430019e-002],
              [ -9.62827443e-003,   1.65234068e-002,   4.18689651e-003],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [  0.00000000e+000,   0.00000000e+000,   0.00000000e+000],
              [ -1.35080460e-002,   3.33361771e-003,  -1.22585519e-002],
              [ -8.26055231e-003,  -1.37492167e-002,  -1.53129065e-002],
              [  1.32799696e-002,  -4.61174276e-003,  -6.72856511e-003],
              [ -1.06221191e-003,  -1.78793839e-002,   1.50772946e-002],
              [  1.57100029e-002,  -3.03022644e-003,   1.28736848e-002],
              [ -1.61885301e-004,   1.84258788e-002,   7.61293906e-003]]))


      def setUp(self):
         self.xml = """<ASAP_params 
     betapol="0.75" 
     cutoff="20.0 18.0 18.0 18.0"
     cutoff_ms="18.0"
     cutoff_coulomb="20.0"
     tolpol="1e-15" 
     yuksmoothlength="20.0" 
     iesr="2 2 2" 
     a_ew="1e-06" 
     n_types="2" 
     gcut="0.0" 
     pred_order="2" 
     maxipol="60" 
     raggio="0.0" 
     tewald="F" 
     yukalpha="0.1"
     tdip_sr="T">

     <per_type_data atomic_num="8" pol="1.0" z="-1.95853" type="1" />
     <per_type_data atomic_num="14" pol="0.0" z="3.91706" type="2" />

     <per_pair_data atnum_i="8"  atnum_j="8"  C_pol="0.46009932"  D_ms="0.0" gamma_ms="11.642637" B_pol="0.87357114" R_ms="8.0465068" />
     <per_pair_data atnum_i="8"  atnum_j="14" C_pol="-1.5091142"  D_ms="0.0" gamma_ms="10.405794" B_pol="1.977039"   R_ms="4.193651"  />
     <per_pair_data atnum_i="14" atnum_j="14" C_pol="0.0"         D_ms="0.0" gamma_ms="4.5784138" B_pol="0.0"        R_ms="13.113727" />
   </ASAP_params>
   """
         if do_compare_p1_p2:
            self.p1 = Potential('IP ASAP', self.xml)
         self.p2 = Potential('IP ASAP2', self.xml)
         self.cutoff = 20.0*BOHR
         self.debug = False


   if do_compare_p1_p2:
      class TestMD(QuippyTestCase):

         def setUp(self):
            self.xml = """<ASAP_params 
        betapol="0.75" 
        cutoff="20.0 18.0 18.0 18.0"
        cutoff_ms="18.0"
        cutoff_coulomb="20.0"
        tolpol="1e-5" 
        yuksmoothlength="10.0" 
        iesr="2 2 2" 
        a_ew="1e-06" 
        n_types="2" 
        gcut="0.0" 
        pred_order="2" 
        maxipol="60" 
        raggio="0.0" 
        tewald="F" 
        yukalpha="0.1"
        tdip_sr="T">

        <per_type_data atomic_num="8" pol="14.131863" z="-1.4295594" type="1" />
        <per_type_data atomic_num="14" pol="0.0" z="2.8591188" type="2" />

        <per_pair_data C_pol="0.44302622" atnum_j="8" atnum_i="8" D_ms="0.00030700577" gamma_ms="12.165654" B_pol="1.1221903" R_ms="7.0252019" />
        <per_pair_data C_pol="-1.5003213" atnum_j="8" atnum_i="14" D_ms="0.0020129372" gamma_ms="11.350477" B_pol="1.973181" R_ms="4.5780828" />
        <per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="0.33967532" gamma_ms="-0.17694797" B_pol="0.0" R_ms="-0.085202834" />
      </ASAP_params>
      """
            if do_compare_p1_p2:
               self.p1 = Potential('IP ASAP', self.xml)
            self.p2 = Potential('IP ASAP2', self.xml)

            self.p1.print_()
            self.p2.print_()

            self.at = alpha_quartz(**sio2.quartz_params['CASTEP_LDA'])
            self.at.cutoff = 20.0*BOHR

            self.at.add_property('efield', 0.0, n_cols=3)
            self.at.add_property('dipoles', 0.0, n_cols=3)
            self.at.add_property('efield_old1', 0.0, n_cols=3)
            self.at.add_property('efield_old2', 0.0, n_cols=3)
            self.at.add_property('efield_old3', 0.0, n_cols=3)

            self.ds1 = DynamicalSystem(self.at)
            self.ds1.rescale_velo(300.0)
            self.ds1.zero_momentum()
            self.ds1.atoms.calc_connect()
            self.ds1.avg_temp = 0.0

            self.ds2 = DynamicalSystem(self.ds1.atoms.copy())
            self.ds2.atoms.calc_connect()
            self.ds2.avg_temp = 0.0


         def test_md_1step(self):

            self.p1.calc(self.ds1.atoms, calc_force=True)
            self.p2.calc(self.ds2.atoms, calc_force=True)

            self.assertArrayAlmostEqual(self.ds1.atoms.force, self.ds2.atoms.force)

            self.ds1.advance_verlet1(1.0, self.ds1.atoms.force)
            self.ds2.advance_verlet1(1.0, self.ds2.atoms.force)

            self.assertArrayAlmostEqual(self.ds1.atoms.pos, self.ds2.atoms.pos)

            self.p1.calc(self.ds1.atoms, calc_force=True)
            self.p2.calc(self.ds2.atoms, calc_force=True)

            self.assertArrayAlmostEqual(self.ds1.atoms.force, self.ds2.atoms.force)

            self.ds1.advance_verlet2(1.0, self.ds1.atoms.force)
            self.ds2.advance_verlet2(1.0, self.ds1.atoms.force)

            self.assertArrayAlmostEqual(self.ds1.atoms.velo, self.ds2.atoms.velo)

            self.ds1.print_status()
            self.ds2.print_status()


         def test_md_2step(self):

            self.p1.calc(self.ds1.atoms, calc_force=True)
            self.p2.calc(self.ds2.atoms, calc_force=True)

            self.assertArrayAlmostEqual(self.ds1.atoms.force, self.ds2.atoms.force)

            for i in range(2):
               self.ds1.advance_verlet1(1.0, self.ds1.atoms.force)
               self.ds2.advance_verlet1(1.0, self.ds2.atoms.force)

               self.assertArrayAlmostEqual(self.ds1.atoms.pos, self.ds2.atoms.pos)

               self.p1.calc(self.ds1.atoms, calc_force=True)
               self.p2.calc(self.ds2.atoms, calc_force=True)

               self.assertArrayAlmostEqual(self.ds1.atoms.force, self.ds2.atoms.force)

               self.ds1.advance_verlet2(1.0, self.ds1.atoms.force)
               self.ds2.advance_verlet2(1.0, self.ds1.atoms.force)

               self.ds1.print_status()
               self.ds2.print_status()

               self.assertArrayAlmostEqual(self.ds1.atoms.velo, self.ds2.atoms.velo)


         def test_md_10step(self):

            traj1 = self.ds1.run(self.p1, dt=0.5, n_steps=10, save_interval=1)
            traj2 = self.ds2.run(self.p2, dt=0.5, n_steps=10, save_interval=1)

            import itertools

            for i, (at1, at2) in enumerate(itertools.izip(traj1, traj2)):
               self.assertArrayAlmostEqual(at1.force, at2.force)
               self.assertArrayAlmostEqual(at1.pos, at2.pos)

         
   


if __name__ == '__main__':
   unittest.main()
