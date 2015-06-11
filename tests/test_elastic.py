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
import quippy._elasticity
from quippy.elasticity import *
from numpy import *
import unittest, itertools, sys, quippy, os, numpy
from quippytest import *


class TestElasticFieldsCubic(QuippyTestCase):

  def setUp(self):

     self.C = farray([[ 151.4276439 ,   76.57244456,   76.57244456,    0.        ,    0.        ,    0.        ],
                      [  76.57244456,  151.4276439 ,   76.57244456,    0.        ,    0.        ,    0.        ],
                      [  76.57244456,   76.57244456,  151.4276439 ,    0.        ,    0.        ,    0.        ],
                      [   0.        ,    0.        ,    0.        ,  109.85498798,    0.        ,    0.        ],
                      [   0.        ,    0.        ,    0.        ,    0.        ,  109.85498798,    0.        ],
                      [   0.        ,    0.        ,    0.        ,    0.        ,   0.        ,  109.85498798]])

     self.C_relaxed = farray([[ 151.28712587,   76.5394162 ,   76.5394162 ,    0.        ,    0.        ,    0.        ],
                              [  76.5394162 ,  151.28712587,   76.5394162 ,    0.        ,    0.        ,    0.        ],
                              [  76.5394162 ,   76.5394162 ,  151.28712587,    0.        ,    0.        ,    0.        ],
                              [   0.        ,    0.        ,    0.        ,   56.32421772,    0.        ,    0.        ],
                              [   0.        ,    0.        ,    0.        ,    0.        ,   56.32421772,    0.        ],
                              [   0.        ,    0.        ,    0.        ,    0.        ,    0.        ,   56.32421772]])


     # Relaxed lattice constant
     self.a = 5.43094977

     self.at0 = diamond(self.a, 14)
     self.eps = 1e-3

     self.xml="""
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

     self.pot = Potential('IP SW', param_str=self.xml)
     verbosity_push(PRINT_SILENT)

  def tearDown(self):
     verbosity_pop()

  def check_strain(self, pattern, C, relax=False):

     # Apply desired strain
     pattern = farray(pattern)
     applied_strain = numpy.where(pattern == 1, self.eps, 0.0)
     strain_m = strain_matrix(applied_strain)
     at = transform(self.at0, strain_m)

     # Compute virial stress for comparison
     at.set_cutoff(self.pot.cutoff())
     at.calc_connect()

     if relax:
        # Minimise internal degrees of freedom
        self.pot.minim(at, 'cg', 1e-6, 100, do_pos=True, do_lat=False)

     self.pot.calc(at, virial=True)
     virial_stress = stress_vector(-at.virial*GPA/at.cell_volume())

     quippy._elasticity.elastic_fields(at, a=self.a, cij=C)

     for i in frange(at.n):
        # Check computed strain matches applied strain
        strain_i = farray([at.s_xx_sub1[i], at.s_yy_sub1[i], at.s_zz_sub1[i], at.s_yz[i], at.s_xz[i], at.s_xy[i]])
        self.assertArrayAlmostEqual(strain_i, applied_strain)

        # Check stress matches dot(C, applied_strain)
        stress_i = farray([at.sig_xx[i], at.sig_yy[i], at.sig_zz[i], at.sig_yz[i], at.sig_xz[i], at.sig_xy[i]])
        self.assertArrayAlmostEqual(stress_i, dot(C, applied_strain))

        # Check calculated stress matches computed virial stress
        self.assertArrayAlmostEqual(stress_i, virial_stress, tol=self.eps)

  def test1(self):
     self.check_strain([1,0,0,0,0,0], self.C, False)

  def test2(self):
     self.check_strain([0,1,0,0,0,0], self.C, False)

  def test3(self):
     self.check_strain([0,0,1,0,0,0], self.C, False)

  def test4(self):
     self.check_strain([0,0,0,1,0,0], self.C, False)

  def test5(self):
     self.check_strain([0,0,0,0,1,0], self.C, False)

  def test6(self):
     self.check_strain([0,0,0,0,0,1], self.C, False)

  def test14(self):
     self.check_strain([1,0,0,1,0,0], self.C, False)

  def test15(self):
     self.check_strain([1,0,0,0,1,0], self.C, False)

  def test16(self):
     self.check_strain([1,0,0,0,0,1], self.C, False)

  def test1_relaxed(self):
     self.check_strain([1,0,0,0,0,0], self.C_relaxed, True)

  def test2_relaxed(self):
     self.check_strain([0,1,0,0,0,0], self.C_relaxed, True)

  def test3_relaxed(self):
     self.check_strain([0,0,1,0,0,0], self.C_relaxed, True)

  def test4_relaxed(self):
     self.check_strain([0,0,0,1,0,0], self.C_relaxed, True)

  def test5_relaxed(self):
     self.check_strain([0,0,0,0,1,0], self.C_relaxed, True)

  def test6_relaxed(self):
     self.check_strain([0,0,0,0,0,1], self.C_relaxed, True)

  def test14_relaxed(self):
     self.check_strain([1,0,0,1,0,0], self.C_relaxed, True)

  def test15_relaxed(self):
     self.check_strain([1,0,0,0,1,0], self.C_relaxed, True)

  def test16_relaxed(self):
     self.check_strain([1,0,0,0,0,1], self.C_relaxed, True)


if __name__ == '__main__':
   unittest.main()
