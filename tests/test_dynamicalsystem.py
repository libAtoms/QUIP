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

"""
This module defines the :class:`DynamicalSystem` and :class:`Dynamics`
classes, for carrying out molecular dynamics simulations.

.. module:: quippy.dynamicalsystem
   :synopsis: Run molecular dynamics simulations
"""

from quippy import *
from numpy import *

import unittest, quippy
from quippytest import *

class TestDynamicalSystem(QuippyTestCase):
   def setUp(self):
       d = diamond(5.44, 14)
       self.at = supercell(d, 2, 2, 2)
       self.ds = DynamicalSystem(self.at)

   def test_atoms(self):
       self.assertEqual(self.ds.atoms, self.at)

   def test_avgpos(self):
       self.assertArrayAlmostEqual(self.ds.atoms.avgpos, self.ds.atoms.pos)

   def test_avg_ke(self):
       self.assertArrayAlmostEqual(self.ds.atoms.avg_ke, 0.5*self.ds.atoms.mass*self.ds.atoms.velo.norm2())

   def test_set_masses(self):
       """
       Test for regressions to GitHub issue #25
       """
       at = self.at.copy()
       at.add_property('mass', 0.0)
       at.mass[:] = at.get_masses()*MASSCONVERT
       mass1 = at.mass.copy()
       ds = DynamicalSystem(at) # should not change masses
       self.assertArrayAlmostEqual(mass1, ds.atoms.mass)

if 'ase' in available_modules:
   import ase
else:
   import quippy.miniase as ase

class GenericTestDynamics(object):
   def common_init(self):
      random.seed(1)
      system_reseed_rng(2065775975)
      self.at1 = diamond(5.44, 14)
      self.at1.rattle(0.01)
      self.at2 = ase.Atoms(self.at1)
      self.dt = 1.0 # fs
      self.skin = 2.0
      self.pot = Potential('IP SW', cutoff_skin=self.skin)
      self.orig_at = self.at1.copy()
      self.at1.set_calculator(self.pot)
      self.at2.set_calculator(self.pot)
      self.ds = DynamicalSystem(self.at1)
      self.dyn = Dynamics(self.at2, self.dt*fs, trajectory=None, logfile=None)


   def test_final_state(self):
      #print 'pos_ref = %r' % self.ds.atoms.get_positions()
      #print
      #print 'velo_ref = %r' % self.ds.atoms.get_velocities()

      self.assertArrayAlmostEqual(self.ds.atoms.get_positions(), self.pos_ref)
      self.assertArrayAlmostEqual(self.dyn.atoms.get_positions(), self.pos_ref)
      self.assertArrayAlmostEqual(self.ds.atoms.velo.T*sqrt(MASSCONVERT), self.velo_ref)
      self.assertArrayAlmostEqual(self.dyn.atoms.get_velocities(), self.velo_ref)



class TestDynamics_1step(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[  4.97540436e-03,  -1.38126114e-03,   6.45355088e-03],
       [  1.37520514e+00,   1.35768022e+00,   1.35766373e+00],
       [  2.73574542e+00,   2.72766112e+00,  -4.70585795e-03],
       [  4.08545164e+00,   4.07535031e+00,   1.35535781e+00],
       [  2.72242689e+00,  -1.90944570e-02,   2.70278269e+00],
       [  4.07437528e+00,   1.34984705e+00,   4.08313749e+00],
       [ -9.05352820e-03,   2.70588380e+00,   2.73460483e+00],
       [  1.35774766e+00,   4.08066033e+00,   4.06579153e+00]])

   velo_ref = array([[  1.67914860e-04,   2.79012558e-05,  -4.74020628e-04],
       [ -5.11497114e-04,   4.41792526e-04,   1.03393102e-04],
       [ -9.48979802e-04,  -2.69231088e-04,  -2.25424519e-04],
       [  5.28719044e-04,  -3.14938172e-04,   3.06548847e-04],
       [  1.47522461e-04,   7.79152955e-04,   6.47960806e-04],
       [ -3.70482022e-05,  -4.99786924e-04,  -1.01582505e-04],
       [  5.43099683e-04,   1.38849570e-04,  -1.04927134e-03],
       [  1.10269070e-04,  -3.03740122e-04,   7.92396240e-04]])

   def setUp(self):
      self.common_init()

      # advance one step with the DynamicalSystem + quippy Atoms
      self.ds.atoms.set_cutoff(self.pot.cutoff(), cutoff_skin=self.skin)
      self.ds.atoms.calc_connect()

      # use initial forces to set initial acceleration
      self.pot.calc(self.ds.atoms, args_str='force')
      self.ds.atoms.acc[...] = self.ds.atoms.force/self.ds.atoms.mass

      self.ds.advance_verlet1(self.dt)
      self.pot.calc(self.ds.atoms, args_str='force')
      self.ds.advance_verlet2(self.dt, self.ds.atoms.force)

      # Advance one step with the Dynamics wrapper
      f0 = self.dyn.atoms.get_forces()
      self.f1 = self.dyn.step(f0)


class TestDynamics_50steps(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[  3.84004892e-03,  -6.83286467e-03,  -3.84833019e-03],
       [  1.35831451e+00,   1.35029822e+00,   1.35533926e+00],
       [  2.72050540e+00,   2.70597514e+00,   1.17276445e-03],
       [  4.07926556e+00,   4.07607935e+00,   1.35518719e+00],
       [  2.71969592e+00,  -1.77990479e-03,   2.72309193e+00],
       [  4.09123448e+00,   1.36339124e+00,   4.07478637e+00],
       [  9.79070646e-03,   2.71688867e+00,   2.71690877e+00],
       [  1.36422729e+00,   4.07258727e+00,   4.07844782e+00]])

   velo_ref = array([[ -5.14487324e-03,  -2.82951901e-03,   7.20871809e-04],
       [  1.53645902e-03,  -1.11169548e-02,   1.46838874e-05],
       [  8.72668825e-03,  -2.62441606e-03,   3.88577807e-03],
       [ -1.07551841e-02,   4.18162354e-03,  -1.83166095e-03],
       [ -1.11052239e-03,  -3.68317014e-03,  -2.84833123e-03],
       [  6.19669513e-03,   1.02391352e-02,  -4.67761315e-04],
       [ -2.05447112e-03,   5.34874634e-03,   4.60664348e-03],
       [  2.60520842e-03,   4.84554900e-04,  -4.08022375e-03]])

   def setUp(self):
      self.common_init()

      # advance 50 steps with the DynamicalSystem + quippy Atoms
      verbosity_push(PRINT_SILENT)
      self.ds.run(self.pot, dt=self.dt, n_steps=50)
      verbosity_pop()

      # Advance 50 steps with the Dynamics wrapper
      self.dyn.run(50)


class TestDynamics_InitTemperature(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[ 0.05803293, -0.04190898,  0.00749261],
       [ 1.34873075,  1.40383754,  1.36193669],
       [ 2.74431837,  2.72345343, -0.01449668],
       [ 4.09026606,  4.06058596,  1.42717428],
       [ 2.63832699, -0.01874375,  2.64817749],
       [ 4.0830447 ,  1.31413232,  4.04176328],
       [ 0.01976061,  2.75224012,  2.78688954],
       [ 1.3643935 ,  4.08301046,  4.04214856]])

   velo_ref = array([[ 0.00601802, -0.01855572,  0.00017378],
       [ 0.00595823, -0.02294915,  0.02800107],
       [-0.02460226,  0.01526761,  0.01087745],
       [-0.00900428,  0.00726496,  0.01150873],
       [-0.01207381, -0.02462606,  0.01587958],
       [ 0.01408524,  0.02512267, -0.02771845],
       [ 0.02072858,  0.02207973, -0.00283827],
       [-0.00110974, -0.00360404, -0.0358839 ]])

   def setUp(self):
      self.common_init()

      # choose some random velocities for one case
      self.dyn.set_temperature(300.)

      # use same random initial velocities
      self.ds.atoms.velo[...] = transpose(self.dyn.atoms.get_velocities()/sqrt(MASSCONVERT))

      # advance 50 steps with the DynamicalSystem + quippy Atoms
      verbosity_push(PRINT_SILENT)
      self.ds.run(self.pot, dt=self.dt, n_steps=50)
      verbosity_pop()

      # Advance 50 steps with the Dynamics wrapper
      self.dyn.run(50)


class TestDynamics_Thermostat(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[  1.04158565e-02,  -9.87095916e-04,   4.63136077e-04],
       [  1.37069442e+00,   1.35732272e+00,   1.36231257e+00],
       [  2.73317994e+00,   2.72473324e+00,  -2.53174942e-03],
       [  4.09431665e+00,   4.07259175e+00,   1.34964666e+00],
       [  2.72208619e+00,  -1.54534377e-02,   2.70315734e+00],
       [  4.07948878e+00,   1.34575567e+00,   4.07878521e+00],
       [ -3.89977809e-03,   2.71244934e+00,   2.72906460e+00],
       [  1.35610951e+00,   4.08095357e+00,   4.06509949e+00]])

   velo_ref = array([[ 0.01003826, -0.00199035, -0.01122275],
       [-0.00438354,  0.00514708,  0.01096993],
       [-0.00706211, -0.00601205,  0.00248372],
       [ 0.01964998, -0.00356976, -0.00803951],
       [-0.0028455 ,  0.00181909, -0.00144312],
       [ 0.00515434, -0.00648112, -0.01035969],
       [ 0.00507569,  0.00888286, -0.00847751],
       [ 0.00127138, -0.00010426,  0.00312442]])

   def setUp(self):
      self.common_init()

      system_reseed_rng(2065775975)
      self.ds.add_thermostat(THERMOSTAT_LANGEVIN, 300.0, tau=500.0)

      # advance 10 steps with the DynamicalSystem + quippy Atoms
      verbosity_push(PRINT_SILENT)
      self.ds.run(self.pot, dt=self.dt, n_steps=10)
      verbosity_pop()

      # we need same sequence of thermostat random forces
      system_reseed_rng(2065775975)
      self.dyn.add_thermostat(THERMOSTAT_LANGEVIN, 300.0, tau=500.0)

      # Advance 10 steps with the Dynamics wrapper
      self.dyn.run(10)


if __name__ == '__main__':
   unittest.main()
