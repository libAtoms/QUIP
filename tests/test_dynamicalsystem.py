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

import ase

class GenericTestDynamics(object):
   def common_init(self):
      random.seed(1)
      system_reseed_rng(2065775975)
      self.at1 = diamond(5.44, 14)
      # Explicitly set the mass as IUPAC updates in
      # ase will not always sync with quippy and can
      # slighly change the final positions.
      self.at1.set_masses([28.085]*len(self.at1))
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
      #print 'velo_ref = %r' % array(self.ds.atoms.velo.T*sqrt(MASSCONVERT))
      #print

      self.assertArrayAlmostEqual(self.ds.atoms.get_positions(), self.pos_ref)
      self.assertArrayAlmostEqual(self.dyn.atoms.get_positions(), self.pos_ref)
      self.assertArrayAlmostEqual(self.ds.atoms.velo.T*sqrt(MASSCONVERT), self.velo_ref)
      self.assertArrayAlmostEqual(self.dyn.atoms.get_velocities(), self.velo_ref)



class TestDynamics_1step(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[  4.97540451e-03,  -1.38126112e-03,   6.45355046e-03],
       [  1.37520514e+00,   1.35768022e+00,   1.35766373e+00],
       [  2.73574542e+00,   2.72766112e+00,  -4.70585816e-03],
       [  4.08545164e+00,   4.07535031e+00,   1.35535781e+00],
       [  2.72242689e+00,  -1.90944563e-02,   2.70278269e+00],
       [  4.07437528e+00,   1.34984705e+00,   4.08313749e+00],
       [ -9.05352772e-03,   2.70588380e+00,   2.73460483e+00],
       [  1.35774766e+00,   4.08066033e+00,   4.06579153e+00]])

   velo_ref = array([[  1.67917869e-04,   2.79017525e-05,  -4.74029119e-04],
       [ -5.11506283e-04,   4.41800437e-04,   1.03394949e-04],
       [ -9.48996805e-04,  -2.69235922e-04,  -2.25428551e-04],
       [  5.28728510e-04,  -3.14943810e-04,   3.06554333e-04],
       [  1.47525100e-04,   7.79166916e-04,   6.47972422e-04],
       [ -3.70488567e-05,  -4.99795865e-04,  -1.01584331e-04],
       [  5.43109420e-04,   1.38852057e-04,  -1.04929014e-03],
       [  1.10271046e-04,  -3.03745565e-04,   7.92410433e-04]])

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

   pos_ref = array([[  3.83982027e-03,  -6.83299035e-03,  -3.84829791e-03],
       [  1.35831458e+00,   1.35029772e+00,   1.35533926e+00],
       [  2.72050579e+00,   2.70597502e+00,   1.17293724e-03],
       [  4.07926508e+00,   4.07607954e+00,   1.35518711e+00],
       [  2.71969587e+00,  -1.78006885e-03,   2.72309180e+00],
       [  4.09123475e+00,   1.36339170e+00,   4.07478635e+00],
       [  9.79061494e-03,   2.71688891e+00,   2.71690897e+00],
       [  1.36422741e+00,   4.07258729e+00,   4.07844764e+00]])

   velo_ref = array([[ -5.14500838e-03,  -2.82954861e-03,   7.20818721e-04],
       [  1.53654932e-03,  -1.11170527e-02,   1.47967446e-05],
       [  8.72682542e-03,  -2.62427183e-03,   3.88571704e-03],
       [ -1.07552557e-02,   4.18162297e-03,  -1.83155420e-03],
       [ -1.11039864e-03,  -3.68321162e-03,  -2.84845292e-03],
       [  6.19661562e-03,   1.02390323e-02,  -4.67696227e-04],
       [ -2.05461501e-03,   5.34887189e-03,   4.60654485e-03],
       [  2.60528742e-03,   4.84557566e-04,  -4.08017401e-03]])

   def setUp(self):
      self.common_init()

      # advance 50 steps with the DynamicalSystem + quippy Atoms
      verbosity_push(PRINT_SILENT)
      self.ds.run(self.pot, dt=self.dt, n_steps=50)
      verbosity_pop()

      # Advance 50 steps with the Dynamics wrapper
      self.dyn.run(50)


class TestDynamics_InitTemperature(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[ 0.05803303, -0.04190969,  0.00749258],
       [ 1.34873105,  1.40383634,  1.36193793],
       [ 2.74431719,  2.72345406, -0.01449615],
       [ 4.09026562,  4.06058634,  1.42717456],
       [ 2.63832671, -0.01874479,  2.64817844],
       [ 4.08304536,  1.3141336 ,  4.04176215],
       [ 0.01976151,  2.752241  ,  2.78688919],
       [ 1.36439344,  4.08301027,  4.04214707]])

   velo_ref = array([[ 0.00601649, -0.01855508,  0.00017399],
       [ 0.00595795, -0.02294921,  0.02800064],
       [-0.02460223,  0.01526762,  0.01087747],
       [-0.00900424,  0.00726531,  0.01150711],
       [-0.01207218, -0.02462552,  0.01588084],
       [ 0.01408512,  0.02512257, -0.02771789],
       [ 0.02072795,  0.02207844, -0.00283951],
       [-0.00110884, -0.00360411, -0.03588265]])

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

   pos_ref = array([[  1.04158945e-02,  -9.87095573e-04,   4.63078398e-04],
       [  1.37069437e+00,   1.35732274e+00,   1.36231260e+00],
       [  2.73317987e+00,   2.72473320e+00,  -2.53175180e-03],
       [  4.09431672e+00,   4.07259172e+00,   1.34964664e+00],
       [  2.72208620e+00,  -1.54533813e-02,   2.70315737e+00],
       [  4.07948880e+00,   1.34575563e+00,   4.07878517e+00],
       [ -3.89972060e-03,   2.71244938e+00,   2.72906452e+00],
       [  1.35610951e+00,   4.08095355e+00,   4.06509952e+00]])

   velo_ref = array([[ 0.01003832, -0.00199037, -0.01122285],
       [-0.0043836 ,  0.00514716,  0.01096999],
       [-0.00706223, -0.00601213,  0.0024837 ],
       [ 0.01965012, -0.00356979, -0.00803952],
       [-0.0028455 ,  0.00181915, -0.00144307],
       [ 0.00515436, -0.00648117, -0.01035976],
       [ 0.00507576,  0.0088829 , -0.00847764],
       [ 0.00127142, -0.00010429,  0.00312452]])

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
