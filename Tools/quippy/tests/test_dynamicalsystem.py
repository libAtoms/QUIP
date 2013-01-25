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


if 'ase' in available_modules:
   import ase
else:
   import quippy.miniase as ase

class GenericTestDynamics(object):
   def common_init(self):
      random.seed(1)
      system_reseed_rng(1)
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

   pos_ref = array([[  4.96714153e-03,  -1.38264301e-03,   6.47688538e-03],
       [  1.37523030e+00,   1.35765847e+00,   1.35765863e+00],
       [  2.73579213e+00,   2.72767435e+00,  -4.69474386e-03],
       [  4.08542560e+00,   4.07536582e+00,   1.35534270e+00],
       [  2.72241962e+00,  -1.91328024e-02,   2.70275082e+00],
       [  4.07437712e+00,   1.34987169e+00,   4.08314247e+00],
       [ -9.08024076e-03,   2.70587696e+00,   2.73465649e+00],
       [  1.35774224e+00,   4.08067528e+00,   4.06575252e+00]])

   velo_ref = array([[  8.41198427e-05,   1.40681531e-05,  -2.37557008e-04],
       [ -2.56141803e-04,   2.21481172e-04,   5.19388261e-05],
       [ -4.75475980e-04,  -1.34614744e-04,  -1.13147122e-04],
       [  2.65109135e-04,  -1.57910324e-04,   1.53764691e-04],
       [  7.40295797e-05,   3.90376214e-04,   3.24470115e-04],
       [ -1.88182438e-05,  -2.50792632e-04,  -5.07363380e-05],
       [  2.71947337e-04,   6.95737223e-05,  -5.25858957e-04],
       [  5.52301328e-05,  -1.52181560e-04,   3.97125793e-04]])

   def setUp(self):
      self.common_init()
      
      # advance one step with the DynamicalSystem + quippy Atoms
      self.ds.atoms.set_cutoff(self.pot.cutoff() + self.skin)
      self.ds.atoms.calc_connect()
      self.ds.advance_verlet1(self.dt)

      self.pot.calc(self.ds.atoms, args_str='force')
      self.ds.advance_verlet2(self.dt, self.ds.atoms.force)

      # Advance one step with the Dynamics wrapper
      f0 = self.dyn.atoms.get_forces()
      f0 = zeros((len(self.at2), 3))
      self.f1 = self.dyn.step(f0)


class TestDynamics_50steps(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[  4.09299760e-03,  -6.69382577e-03,  -3.88404207e-03],
                    [  1.35823869e+00,   1.35084502e+00,   1.35533860e+00],
                    [  2.72007596e+00,   2.70610397e+00,   9.81540013e-04],
                    [  4.07979455e+00,   4.07587373e+00,   1.35527734e+00],
                    [  2.71975065e+00,  -1.59847278e-03,   2.72323237e+00],
                    [  4.09092971e+00,   1.36288737e+00,   4.07480938e+00],
                    [  9.89206367e-03,   2.71662593e+00,   2.71668154e+00],
                    [  1.36409929e+00,   4.07256340e+00,   4.07864905e+00]])

   velo_ref = array([[-0.0050459 , -0.0028246 ,  0.00078659],
                     [ 0.00145153, -0.01111762, -0.00010991],
                     [ 0.00866064, -0.0028097 ,  0.00399148],
                     [-0.01078164,  0.00422337, -0.00196776],
                     [-0.00125813, -0.00367343, -0.00274162],
                     [ 0.00634546,  0.01045345, -0.0005444 ],
                     [-0.00191543,  0.00526222,  0.00476076],
                     [ 0.00254349,  0.00048632, -0.00417513]])
   
   def setUp(self):
      self.common_init()
      
      # advance 50 steps with the DynamicalSystem + quippy Atoms
      self.ds.atoms.set_cutoff(self.pot.cutoff() + self.skin)
      self.ds.atoms.calc_connect()
      for i in range(50):
         self.ds.advance_verlet1(self.dt)
         self.pot.calc(self.ds.atoms, args_str='force')
         self.ds.advance_verlet2(self.dt, self.ds.atoms.force)
         if i % 10 == 0:
            self.ds.atoms.calc_connect()

      # Advance 50 steps with the Dynamics wrapper
      self.dyn.run(50)
      

class TestDynamics_InitTemperature(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[ 0.05828672, -0.04175836,  0.00744327],
                    [ 1.34867502,  1.40435697,  1.36197131],
                    [ 2.74389059,  2.72356746, -0.01471436],
                    [ 4.09077239,  4.06040704,  1.42727504],
                    [ 2.63838305, -0.01856504,  2.64833862],
                    [ 4.08272082,  1.31363493,  4.04177582],
                    [ 0.01987994,  2.75198344,  2.78664096],
                    [ 1.36426537,  4.08298068,  4.04235511]])
 
   velo_ref = array([[ 0.00612828, -0.01856402,  0.00025044],
                     [ 0.00586698, -0.02294197,  0.02786251],
                     [-0.02464898,  0.01508121,  0.01099344],
                     [-0.00902715,  0.00729644,  0.01134272],
                     [-0.01221629, -0.02461637,  0.01599742],
                     [ 0.01423147,  0.02536111, -0.02778432],
                     [ 0.02085133,  0.02196878, -0.00269975],
                     [-0.00118563, -0.00358519, -0.03596246]])
   
   def setUp(self):
      self.common_init()

      # choose some random velocities for one case
      self.dyn.set_temperature(300.)
      
      # use same random initial velocities
      self.ds.atoms.velo[...] = transpose(self.dyn.atoms.get_velocities()/sqrt(MASSCONVERT))

      # advance 50 steps with the DynamicalSystem + quippy Atoms
      self.ds.atoms.set_cutoff(self.pot.cutoff() + self.skin)
      self.ds.atoms.calc_connect()
      for i in range(50):
         self.ds.advance_verlet1(self.dt)
         self.pot.calc(self.ds.atoms, args_str='force')
         self.ds.advance_verlet2(self.dt, self.ds.atoms.force)
         if self.ds.nsteps % 10 == 0:
            self.ds.atoms.calc_connect()

      # Advance 50 steps with the Dynamics wrapper
      self.dyn.run(50)


class TestDynamics_Thermostat(GenericTestDynamics, QuippyTestCase):

   pos_ref = array([[  1.03441912e-02,  -9.93765959e-04,   6.61033474e-04],
       [  1.37091923e+00,   1.35714342e+00,   1.36227657e+00],
       [  2.73358216e+00,   2.72486376e+00,  -2.44797826e-03],
       [  4.09410455e+00,   4.07271864e+00,   1.34952688e+00],
       [  2.72203064e+00,  -1.57845266e-02,   2.70287199e+00],
       [  4.07948929e+00,   1.34594514e+00,   4.07883750e+00],
       [ -4.13997666e-03,   2.71239021e+00,   2.72950169e+00],
       [  1.35606148e+00,   4.08108287e+00,   4.06476956e+00]])

   velo_ref = array([[  9.98677395e-03,  -1.98364045e-03,  -1.10884055e-02],
       [ -4.20611633e-03,   5.03863782e-03,   1.09616036e-02],
       [ -6.77615246e-03,  -5.88272275e-03,   2.51711678e-03],
       [  1.95264142e-02,  -3.49415152e-03,  -8.10196389e-03],
       [ -2.86898483e-03,   1.58132087e-03,  -1.66973843e-03],
       [  5.12033588e-03,  -6.39646390e-03,  -1.03023615e-02],
       [  4.88230847e-03,   8.83891446e-03,  -8.18347793e-03],
       [  1.23393064e-03,  -1.04109955e-05,   2.90271098e-03]])
   
   def setUp(self):
      self.common_init()

      system_reseed_rng(1)
      self.ds.add_thermostat(THERMOSTAT_LANGEVIN, 300.0, tau=500.0)
      
      # advance 10 steps with the DynamicalSystem + quippy Atoms
      self.ds.atoms.set_cutoff(self.pot.cutoff() + self.skin)
      self.ds.atoms.calc_connect()
      for i in range(10):
         self.ds.advance_verlet1(self.dt)
         self.pot.calc(self.ds.atoms, args_str='force')
         self.ds.advance_verlet2(self.dt, self.ds.atoms.force)
         if self.ds.nsteps % 10 == 0:
            self.ds.atoms.calc_connect()

      # we need same sequence of thermostat random forces
      system_reseed_rng(1)
      self.dyn.add_thermostat(THERMOSTAT_LANGEVIN, 300.0, tau=500.0)

      # Advance 10 steps with the Dynamics wrapper
      self.dyn.run(10)


if __name__ == '__main__':
   unittest.main()
