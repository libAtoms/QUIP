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

if 'ase' in available_modules:
   import ase
else:
   import quippy.miniase as ase

class TestConvert_Simple(QuippyTestCase):

   def setUp(self):
      self.dia_quippy   = diamond(5.44, 14)
      self.dia_ase      = ase.Atoms(self.dia_quippy)
      self.dia_quippy_2 = quippy.Atoms(self.dia_ase)

   def test_equal(self):
      self.assertEqual(self.dia_quippy, self.dia_quippy_2)

class TestConvert_ExtraProperties(QuippyTestCase):

   def setUp(self):
      self.dia_quippy   = diamond(5.44, 14)
      self.dia_quippy.add_property('magmoms', np.random.uniform(-1,1,size=3*8).reshape(3,8))
      self.dia_quippy.add_property('charge', [1.0]*8)
      self.dia_ase      = ase.Atoms(self.dia_quippy)
      self.dia_quippy_2 = quippy.Atoms(self.dia_ase)

   def test_equal(self):
      self.assertEqual(self.dia_quippy, self.dia_quippy_2)

class TestConvert_UnsupportedProperties(QuippyTestCase):

   def setUp(self):
      self.dia_quippy   = diamond(5.44, 14)
      self.dia_quippy.add_property('new', -1.0)
      self.dia_ase      = ase.Atoms(self.dia_quippy)
      self.dia_quippy_2 = quippy.Atoms(self.dia_ase)

   def test_not_equal(self):
      self.assertNotEqual(self.dia_quippy, self.dia_quippy_2)

   def test_equal_after_remove(self):
      self.dia_quippy.remove_property('new')
      self.assertEqual(self.dia_quippy, self.dia_quippy_2)


class TestCalculator_SW_Potential(QuippyTestCase):
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

      system_reseed_rng(2065775975)
      self.pot = Potential('IP SW', param_str=xml)

      self.at = diamond(5.44, 14)
      randomise(self.at.pos, 0.1)
      self.at.set_cutoff(self.pot.cutoff())
      self.at.calc_connect()

      self.e, = fvar('e')
      self.f = fzeros((3,self.at.n))
      self.df = fzeros((3,self.at.n))
      self.v = fzeros((3,3))

      self.e_ref = -34.5038375509

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

      self.ase_at = ase.Atoms(self.at)
      self.ase_at.set_calculator(self.pot)

   def test_energy(self):
      self.assertAlmostEqual(self.ase_at.get_potential_energy(), self.e_ref)

   def test_forces(self):
      self.assertArrayAlmostEqual(self.ase_at.get_forces(), self.f_ref.T)

   def test_stress(self):
      from quippy.elasticity import stress_vector
      self.assertArrayAlmostEqual(self.ase_at.get_stress(), stress_vector(-self.v_ref/self.at.cell_volume()))

   def test_numeric_forces(self):
      self.assertArrayAlmostEqual(self.pot.get_numeric_forces(self.ase_at), self.f_ref.T, tol=1e-4)


class TestNumbersSpecies(QuippyTestCase):
    def setUp(self):
        self.ase_at = ase.Atoms(numbers=[1,2,3], positions=np.zeros((3,3)))
        self.at = Atoms(numbers=[1,2,3], positions=np.zeros((3,3)))

    def test_get_atomic_numbers_ase(self):
        self.assertEqual(list(self.ase_at.get_atomic_numbers()), [1,2,3])

    def test_get_atomic_numbers_quippy(self):
        self.assertEqual(list(self.at.get_atomic_numbers()), [1,2,3])

    def test_get_chemical_symbols_ase(self):
        self.assertEqual(self.ase_at.get_chemical_symbols(), ['H', 'He', 'Li'])

    def test_get_chemical_symbols_quippy(self):
        self.assertEqual(self.at.get_chemical_symbols(), ['H', 'He', 'Li'])

    def test_get_species(self):
        sp = list([''.join(s).strip() for s in self.at.species])
        self.assertEqual(sp, ['H', 'He', 'Li'])

    def test_set_atomic_numbers(self):
        self.at.set_atomic_numbers([4,5,6])
        self.assertEqual(list(self.at.z), [4,5,6])
        self.assertEqual(list(self.at.get_atomic_numbers()), [4,5,6])
        self.assertEqual(list(self.at.get_chemical_symbols()), ['Be', 'B', 'C'])
        sp = list([''.join(s).strip() for s in self.at.species])        
        self.assertEqual(sp, ['Be', 'B', 'C'])

    def test_set_chemical_symbols(self):
        self.at.set_chemical_symbols(['Be', 'B', 'C'])
        self.assertEqual(list(self.at.z), [4,5,6])        
        self.assertEqual(list(self.at.get_atomic_numbers()), [4,5,6])
        self.assertEqual(list(self.at.get_chemical_symbols()), ['Be', 'B', 'C'])
        sp = list([''.join(s).strip() for s in self.at.species])        
        self.assertEqual(sp, ['Be', 'B', 'C'])

    def test_set_atomic_numbers_set_species_false(self):
        self.at.set_atomic_numbers([4,5,6], set_species=False)
        self.assertEqual(list(self.at.z), [4,5,6])
        self.assertEqual(list(self.at.get_atomic_numbers()), [4,5,6])
        self.assertEqual(list(self.at.get_chemical_symbols()), ['Be', 'B', 'C'])
        sp = list([''.join(s).strip() for s in self.at.species])        
        self.assertEqual(sp, ['H', 'He', 'Li']) # still previous value

    def test_set_chemical_symbols_set_species_false(self):
        self.at.set_chemical_symbols(['Be', 'B', 'C'], set_species=False)
        self.assertEqual(list(self.at.z), [4,5,6])        
        self.assertEqual(list(self.at.get_atomic_numbers()), [4,5,6])
        self.assertEqual(list(self.at.get_chemical_symbols()), ['Be', 'B', 'C'])
        sp = list([''.join(s).strip() for s in self.at.species])
        self.assertEqual(sp, ['H', 'He', 'Li']) # still previous value

    def test_constructor_set_species_false(self):
        at = Atoms(self.ase_at, set_species=False)
        self.assert_(not at.has_property('species'))
        

class TestOperators(QuippyTestCase):

    def setUp(self):

        self.a = Atoms("""2
Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3:Z:I:1 nneightol=1.2 cutoff=0.0
Si       0.00000000       0.00000000       0.00000000       14 
Si       2.00000000       2.00000000       2.00000000       14""", format='string')
        self.b = Atoms("""2
Lattice="20.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3:Z:I:1 nneightol=1.2 cutoff=0.0
C       1.00000000       0.00000000       0.00000000        6 
C       3.00000000       3.00000000       3.00000000        6""", format='string')

    def test_operator_plus(self):
        at = self.a + self.b
        self.assertEquals(list(at.z), list(self.a.z) + list(self.b.z))
        self.assertArrayAlmostEqual(at.get_positions(), np.r_[self.a.get_positions(),
                                                              self.b.get_positions()])

    def test_operator_plus_equals(self):
        a0 = self.a.copy()        
        self.a += self.b
        self.assertEquals(list(self.a.z), list(a0.z) + list(self.b.z))
        self.assertArrayAlmostEqual(self.a.get_positions(), np.r_[a0.get_positions(),
                                                                  self.b.get_positions()])
        

    def test_operator_times(self):
        at = self.a * (3, 3, 3)
        self.assertEquals(list(at.z), list(self.a.z)*(3*3*3))

    def test_operator_times_equals(self):
        a0 = self.a.copy()
        self.a *= (3, 3, 3)
        self.assertEquals(list(self.a.z), list(a0.z)*(3*3*3))


if __name__ == '__main__':
   unittest.main()
