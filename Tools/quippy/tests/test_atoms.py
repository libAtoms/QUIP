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
from quippy.atoms import *
from numpy import *

import unittest, quippy
from quippytest import *


class TestAtoms_LowLevel(QuippyTestCase):
   def setUp(self):
      self.properties = Dictionary()

      self.properties['pos'] = fzeros((3,10))
      self.properties['Z'] = [14]*10
      self.properties['species'] = s2a(['Si']*10, 10)
      
      self.at = Atoms(n=10, lattice=10.0*fidentity(3), properties=self.properties)
      self.at_def = Atoms(n=10, lattice=10.0*fidentity(3))

      self.dia = diamond(5.44, 14)
      self.dia.pos[:] += numpy.random.uniform(-0.1,0.1,self.dia.n*3).reshape(3,self.dia.n)
      self.dia.set_cutoff_factor(1.3)
      self.dia.calc_connect()

   def test_n(self):
      self.assertEqual(self.at.n, 10)

   def test_lattice(self):
      self.assertArrayAlmostEqual(self.at.lattice, [[10.0, 0.0, 0.0],
                                                    [0.0, 10.0, 0.0],
                                                    [0.0, 0.0, 10.0]])

   def test_g(self):
      self.assertArrayAlmostEqual(self.at.g, [[0.1, 0.0, 0.0],
                                              [0.0, 0.1, 0.0],
                                              [0.0, 0.0, 0.1]])

   def test_properties(self):
      self.assertEqual(self.at.properties, self.properties)

   def test_params(self):
      params = {'k1': 1, 'k2': 2.0, 'k3': 's', 'k4': fidentity(3)}
      self.at.params.update(params)

      for k in self.at.params:
         if isinstance(params[k], FortranArray):
            self.assertArrayAlmostEqual(params[k], self.at.params[k])
         else:
            self.assertEqual(params[k], self.at.params[k])

   def test_z(self):
      self.assert_((self.at.z == 14).all())

   def test_pos(self):
      self.assertAlmostEqual(abs(self.at.pos - 0.0).max(), 0.0)

   def test_species(self):
      self.assertEqual(a2s(self.at.species), ['Si']*self.at.n)

   def test_zero_1(self):
      self.at.pos[...] = 0.0
      self.at.z[...] = 0
      self.at.species[...] = ' '
      self.assertEqual(self.at, self.at_def)

   def test_zero_2(self):
      self.at.zero()
      self.assertEqual(self.at, self.at_def)

   def test_assignment(self):
      at2 = Atoms()
      at2.assignment(self.at)
      self.assertEqual(at2, self.at)
      
   def test_copy_without_connect(self):
      at2 = Atoms()
      self.at.calc_connect()
      self.assert_(self.at.connect.initialised)
      at2.copy_without_connect(self.at)
      self.assert_(not at2.connect.initialised)
      self.assertEqual(at2, self.at)

   def test_remove_property_1(self):
      self.at.remove_property('pos')
      self.assertRaises(KeyError, self.at.properties.__getitem__, 'pos')

   def test_remove_property_2(self):
      self.at.pos # access pos before removing to put it in cache
      self.at.remove_property('pos')
      self.assertRaises(AttributeError, getattr, self.at, 'pos')

   def test_has_property_1(self):
      self.assert_(all([self.at.has_property(p) for p in self.at.properties.keys()]))

   def test_has_property_2(self):
      self.assert_(not self.at.has_property('pos_'))

   def test_set_cutoff(self):
      self.at.set_cutoff(5.0, 6.0)
      self.assertEqual((self.at.cutoff, self.at.cutoff_break, self.at.use_uniform_cutoff), (5.0, 6.0, True))
      
   def test_set_cutoff_factor(self):
      self.at.set_cutoff_factor(1.3, 1.4)
      self.assertEqual((self.at.cutoff, self.at.cutoff_break, self.at.use_uniform_cutoff), (1.3, 1.4, False))

   def test_set_cutoff_factor_default(self):
      self.at.set_cutoff_factor(0.0, 0.0)
      self.assertEqual((self.at.cutoff, self.at.cutoff_break, self.at.use_uniform_cutoff), (DEFAULT_NNEIGHTOL, DEFAULT_NNEIGHTOL, False))

   def test_set_cutoff_minimum_1(self):
      self.at.set_cutoff(5.0, 6.0)
      self.at.set_cutoff_minimum(6.0, 7.0)
      self.assertEqual((self.at.cutoff, self.at.cutoff_break, self.at.use_uniform_cutoff), (6.0, 7.0, True))      

   def test_set_cutoff_minimum_2(self):
      self.at.set_cutoff(5.0, 6.0)
      self.at.set_cutoff_minimum(4.0, 5.0)
      self.assertEqual((self.at.cutoff, self.at.cutoff_break, self.at.use_uniform_cutoff), (5.0, 6.0, True))      

   def test_cutoff_1(self):
      self.at.set_cutoff(5.0)
      self.assertEqual(self.at.cutoff_(14, 14), 5.0)

   def test_cutoff_2(self):
      self.at.set_cutoff_factor(1.3)
      self.assertEqual(self.at.cutoff_(14, 14), 1.3*bond_length(14, 14))
      
   def test_cutoff_break_1(self):
      self.at.set_cutoff(5.0, 6.0)
      self.assertEqual(self.at.cutoff_break_(14, 14), 6.0)

   def test_cutoff_break_2(self):
      self.at.set_cutoff_factor(1.3, 1.4)
      self.assertEqual(self.at.cutoff_break_(14, 14), 1.4*bond_length(14, 14))

   def test_set_lattice_1(self):
      L = 5.0*fidentity(3)
      self.at.set_lattice(L, scale_positions=False)
      self.assertArrayAlmostEqual(self.at.lattice, L)
      
   def test_set_lattice_2(self):
      L = 5.0*fidentity(3)
      g = linalg.inv(L)
      self.at.set_lattice(L, scale_positions=False)
      self.assertArrayAlmostEqual(self.at.g, g)

   def test_set_lattice_scale_positions_1(self):
      L = 5.0*fidentity(3)
      p0 = fzeros((3,self.dia.n))
      p0[:] = self.dia.pos[:]
      self.dia.set_lattice(L, scale_positions=True)
      self.assertArrayAlmostEqual(self.dia.lattice, L)
      self.assert_(all(abs(p0/self.dia.pos - 5.44/5.0) < 1e-7))

   def test_set_lattice_scale_positions_2(self):
      L = 5.0*fidentity(3)
      p0 = fzeros((3,self.dia.n))
      p0[:] = self.dia.pos[:]
      self.dia.set_lattice(L, scale_positions=False)
      self.assertArrayAlmostEqual(self.dia.lattice, L)
      self.assert_(all(abs(p0/self.dia.pos - 1.0) < 1e-7))      

   def test_set_atoms_1(self):
      self.at.set_atoms(6)
      self.assertEqual(list(self.at.z), [6]*self.at.n)

   def test_set_atoms_2(self):
      self.at.set_atoms(6)
      self.assertEqual(a2s(self.at.species), ['C']*self.at.n)

   def test_set_atoms_3(self):
      new_z = [6]*5 + [14]*5
      self.at.set_atoms(new_z)
      self.assertEqual(list(self.at.z), new_z)

   def test_set_atoms_4(self):
      new_z = [6]*5 + [14]*5
      new_species = [ElementName[z] for z in new_z]
      self.at.set_atoms(new_z)
      self.assertEqual(a2s(self.at.species), new_species)

   def test_add_atom_1(self):
      self.at.add_atoms([1.0, 1.0, 1.0], 6)
      self.assertEqual(list(self.at.z), [14]*10 + [6])

   def test_add_atom_2(self):
      self.at.add_atoms([1.0, 1.0, 1.0], 6)
      self.assertArrayAlmostEqual(self.at.pos[-1], [1.0, 1.0, 1.0])

   def test_add_atom_3(self):
      self.at.add_atoms(farray([[1.0, 1.0, 1.0],
                                [2.0, 2.0, 2.0]]).T, [6, 6])
      self.assertEqual(list(self.at.z), [14]*10 + [6]*2)
      
   def test_add_atom_4(self):
      self.at.add_atoms(farray([[1.0, 1.0, 1.0],
                                [2.0, 2.0, 2.0]]).T, [6, 6])
      self.assertArrayAlmostEqual(self.at.pos[:,-2:].T, [[1.0, 1.0, 1.0],
                                                         [2.0, 2.0, 2.0]])
   def test_add_atom_5(self):
      self.at.add_atoms(farray([[1.0, 1.0, 1.0],
                                [2.0, 2.0, 2.0]]).T, [6, 6], mass=[10.0, 10.0])
      self.assertArrayAlmostEqual(self.at.pos[:,-2:].T, [[1.0, 1.0, 1.0],
                                                         [2.0, 2.0, 2.0]])
   def test_add_atom_6(self):
      self.at.add_atoms(farray([[1.0, 1.0, 1.0],
                                [2.0, 2.0, 2.0]]).T, [6, 6], mass=[10.0, 10.0])
      self.assertArrayAlmostEqual([self.at.z[-2], self.at.z[-1]], [6, 6])

   def test_add_atom_7(self):
      self.at.add_atoms(farray([[1.0, 1.0, 1.0],
                                [2.0, 2.0, 2.0]]).T, [6, 6], mass=[10.0, 10.0])
      self.assertArrayAlmostEqual([self.at.mass[-2], self.at.mass[-1]], [10.0, 10.0])

   def test_remove_atoms_1(self):
      cp = self.dia.copy()
      cp.remove_atoms(1)
      self.assertEqual(cp[1], self.dia[8])

   def test_remove_atoms_2(self):
      cp1 = self.dia.copy()
      cp1.remove_atoms(1)
      cp2 = self.dia.select(list=[8,2,3,4,5,6,7],orig_index=False)
      self.assertEqual(cp1, cp2)

   def test_remove_atoms_3(self):
      cp1 = self.dia.copy()
      cp1.remove_atoms([1,5])
      cp2 = self.dia.select(list=[8,2,3,4,7,6],orig_index=False)
      self.assertEqual(cp1, cp2)

   def test_add_atoms_fixed_size(self):
      a = Atoms(n=10,fixed_size=True)
      self.assertRaises(RuntimeError, a.add_atoms, [1.0, 2.0, 3.0], 6)

   def test_remove_atoms_fixed_size(self):
      a = Atoms(n=10,fixed_size=True)
      self.assertRaises(RuntimeError, a.remove_atoms, 10)

   def test_remove_atoms_cache_error(self):
      cp = self.dia.copy()
      self.assertEqual(cp.pos.shape, (3, 8))
      cp.remove_atoms(1)
      self.assertEqual(cp.pos.shape, (3, 7))

   def test_map_into_cell_range(self):
      self.dia.map_into_cell()
      t = dot(self.dia.g, self.dia.pos)
      self.assert_(t.max() < 0.5 and t.min() >= -0.5)

   def test_map_into_cell_shift(self):
      self.dia.map_into_cell()
      t1 = dot(self.dia.g, self.dia.pos)
      self.dia.pos[:] = dot(self.dia.lattice, t1 + 1.0)
      self.dia.map_into_cell()
      t2 = dot(self.dia.g, self.dia.pos)
      self.assertArrayAlmostEqual(t1, t2)
      
   def test_map_into_cell_idempotence(self):
      cp1 = self.dia.copy()
      cp1.map_into_cell()
      cp2 = cp1.copy()
      cp2.map_into_cell()
      self.assertEqual(cp1, cp2)

   def test_calc_dists(self):
      d12_1 = fzeros(3); d12_2 = fzeros(3)
      delta = farray([0.1, 0.1, 0.1])
      self.dia.neighbour(1, 1, diff=d12_1)
      self.dia.pos[:,1] += delta
      self.dia.calc_dists()
      self.dia.neighbour(1, 1, diff=d12_2)
      self.assertAlmostEqual((d12_1 - d12_2).norm(), delta.norm())
      
   def test_diff(self):
      d1 = fzeros(3)
      shift = fzeros(3, dtype=numpy.int32)
      j = self.dia.neighbour(1, 1, diff=d1, shift=shift)
      d2 = self.dia.diff(1, j, shift)
      self.assertArrayAlmostEqual(d1, d2)

   def test_diff_min_image(self):
      d1 = fzeros(3)
      j = self.dia.neighbour(1, 1, diff=d1)
      d2 = self.dia.diff_min_image(1, j)
      self.assertArrayAlmostEqual(d1, d2)

   def test_diff_min_image_atom_vec(self):
      d1 = fzeros(3)
      j = self.dia.neighbour(1, 1, diff=d1)
      d2 = self.dia.diff_min_image(1, self.dia.pos[j])
      self.assertArrayAlmostEqual(d1, d2)

   def test_diff_min_image_vec_atom(self):
      d1 = fzeros(3)
      j = self.dia.neighbour(1, 1, diff=d1)
      d2 = self.dia.diff_min_image(self.dia.pos[1], j)
      self.assertArrayAlmostEqual(d1, d2)
   
   def test_diff_min_image_vec_vec(self):
      d1 = fzeros(3)
      j = self.dia.neighbour(1, 1, diff=d1)
      d2 = self.dia.diff_min_image(self.dia.pos[1], self.dia.pos[j])
      self.assertArrayAlmostEqual(d1, d2)
   
   def test_distance(self):
      d1 = farray(0.0)
      shift = fzeros(3, dtype=numpy.int32)
      j = self.dia.neighbour(1, 1, distance=d1, shift=shift)
      d2 = self.dia.distance(1, j, shift)
      self.assertAlmostEqual(d1, d2)
      
   def test_distance_min_image(self):
      d1 = farray(0.0)
      j = self.dia.neighbour(1, 1, distance=d1)
      d2 = self.dia.distance_min_image(1, j)
      self.assertAlmostEqual(d1, d2)

   def test_distance_min_image_atom_vec(self):
      d1 = farray(0.0)
      j = self.dia.neighbour(1, 1, distance=d1)
      d2 = self.dia.distance_min_image(1, self.dia.pos[j])
      self.assertAlmostEqual(d1, d2)

   def test_distance_min_image_vec_atom(self):
      d1 = farray(0.0)
      j = self.dia.neighbour(1, 1, distance=d1)
      d2 = self.dia.distance_min_image(self.dia.pos[1], j)
      self.assertAlmostEqual(d1, d2)
   
   def test_distance_min_image_vec_vec(self):
      d1 = farray(0.0)
      j = self.dia.neighbour(1, 1, distance=d1)
      d2 = self.dia.distance_min_image(self.dia.pos[1], self.dia.pos[j])
      self.assertAlmostEqual(d1, d2)

   def test_realpos(self):
      t_dia = self.dia.copy()
      t_dia.map_into_cell()
      rp1 = t_dia.realpos(1)
      t_dia.travel[:] = 1
      rp2 = t_dia.realpos(1)
      self.assertArrayAlmostEqual(rp2 - rp1, dot(t_dia.lattice, [1,1,1]))

   def test_cosine(self):
      c1 = fzeros(3)
      c2 = fzeros(3)
      j = self.dia.neighbour(1, 1, cosines=c1)
      k = self.dia.neighbour(1, 2, cosines=c2)
      c3 = self.dia.cosine(1, j, k)
      self.assertArrayAlmostEqual(dot(c1, c2), c3)

   def test_cosine_neighbour(self):
      c1 = fzeros(3)
      c2 = fzeros(3)
      self.dia.neighbour(1, 1, cosines=c1)
      self.dia.neighbour(1, 2, cosines=c2)
      c3 = self.dia.cosine_neighbour(1, 1, 2)
      self.assertArrayAlmostEqual(dot(c1, c2), c3)
      
   def test_direction_cosines(self):
      c1 = fzeros(3)
      shift = fzeros(3, dtype=int32)
      j = self.dia.neighbour(1, 1, cosines=c1, shift=shift)
      c2 = self.dia.direction_cosines(1, j, shift=shift)
      self.assertArrayAlmostEqual(c1, c2)
      
   def test_direction_cosines(self):
      c1 = fzeros(3)
      j = self.dia.neighbour(1, 1, cosines=c1)
      c2 = self.dia.direction_cosines_min_image(1, j)
      self.assertArrayAlmostEqual(c1, c2)

   def test_prop_names_string(self):
      self.assertEqual(self.at.prop_names_string().strip(), ':'.join(self.at.properties.keys()))

   def test_bond_length(self):
      self.assertAlmostEqual(bond_length(14,14), 2*ElementCovRad[14])

   def test_termination_bond_rescale(self):
      self.assertAlmostEqual(termination_bond_rescale(14, 14), (ElementCovRad[14]+ElementCovRad[1])/(2*ElementCovRad[14]))

   def test_complement(self):
      L1 = Table(1,0,0,0)
      L1.append([1,2,3])
      L2 = self.at.complement(L1)
      self.assertEqual(list(L2.int[1,:]), [i for i in frange(self.at.n) if i not in L1.int[1,:]])

   def test_difference(self):
      L1 = Table(1,0,0,0)
      L1.append([1,2,3,4,5,6,7])
      L2 = Table(1,0,0,0)
      L2.append([2,4,6])
      L3 = difference(L1, L2)
      self.assertEqual(list(L3.int[1,:]), list(set(L1.int[1,:]).difference(L2.int[1,:])))


class TestAtoms_Rotate(QuippyTestCase):

   def setUp(self):
      self.a = diamond(5.44, 14)
      self.a.add_property('force', 0., n_cols=3)
      self.a.force[:]  = numpy.random.uniform(-0.1,0.1,self.a.n*3).reshape(3,self.a.n)
      self.a.add_property('not_a_vector', 0., n_cols=3)
      self.a.not_a_vector[:]  = numpy.random.uniform(-0.1,0.1,self.a.n*3).reshape(3,self.a.n)
      self.a.params['virial'] = numpy.random.uniform(-0.1,0.1,9).reshape(3,3)
      self.c = self.a.copy()

      theta = pi/3.
      self.c.rotate([1., 0, 0], theta)
      self.R = farray([[1.,   0.,          0.,        ],
                       [0.,   cos(theta), -sin(theta) ],
                       [0.,   sin(theta),  cos(theta) ]])

   def test_lattice(self):
      # Lattice transformed
      self.assertArrayAlmostEqual(dot(self.R, self.a.lattice), self.c.lattice)

   def test_pos(self):
      # Vector properties transformed
      self.assertArrayAlmostEqual(dot(self.R, self.a.pos), self.c.pos)

   def test_force(self):
      self.assertArrayAlmostEqual(dot(self.R, self.a.force), self.c.force)

   def test_non_vector(self):
      # Non-vector properties shouldn't be changed
      self.assertArrayAlmostEqual(self.a.not_a_vector, self.c.not_a_vector)

   def test_frac(self):
      # Fractional coordinates should match
      self.assertArrayAlmostEqual(dot(self.a.g, self.a.pos), dot(self.c.g, self.c.pos))

   def test_virial(self):
      # Trace should be preserved
      self.assertAlmostEqual(self.a.virial.trace(), self.c.virial.trace())
   
   

class TestAtoms_Extras(QuippyTestCase):
   def setUp(self):
      self.a = 5.44
      self.at = diamond(self.a, 14)
      self.params = {'k1': 1, 'k2': 2.0, 'k3': 's', 'k4': fidentity(3)}
      self.at.params.update(self.params)

   def test_copy(self):
      cp = self.at.copy()
      self.assertEqual(cp, self.at)
      self.assertNotEqual(id(cp), id(self.at))
      self.assert_(not all(cp._fpointer == self.at._fpointer))

   def test_copy_constructor(self):
      cp = Atoms(self.at)
      self.assertEqual(cp, self.at)
      self.assertNotEqual(id(cp), id(self.at))
      self.assert_(not all(cp._fpointer == self.at._fpointer))

   def test_len(self):
      self.assertEqual(len(self.at), self.at.n)

   def test_getattr(self):
      for k in self.at.params:
         if isinstance(self.params[k], FortranArray):
            self.assertArrayAlmostEqual(self.params[k], getattr(self.at, k))
         else:
            self.assertEqual(self.params[k], getattr(self.at, k))

   def test_getitem(self):
      d = self.at[1]
      self.assertEqual(d, {'Z': 14, 'pos': [0.0, 0.0, 0.0], 'species': 'Si'})

   def test_eq1(self):
      cp = self.at.copy()
      cp.add_atoms([0., 0., 0.], 14)
      self.assertNotEqual(cp, self.at)

   def test_eq2(self):
      d = diamond(self.a, 14)
      d.params.update(self.at.params)
      self.assertEqual(d, self.at)

   def test_eq3(self):
      cp = self.at.copy()
      cp.params['f'] = 1
      self.assertNotEqual(cp, self.at)

   def test_eq4(self):
      cp = self.at.copy()
      cp.add_property('f', 0)
      self.assertNotEqual(cp, self.at)

   def test_select_noargs(self):
      self.assertRaises(ValueError, self.at.select)

   def test_select_mask_all(self):
      sub_at = self.at.select(mask=[1]*self.at.n, orig_index=False)
      self.assertEqual(self.at, sub_at)

   def test_select_mask_one(self):
      sub_at = self.at.select(mask=[1,0,0,0,0,0,0,0], orig_index=False)
      self.assertEqual(sub_at[1], self.at[1])

   def test_select_list(self):
      sub_at = self.at.select(list=[2], orig_index=False)
      self.assertEqual(sub_at[1], self.at[2])

   def test_select_list_array(self):
      sub_at = self.at.select(list=array([2]), orig_index=False)
      self.assertEqual(sub_at[1], self.at[2])

   def test_select_list_farray(self):
      sub_at = self.at.select(list=farray([2]), orig_index=False)
      self.assertEqual(sub_at[1], self.at[2])

   def test_select_list_2(self):
      sub_at = self.at.select(list=[2,3], orig_index=False)
      self.assertEqual((sub_at[1],  sub_at[2]),
                       (self.at[2], self.at[3]))

   def check_property(self, name, type, n_cols=None):
      t, l, (l1,l2) = self.at.properties.get_type_and_size(name)
      if type in (T_INTEGER_A, T_REAL_A, T_LOGICAL_A):
         s = self.at.n
         s1 = 0
         s2 = 0
      elif type in (T_INTEGER_A2, T_REAL_A2):
         s = 0
         s1 = n_cols
         s2 = self.at.n
      elif type == T_CHAR_A:
         s = 0
         s1 = TABLE_STRING_LENGTH
         s2 = self.at.n
      else:
         raise TypeError('bad property type %d' % type)
      self.assertEqual((t, l, l1, l2), (type, s, s1, s2))
      
   def test_add_property_int_from_scalar(self):
      self.at.add_property('int', 0)
      self.check_property('int', T_INTEGER_A)
      self.assertEqual(list(self.at.int), [0]*8)

   def test_add_property_int_from_list(self):
      self.at.add_property('int', [0]*8)
      self.check_property('int', T_INTEGER_A)
      self.assertEqual(list(self.at.int), [0]*8)

   def test_add_property_int_from_list_bad_length(self):
      p0 = self.at.properties.copy()
      self.assertRaises(RuntimeError, self.at.add_property, 'int', [0]*9)
      self.assertEqual(self.at.properties, p0)

   def test_add_property_int_from_array(self):
      self.at.add_property('int', fzeros(8,dtype=int))
      self.check_property('int', T_INTEGER_A)
      self.assertEqual(list(self.at.int), [0]*8)      
      
   def test_add_property_int2d_from_scalar(self):
      self.at.add_property('int2d', 0, n_cols=2)
      self.check_property('int2d', T_INTEGER_A2, 2)
      self.assertArrayAlmostEqual(self.at.int2d, [[0]*8]*2)

   def test_add_property_int2d_from_list(self):
      self.at.add_property('int2d', [[0]*8]*2)
      self.check_property('int2d', T_INTEGER_A2, 2)
      self.assertArrayAlmostEqual(self.at.int2d, [[0]*8]*2)
   
   def test_add_property_int2d_from_list_bad_length(self):
      p0 = self.at.properties.copy()
      self.assertRaises(RuntimeError, self.at.add_property, 'int2d', [[0]*9]*2)
      self.assertEqual(self.at.properties, p0)

   def test_add_property_int2d_from_array(self):
      self.at.add_property('int2d', fzeros((2,8),dtype=int))
      self.check_property('int2d', T_INTEGER_A2, 2)
      self.assertArrayAlmostEqual(list(self.at.int2d), [[0,0]]*8)

   def test_add_property_int_from_random_array(self):
      r = numpy.random.randint(1, 100, 8)
      self.at.add_property('int', r)
      self.check_property('int', T_INTEGER_A)
      self.assertArrayAlmostEqual(list(self.at.int), r)

   def test_add_property_real_from_scalar(self):
      self.at.add_property('real', 0.0)
      self.check_property('real', T_REAL_A)
      self.assertEqual(list(self.at.real), [0.0]*8)

   def test_add_property_real_from_list(self):
      self.at.add_property('real', [0.0]*8)
      self.check_property('real', T_REAL_A)
      self.assertEqual(list(self.at.real), [0.0]*8)

   def test_add_property_real_from_list_bad_length(self):
      p0 = self.at.properties.copy()
      self.assertRaises(RuntimeError, self.at.add_property, 'real', [0]*9)
      self.assertEqual(self.at.properties, p0)

   def test_add_property_real_from_array(self):
      self.at.add_property('real', fzeros(8))
      self.check_property('real', T_REAL_A)
      self.assertEqual(list(self.at.real), [0.]*8)      
      
   def test_add_property_real2d_from_scalar(self):
      self.at.add_property('real2d', 0.0, n_cols=2)
      self.check_property('real2d', T_REAL_A2, 2)
      self.assertArrayAlmostEqual(self.at.real2d, [[0.]*8]*2)

   def test_add_property_real2d_from_list(self):
      self.at.add_property('real2d', [[0.]*8]*2)
      self.check_property('real2d', T_REAL_A2, 2)
      self.assertArrayAlmostEqual(self.at.real2d, [[0.]*8]*2)
   
   def test_add_property_real2d_from_list_bad_length(self):
      p0 = self.at.properties.copy()
      self.assertRaises(RuntimeError, self.at.add_property, 'real2d', [[0,0]]*9)
      self.assertEqual(self.at.properties, p0)

   def test_add_property_real2d_from_array(self):
      self.at.add_property('real2d', fzeros((2,8)))
      self.check_property('real2d', T_REAL_A2, 2)
      self.assertArrayAlmostEqual(list(self.at.real2d), [[0.,0.]]*8)

   def test_add_property_real_from_random_array(self):
      r = numpy.random.uniform(1.0, 100.0, 8)
      self.at.add_property('real', r)
      self.check_property('real', T_REAL_A)      
      self.assertArrayAlmostEqual(list(self.at.real), r)
   
   def test_add_property_logical_from_scalar(self):
      self.at.add_property('logical', True)
      self.check_property('logical', T_LOGICAL_A)      
      self.assertEqual(list(self.at.logical), [True]*8)

   def test_add_property_logical_from_list(self):
      self.at.add_property('logical', [True]*8)
      self.check_property('logical', T_LOGICAL_A)      
      self.assertEqual(list(self.at.logical), [True]*8)

   def test_add_property_logical_from_list_bad_length(self):
      p0 = self.at.properties.copy()
      self.assertRaises(RuntimeError, self.at.add_property, 'logical', [True]*9)
      self.assertEqual(self.at.properties, p0)

   def test_add_property_logical_from_array(self):
      self.at.add_property('logical', fzeros(8,dtype=bool))
      self.check_property('logical', T_LOGICAL_A)            
      self.assertEqual(list(self.at.logical), [False]*8)      
      
   def test_add_property_logical_from_random_array(self):
      r = farray(numpy.random.randint(0, 2, 8), dtype=bool)
      self.at.add_property('logical', r)
      self.check_property('logical', T_LOGICAL_A)            
      self.assertArrayAlmostEqual(list(self.at.logical), r)

   def test_add_property_str_from_scalar(self):
      self.at.add_property('str', ' '*TABLE_STRING_LENGTH)
      self.check_property('str', T_CHAR_A)            
      self.assertEqual(self.at.str.shape, (TABLE_STRING_LENGTH, 8))

   def test_add_property_str_from_list(self):
      self.at.add_property('str', [' '*TABLE_STRING_LENGTH]*8)
      self.check_property('str', T_CHAR_A)      
      self.assertEqual(self.at.str.shape, (TABLE_STRING_LENGTH, 8))

   def test_add_property_str_from_list_bad_length(self):
      p0 = self.at.properties.copy()
      self.assertRaises(RuntimeError, self.at.add_property, 'str', ['']*9)
      self.assertEqual(self.at.properties, p0)

   def test_add_property_str_from_array(self):
       # array tranpose required because f2py tranposes character(len != 1) arrays
      self.at.add_property('str', fzeros((TABLE_STRING_LENGTH,8),dtype=str).T) 
      self.check_property('str', T_CHAR_A)      
      self.assertEqual(self.at.str.shape, (TABLE_STRING_LENGTH, 8))

   def test_add_property_property_type(self):
      self.at.add_property('logical', [1]*self.at.n, property_type=T_LOGICAL_A)
      self.check_property('logical', T_LOGICAL_A)
      self.assertEqual(list(self.at.logical), [True]*8)

   def test_add_property_no_property_type(self):
      # without propety_type, ambiguity between types PROPERTY_INT and PROPETY_LOGICAL
      self.at.add_property('int', [1]*self.at.n)
      self.check_property('int', T_INTEGER_A)
      self.assertEqual(list(self.at.int), [1]*8)

   # bcast() test relies on SIZEOF_ATOMS parameter, which is compiler and architecture dependant
   if QUIP_MAKEFILE['QUIPPY_FCOMPILER'] == 'gnu95' and '-DDARWIN' in QUIP_MAKEFILE['DEFINES']:
      def test_bcast(self):
         from quippy import MPI_context, atoms_bcast
         mpi = MPI_context()
         atoms_bcast(mpi, self.at)


class TestAtoms_Neighbour(QuippyTestCase):

   def setUp(self):
      self.a = 5.44
      self.bond_length = self.a * sqrt(3.0)/4
      
      self.at = diamond(self.a, 14)
      self.at.set_cutoff_factor(1.2)
      self.at.calc_connect()

   def test_n_neighbours(self):
      self.assertEqual([self.at.n_neighbours(i) for i in frange(self.at.n) ], [4]*8)

   def test_n_neighbours_max_dist_1(self):
      self.assertEqual([self.at.n_neighbours(i, max_dist=1.0) for i in frange(self.at.n) ], [0]*8)

   def test_n_neighbours_max_dist_2(self):
      self.assertEqual([self.at.n_neighbours(i, max_dist=3.0) for i in frange(self.at.n) ], [4]*8)

   def test_n_neighbours_max_factor_1(self):
      self.assertEqual([self.at.n_neighbours(i, max_factor=1.0) for i in frange(self.at.n) ], [0]*8)

   def test_n_neighbours_max_factor_2(self):
      self.assertEqual([self.at.n_neighbours(i, max_factor=1.2) for i in frange(self.at.n) ], [4]*8)

   def test_j(self):
      self.assertEqual(sorted([self.at.neighbour(1, i) for i in frange(self.at.n_neighbours(1))]), [2,4,6,8])

   def test_j_max_dist(self):
      self.assertEqual(sorted([self.at.neighbour(1, i, max_dist=1.0) for i in frange(self.at.n_neighbours(1))]), [0,0,0,0])

   def test_12(self):
      self.assertEqual(self.at.neighbour(1,4), 2)

   def test_distance_12(self):
      distance = farray(0.0)
      self.at.neighbour(1,4,distance=distance)
      self.assertAlmostEqual(distance, self.bond_length)

   def test_diff_norm_12(self):
      diff = fzeros(3)
      self.at.neighbour(1,4,diff=diff)
      self.assertAlmostEqual(diff.norm(), self.bond_length)

   def test_diff_12(self):
      # r_12 diff vector should be a*(0.25, 0.25, 0.25)
      diff = fzeros(3)
      self.at.neighbour(1,4,diff=diff)
      self.assertArrayAlmostEqual(diff, self.a*farray([0.25, 0.25, 0.25]))
                      
   def test_cosines_12(self):
      # r_12 cosines vector should be (1, 1, 1)/sqrt(3)
      cosines = fzeros(3)
      self.at.neighbour(1,4,cosines=cosines)
      self.assertArrayAlmostEqual(cosines, farray([1., 1., 1.])/sqrt(3.))

   def test_shift(self):
      # r_12 shift should be (0, 0, 0)
      shift = fzeros(3,dtype=numpy.int32)
      self.at.neighbour(1,4,shift=shift)
      self.assertArrayAlmostEqual(shift, [0.0, 0.0, 0.0])

   def test_is_nearest_neighbour_1(self):
      self.assert_(all([self.at.is_nearest_neighbour(1, n) for n in frange(self.at.n_neighbours(1))]))

   def test_is_nearest_neighbour_2(self):
      self.at.nneightol = 1.0
      self.assert_(not any([self.at.is_nearest_neighbour(1, n) for n in frange(self.at.n_neighbours(1))]))

   def test_high_level_n(self):
      self.assertEqual([len(self.at.neighbours[i]) for i in frange(self.at.n)],[4]*8)
      
   def test_high_level_j(self):
      self.assertEqual(sorted([n.j for n in self.at.neighbours[1]]), [2,4,6,8])

   def test_high_level_12(self):
      self.assertEqual(self.at.neighbours[1][4].j, 2)

   def test_high_level_distance_12(self):
      self.assertAlmostEqual(self.at.neighbours[1][4].distance, self.bond_length)

   def test_high_level_diff_12(self):
      self.assertArrayAlmostEqual(self.at.neighbours[1][4].diff, self.a*farray([0.25, 0.25, 0.25]))
      
   def test_high_level_cosines_12(self):
      self.assertArrayAlmostEqual(self.at.neighbours[1][4].cosines, farray([1., 1., 1.])/sqrt(3.))

   def test_high_level_shift_12(self):
      self.assertArrayAlmostEqual(self.at.neighbours[1][4].shift, [0, 0, 0])
      

class TestAtoms_CalcConnect_Hysteretic(QuippyTestCase):
   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 3, 3, 3)
      self.at.set_cutoff_factor(1.2)
      self.at.calc_connect()

   def test_equiv(self):
      at_h = self.at.copy()
      at_h.calc_connect_hysteretic()

      self.assertEqual([[n.j for n in self.at.neighbours[i]] for i in frange(self.at.n)],
                       [[n.j for n in at_h.neighbours[i]] for i in frange(at_h.n)] )

   def test_move_atom(self):
      move_atom = 66
      at_move = self.at.copy()
      at_move.pos[2,move_atom] -= 0.5
      at_move.calc_connect()

      self.assertNotEqual([[n.j for n in self.at.neighbours[i]] for i in frange(self.at.n)],
                          [[n.j for n in at_move.neighbours[i]] for i in frange(at_move.n)])

   def test_move_atom_hysteretic(self):
      move_atom = 66
      at_h = self.at.copy()
      at_h.set_cutoff_factor(1.2, 1.4)
      
      at_h.calc_connect_hysteretic()
      at_h.pos[2,move_atom] -= 0.5
      at_h.calc_connect_hysteretic()

      self.assertEqual([[n.j for n in self.at.neighbours[i]] for i in frange(self.at.n)],
                       [[n.j for n in at_h.neighbours[i]] for i in frange(at_h.n)] )

   def test_move_atom_hysteretic_alt_connect(self):
      move_atom = 66
      at_h = self.at.copy()
      at_h.set_cutoff_factor(1.2, 1.4)

      at_h.calc_connect_hysteretic(at_h.hysteretic_connect)
      at_h.pos[2,move_atom] -= 0.5
      at_h.calc_connect_hysteretic(at_h.hysteretic_connect)

      self.assertEqual([[n.j for n in self.at.neighbours[i]] for i in frange(self.at.n)],
                       [[n.j for n in at_h.hysteretic_neighbours[i]] for i in frange(at_h.n)])

   def test_move_atom_hysteretic_origin_extent(self):
      move_atom = 66
      at_h = self.at.copy()
      at_h.set_cutoff_factor(1.2, 1.4)

      active_mask = fzeros(at_h.n, dtype=bool)
      active_mask[move_atom] = True

      origin, extent = estimate_origin_extent(at_h, active_mask, 3.0)
      extent_inv = linalg.inv(extent)

      at_h.calc_connect_hysteretic(at_h.hysteretic_connect, origin=origin, extent=extent)
      at_h.pos[2,move_atom] -= 0.5
      at_h.calc_connect_hysteretic(at_h.hysteretic_connect, origin=origin, extent=extent)

      self.assertEqual(sorted([n.j for n in self.at.neighbours[move_atom]]),
                       sorted([n.j for n in at_h.hysteretic_neighbours[move_atom]]))
      


class TestGeometry(QuippyTestCase):

   def test_make_lattice_cubic(self):
      a = b = c = 1.0
      alpha = beta = gamma = 90.0*PI/180.0
      L = make_lattice(a,b,c,alpha,beta,gamma)
      self.assertEqual(get_lattice_params(L), (a,b,c,alpha,beta,gamma))
   
   def test_make_lattice_trigonal(self):
      a = b = 1.0
      c = 2.0
      alpha = beta = 90.0*PI/180.0
      gamma = 120.0*PI/180.0
      L = make_lattice(a,b,c,alpha,beta,gamma)
      self.assertArrayAlmostEqual(get_lattice_params(L), (a,b,c,alpha,beta,gamma))

   def test_make_lattice_triclinic(self):
      a = 1.0
      b = 2.0
      c = 3.0
      alpha = 60.0*PI/180.0
      beta  = 70.0*PI/180.0
      gamma = 80.0*PI/180.0
      L = make_lattice(a,b,c,alpha,beta,gamma)
      self.assertArrayAlmostEqual(get_lattice_params(L), (a,b,c,alpha,beta,gamma))

   def test_make_lattice_cubic_2(self):
      a = b = c = 1.0
      alpha = beta = gamma = 90.0*PI/180.0
      L = quippy._quippy.qp_make_lattice(a,b,c,alpha,beta,gamma)
      fvar(['a2','b2','c2','alpha2','beta2','gamma2'])
      quippy._quippy.qp_get_lattice_params(L,a2,b2,c2,alpha2,beta2,gamma2)
      self.assertArrayAlmostEqual((a2,b2,c2,alpha2,beta2,gamma2), (a,b,c,alpha,beta,gamma))
   
   def test_make_lattice_trigonal_2(self):
      a = b = 1.0
      c = 2.0
      alpha = beta = 90.0*PI/180.0
      gamma = 120.0*PI/180.0
      L = quippy._quippy.qp_make_lattice(a,b,c,alpha,beta,gamma)
      fvar(['a2','b2','c2','alpha2','beta2','gamma2'])
      quippy._quippy.qp_get_lattice_params(L,a2,b2,c2,alpha2,beta2,gamma2)
      self.assertArrayAlmostEqual((a2,b2,c2,alpha2,beta2,gamma2), (a,b,c,alpha,beta,gamma))

   def test_make_lattice_triclinic_2(self):
      a = 1.0
      b = 2.0
      c = 3.0
      alpha = 60.0*PI/180.0
      beta  = 70.0*PI/180.0
      gamma = 80.0*PI/180.0
      L = quippy._quippy.qp_make_lattice(a,b,c,alpha,beta,gamma)
      fvar(['a2','b2','c2','alpha2','beta2','gamma2'])
      quippy._quippy.qp_get_lattice_params(L,a2,b2,c2,alpha2,beta2,gamma2)
      self.assertArrayAlmostEqual((a2,b2,c2,alpha2,beta2,gamma2), (a,b,c,alpha,beta,gamma))

   def test_max_cutoff_1(self):
      dia = diamond(5.44, 14)
      sup = supercell(dia, 3, 1, 1)
      c = max_cutoff(sup.lattice)
      sup.set_cutoff(c)
      sup.calc_connect()
      self.assert_(sup.connect.is_min_image.all())

   def test_cell_volume(self):
      dia = diamond(5.44, 14)
      sup = supercell(dia, 3, 1, 1)
      self.assertAlmostEqual(sup.cell_volume(), abs(dot(sup.lattice[:,1],cross(sup.lattice[:,2],sup.lattice[:,3]))))

   def test_fit_box_in_cell(self):
      dia = diamond(5.44, 14)
      r  = 10.0
      na, nb, nc = fit_box_in_cell(r, r, r, dia.lattice)
      sup = supercell(dia, na, nb, nc)
      a = sup.lattice[:,1].norm()
      b = sup.lattice[:,2].norm()
      c = sup.lattice[:,3].norm()
      self.assert_(a >= 2.0*r and b >= 2.0*r and c >= 2.0*r)

   def test_centre_of_mass(self):
      dia = diamond(5.44, 14)
      dia.add_property('mass', [ElementMass[z] for z in dia.z])
      self.assertArrayAlmostEqual(dia.centre_of_mass(),
                                  farray([dia.diff_min_image(1,i) for i in frange(dia.n)]).T.mean(axis=2))
      

if __name__ == '__main__':
   unittest.main()
