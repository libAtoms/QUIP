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
from numpy import *
import unittest, quippy, os
from quippytest import *

def sort_for_comparison(t, t_padded):
      # sort subset list
      t.sort(int_cols=(1))

      # copy and sort original part of padded list
      t_padded_orig = Table()
      if (t.n >= 1):
	t_padded_orig.select(t_padded, row_list=frange(1,t.n))
	t_padded_orig.sort(int_cols=(1))

      # copy and sort new part of padded list
      t_padded_new = Table()
      if (t_padded.n >= t.n+1):
	t_padded_new.select(t_padded, row_list=frange(t.n+1,t_padded.n))
	t_padded_new.sort(int_cols=(1))

      # concatenate back into padded list
      t_padded.wipe()
      if (t_padded_orig.n > 0):
	t_padded.append(t_padded_orig)
      if (t_padded_new.n > 0):
	t_padded.append(t_padded_new)

class TestBFSStep(QuippyTestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 4, 4, 4)
      self.dia.set_cutoff(5.0)
      self.dia.calc_connect()
      self.at.set_cutoff(5.0)
      self.at.calc_connect()

      self.seed = Table(4,0,0,0)
      self.seed.append([1,0,0,0])

      self.max_r = 4.0

   def test_nneighb_only_true(self):
      # Recalc connectivity, only including nearest neighbours
      self.at.set_cutoff(3.0)
      self.at.calc_connect()

      crust = self.at.bfs_step(self.seed)
      
      # Results should be equal to the nearest neighbours of atom #1
      self.assertEqual(sorted(crust.int[1,:]),
                       sorted([neighb.j for neighb in self.at.neighbours[1]]))

   def test_nneighb_only_false(self):
      crust = self.at.bfs_step(self.seed, nneighb_only=False)

      # Results should be equal to all the neighbours of atom #1
      self.assertEqual(sorted(crust.int[1,:]),
                       sorted([neighb.j for neighb in self.at.neighbours[1]]))

   def test_multiple_images_1(self):
      t = Table(1,0,0,0)
      t.append([1,2,3])
      self.assert_(not multiple_images(t))

   def test_multiple_images_2(self):
      t = Table(1,0,0,0)
      t.append([1,2,1])
      self.assert_(multiple_images(t))

   def test_min_images_only_false(self):
      # min_images_only False -> repeated atom indices
      crust1 = self.dia.bfs_step(self.seed, nneighb_only=False, min_images_only=False)
      self.assert_(multiple_images(crust1))

   def test_min_images_only_true(self):
      # min_images_only True -> no repeated atom indices
      crust2 = self.dia.bfs_step(self.seed, nneighb_only=False, min_images_only=True)
      self.assert_(not multiple_images(crust2))

   def test_discard_non_min_images(self):
      crust1 = self.dia.bfs_step(self.seed, nneighb_only=False, min_images_only=False)
      crust2 = self.dia.bfs_step(self.seed, nneighb_only=False, min_images_only=True)
      
      # After discarding non-minimum-images, should have crust1 == crust2
      discard_non_min_images(crust1)
      self.assertArrayAlmostEqual(crust1.int, crust2.int)

   def test_max_r(self):
      crust = self.at.bfs_step(self.seed, nneighb_only=False, max_r=self.max_r)
      d = [self.at.distance(1, j, shift) for j, shift in zip(crust.int[1,:], crust.int[2:,:])]
      self.assert_(max(d) <= self.max_r)
      
   def test_property(self):
      mask = fzeros(self.at.n)
      mask[dot(self.at.g, self.at.pos)[1,:] < 0] = 1 # left half of cell

      # Assert hopping has been restricted to atoms where mask == 1
      crust = self.at.bfs_step(self.seed, nneighb_only=False, property=mask)
      self.assert_(mask[crust.int[1,:]].all())
      
class TestBFSGrow(QuippyTestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 4, 4, 4)
      self.dia.set_cutoff(5.0)
      self.dia.calc_connect()
      self.at.set_cutoff(5.0)
      self.at.calc_connect()

      self.seed = Table(4,0,0,0)
      self.seed.append([1,0,0,0])

   def test_bfs_grow_single_1hop(self):
      crust = self.at.bfs_step(self.seed)
      res   = self.at.bfs_grow_single(atom=1, n=1)

      # Assert res = [1] + crust
      self.assertEqual(list(res.int[1,:]), [1] + list(crust.int[1,:]))

   def test_bfs_grow_single_2hops(self):
      res1 = self.at.bfs_grow_single(atom=1, n=2)

      res2 = self.seed.copy()
      res2.append(self.at.bfs_step(res2))
      res2.append(self.at.bfs_step(res2))
      
      self.assertEqual(res1, res2)

   def test_bfs_grow_single_2hops_nneigb_only_false(self):
      res1 = self.at.bfs_grow_single(atom=1, n=2, nneighb_only=False)

      res2 = self.seed.copy()
      res2.append(self.at.bfs_step(res2, nneighb_only=False))
      res2.append(self.at.bfs_step(res2, nneighb_only=False))
      
      self.assertEqual(res1, res2)

   def test_bfs_grow_single_2hops_min_image_only_true(self):
      res1 = self.dia.bfs_grow_single(atom=1, n=2, min_images_only=True)

      res2 = self.seed.copy()
      res2.append(self.dia.bfs_step(res2, min_images_only=True))
      res2.append(self.dia.bfs_step(res2, min_images_only=True))
     
      self.assertEqual(res1, res2)

   def test_bfs_grow_list(self):
      res1 = self.at.bfs_grow_single(atom=1, n=2)
      res2 = self.seed.copy()
      self.at.bfs_grow_list(res2, n=2) # grows res2 in place
      self.assertEqual(res1, res2)


class TestCluster_TerminateFalse(QuippyTestCase):
   
   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 4, 4, 4)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=2)
      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.t = create_cluster_info_from_mark(self.at, "terminate=F cluster_allow_modification=F")

      sort_for_comparison(self.embed, self.t)

   def test_int_shape(self):
      self.assertEqual(self.t.int.shape,  (6, self.embed.n))

   def test_real_shape(self):
      self.assertEqual(self.t.real.shape, (4, self.embed.n))

   def test_i_shift(self):
      # First four int columns are (i, shift), and should match input
      self.assertArrayAlmostEqual(self.t.int[1:4,1:self.embed.n], self.embed.int[1:4,1:self.embed.n])

   def test_z(self):
      # Next int column is atomic number
      self.assertArrayAlmostEqual(self.t.int[5,1:self.embed.n], self.at.z[self.t.int[1,1:self.embed.n]])

   def test_term_index(self):
      # Final int column is term index
      self.assert_((self.t.int[6,1:self.embed.n] == 0).all())

   def test_pos(self):
      # First three real columns are positions, relative to atom #1 in cluster
      # for this test case, atom #1 is at origin
      self.assertArrayAlmostEqual(self.t.real[1:3,1:self.embed.n], self.at.pos[self.embed.int[1,1:self.embed.n]])

   def test_rescale(self):
      # Final real column is rescale factor, should be all 1.0 for terminate=False
      self.assertAlmostEqual(abs(self.t.real[4,1:self.embed.n] - 1.0).max(), 0.0)


class TestCluster_TerminateTrue(QuippyTestCase):
   
   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 4, 4, 4)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=2)
      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.cut_bonds = Table(8,0,0,0)
      self.t = create_cluster_info_from_mark(self.at, "terminate=T cluster_allow_modification=F", self.cut_bonds)

      sort_for_comparison(self.embed, self.t)

      self.n_term = self.t.n - self.embed.n
      self.term_atoms =  [122, 508, 510,  26, 124,  30,  98, 126, 100, 410, 508, 512,  26, 412,  32, 386, 416, 388,
                          482, 510, 512,  98, 486, 104, 386, 488, 390,  30,  32,   4, 100, 104,   6, 388, 390,   8]
      self.term_atoms.sort()

      self.cluster = carve_cluster(self.at, cluster_info=self.t, args_str="")
      self.cluster.calc_connect_hysteretic(1.2, 1.2)

   def test_i_shift(self):
      self.assertArrayAlmostEqual(self.t.int[1:4,1:self.embed.n], self.embed.int[1:4,1:self.embed.n])
      # Check atoms being terminated against an explicit list
      self.assertEqual(list(self.t.int[1,self.embed.n+1:]), self.term_atoms)

   def test_z(self):
      self.assertArrayAlmostEqual(self.t.int[5,1:self.embed.n], self.at.z[self.t.int[1,1:self.embed.n]])
      # All term atoms should have Z=1
      self.assertEqual(list(self.t.int[5,self.embed.n+1:]), [1]*self.n_term)

   def test_term_index(self):
      self.assert_((self.t.int[6,1:self.embed.n] == 0).all())

      # Termination indices should be in range 1 <= t <= n_embed
      self.assert_(self.t.int[6,self.embed.n+1:].min() >= 1)
      self.assert_(self.t.int[6,self.embed.n+1:].max() <= self.embed.n)

   def test_pos(self):
      self.assertArrayAlmostEqual(self.t.real[1:3,1:self.embed.n], self.at.pos[self.embed.int[1,1:self.embed.n]])

   def test_rescale(self):
      self.assertAlmostEqual(abs(self.t.real[4,1:self.embed.n] - 1.0).max(), 0.0)

      # rescale factors should all be equal to length(Si-H)/length(Si-Si)
      self.assertAlmostEqual(abs(self.t.real[4,self.embed.n+1:] - bond_length(14, 1)/bond_length(14, 14)).max(), 0)

   def test_cut_bonds_shape(self):
      self.assertEqual(self.cut_bonds.int.shape, (8,self.n_term))

   def test_cut_bonds(self):
      # cut_bonds is list of pairs (i,j) straddling boundary
      # all of (i,ishift) must be IN
      # all of (j,jshift) must be OUT
      for i, j, ishift, jshift in zip(self.cut_bonds.int[1,:],
                                      self.cut_bonds.int[2,:],
                                      self.cut_bonds.int[3:5,:],
                                      self.cut_bonds.int[6:8,:]):
         
         self.assert_(self.embed.find([i] + list(ishift)) != 0)
         self.assert_(self.embed.find([j] + list(jshift)) == 0)

   def test_cluster_coordination_si(self):
      # Check all Si atoms have 4 neighbours
      si_n_neighb = farray([len(self.cluster.neighbours[i]) for i in frange(self.cluster.n) if self.cluster.z[i] == 14])
      self.assert_((si_n_neighb == 4).all())

   def test_cluster_coordination_h(self):
      # Check all H atoms have 1 neighbour
      h_n_neighb  = farray([len(self.cluster.neighbours[i]) for i in frange(self.cluster.n) if self.cluster.z[i] == 1])
      self.assert_((h_n_neighb == 1).all())

   def test_cluster_hybrid_mark(self):
      # cluster hybrid_mark property should match that of original atoms object
      self.assert_((self.at.hybrid_mark[self.embed.int[1,:]] == self.cluster.hybrid_mark[1:self.embed.n]).all())

   def test_cluster_index(self):
      self.assert_((self.cluster.index == self.t.int[1,:]).all())

   def test_cluster_termindex(self):
      self.assert_((self.cluster.termindex == self.t.int[6,:]).all())

   def test_cluster_rescale(self):
      self.assertAlmostEqual(abs(self.cluster.rescale - self.t.real[4,:]).max(), 0.0)

   def test_cluster_shift(self):
      self.assert_((self.cluster.shift == self.t.int[2:4,:]).all())

   def test_cluster_ident(self):
      self.assert_((self.cluster.cluster_ident == self.t.str[1,:]).all())

                      

class TestCluster_Periodic(QuippyTestCase):
   
   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 4, 4, 1)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()

      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=2, min_images_only=True)
      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.t = create_cluster_info_from_mark(self.at, "terminate=T cluster_allow_modification=F cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T")
      self.cluster = carve_cluster(self.at, cluster_info=self.t, cluster_periodic_z=True)
      self.cluster.calc_connect_hysteretic(1.2, 1.2)
      
   def test_cluster_coordination_si(self):
      # Check all Si atoms have 4 neighbours
      si_n_neighb = farray([len(self.cluster.neighbours[i]) for i in frange(self.cluster.n) if self.cluster.z[i] == 14])
      self.assert_((si_n_neighb == 4).all())

   def test_cluster_coordination_h(self):
      # Check all H atoms have 1 neighbour
      h_n_neighb  = farray([len(self.cluster.neighbours[i]) for i in frange(self.cluster.n) if self.cluster.z[i] == 1])
      self.assert_((h_n_neighb == 1).all())

   def test_lattice(self):
      # cluster lattice should match underlying atoms in z direction
      self.assertAlmostEqual(self.cluster.lattice[3,3], self.at.lattice[3,3])


class TestCluster_EvenElectrons(QuippyTestCase):

   def setUp(self):
      dia = diamond(5.44, 14)
      self.at = supercell(dia, 4, 4, 4)
      self.at.set_cutoff(1.2*bond_length(14, 14))
      self.at.calc_connect()

      self.embedlist = Table()
      self.embedlist.append((1,0,0,0))
      self.at.bfs_grow_list(self.embedlist, 2, nneighb_only=True)

      self.at.add_property('hybrid',0)
      self.at.hybrid[self.embedlist.int[1,:]] = 1

      self.bufferlist = self.embedlist.copy()
      self.at.bfs_grow_list(self.bufferlist, 2, nneighb_only=True)

      self.at.add_property('hybrid_mark', 0)
      self.at.hybrid_mark[self.at.hybrid == 1] = HYBRID_ACTIVE_MARK
      for i in self.bufferlist.int[1,:]:
         if not self.at.hybrid_mark[i]:
            self.at.hybrid_mark[i] = HYBRID_BUFFER_MARK

      self.args = {'even_electrons':True}

   def test_14(self):
      # with z[1] = 14, no change needed
      cluster_info = create_cluster_info_from_mark(self.at, args_str(self.args))
      cluster = carve_cluster(self.at, args_str(self.args), cluster_info)
      self.assert_(cluster.z.sum() % 2 == 0)
      
   def test_15(self):
      # with z[1] = 15, need to remove an H atom
      self.at.z[1] = 15
      cluster_info = create_cluster_info_from_mark(self.at, args_str(self.args))
      cluster = carve_cluster(self.at, args_str(self.args), cluster_info)
      self.assert_(cluster.z.sum() % 2 == 0)

   def test_13(self):
      # with z[1] = 13, need to remove an H atom
      self.at.z[1] = 13
      cluster_info = create_cluster_info_from_mark(self.at, args_str(self.args))
      cluster = carve_cluster(self.at, args_str(self.args), cluster_info)
      self.assert_(cluster.z.sum() % 2 == 0)


class TestCluster_HollowSection(QuippyTestCase):

   def setUp(self):

      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 4, 4, 4)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=2)

      # Introduce a hollow section by unmarking atom 2
      self.hollow_atom = 2
      self.embed.delete(self.embed.find((self.hollow_atom,0,0,0)))

      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.t = create_cluster_info_from_mark(self.at, "terminate=F cluster_allow_modification=T")

      sort_for_comparison(self.embed, self.t)

   def test_int_shape(self):
      self.assertEqual(self.t.int.shape, (6, self.embed.n + 1))

   def test_real_shape(self):
      self.assertEqual(self.t.real.shape, (4, self.embed.n + 1))

   def test_i_shift(self):
      self.assertArrayAlmostEqual(self.t.int[1:4,1:self.embed.n], self.embed.int[1:4,1:self.embed.n])

   def test_i_shift_last(self):
      # Last entry should be [hollow_atom, 0,0,0]
      self.assertEqual(list(self.t.int[1:4,-1]), [self.hollow_atom, 0, 0, 0])

   def test_mark_last(self):
      # Should be marked as having been added by n_cut_bond section detector
      self.assertEqual(''.join(self.t.str[1,-1]).strip(), 'n_cut_bond')

   def test_z(self):
      self.assertArrayAlmostEqual(self.t.int[5,1:self.embed.n], self.at.z[self.t.int[1,1:self.embed.n]])

   def test_z_last(self):
      # Last atom in cluster should be hollow_atom
      self.assertEqual(self.t.int[5,-1], self.at.z[self.hollow_atom])

   def test_term_index(self):
      # Final int column is term index
      self.assert_((self.t.int[6,:] == 0).all())

   def test_pos(self):
      # First three real columns are positions, relative to atom #1 in cluster
      # for this test case, atom #1 is at origin
      self.assertArrayAlmostEqual(self.t.real[1:3,1:self.embed.n], self.at.pos[self.embed.int[1,1:self.embed.n]])

   def test_pos_last(self):
      self.assertArrayAlmostEqual(self.t.real[1:3,-1], self.at.pos[self.hollow_atom])

   def test_rescale(self):
      self.assertAlmostEqual(abs(self.t.real[4,:] - 1.0).max(), 0.0)

   def test_make_convex(self):
      # make_convex() should give same result as hollow section detector within create_cluster_info()
      convex = self.embed.copy()
      self.at.make_convex(convex)
      self.assertArrayAlmostEqual(convex.int, self.t.int[1:4,:])
      
      

class TestCluster_TerminationClash(QuippyTestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 3, 3, 3)

      # Introduce a termination clash by moving atoms 51 and 55
      # atoms 50 should then be added by the clash detector
      d = self.at.diff_min_image(51, 55)
      self.at.pos[51] = self.at.pos[51] + d/4.0
      self.at.pos[55] = self.at.pos[55] - d/4.0

      self.clash_atom = 50

      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=2)

      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.t = create_cluster_info_from_mark(self.at, "terminate=T reduce_n_cut_bonds=F cluster_allow_modification=T")

      sort_for_comparison(self.embed, self.t)

      self.clash_ind = -1
      for i in frange(1,self.t.n):
	if (str(self.t.str[1,i].stripstrings()) == 'clash'):
	  self.assertEqual(self.clash_ind,-1)
	  self.clash_ind = i

   def test_i_shift(self):
      self.assertArrayAlmostEqual(self.t.int[1:4,1:self.embed.n], self.embed.int[1:4,1:self.embed.n])

   def test_i_shift_clash(self):
      # Last entry should be [clash_atom, 0,-1,0]
      self.assertNotEqual(self.clash_ind,-1)
      # now 0 -1 0, because of no map_into_cell()
      self.assertEqual(list(self.t.int[1:4,self.clash_ind]), [self.clash_atom, 0, -1, 0])

# clash marking implicitly tested by test_i_shift_clash
#   def test_mark_clash(self):
#      # Should be marked as having been added by terminatin clash detector
#      self.assertEqual(str(self.t.str[1,self.embed.n+1]), 'clash')

   def test_z(self):
      self.assertArrayAlmostEqual(self.t.int[5,1:self.embed.n], self.at.z[self.t.int[1,1:self.embed.n]])

   def test_z_clash(self):
      # Last atom in cluster should be clash_atom
      self.assertEqual(self.t.int[5,self.clash_ind], self.at.z[self.clash_atom])

   def test_term_index(self):
      # Final int column is term index
      self.assert_((self.t.int[6,1:self.embed.n] == 0).all())
      self.assert_((self.t.int[6,self.clash_ind] == 0).all())

   def test_pos(self):
      # First three real columns are positions, relative to atom #1 in cluster
      # for this test case, atom #1 is at origin
      self.assertArrayAlmostEqual(self.t.real[1:3,1:self.embed.n], self.at.pos[self.embed.int[1,1:self.embed.n]])

   def test_pos_clash(self):
      self.assertArrayAlmostEqual(self.t.real[1:3,self.clash_ind], self.at.pos[self.clash_atom])

   def test_rescale(self):
      self.assertAlmostEqual(abs(self.t.real[4,1:self.embed.n] - 1.0).max(), 0.0)
      self.assertAlmostEqual(abs(self.t.real[4,self.clash_ind] - 1.0).max(), 0.0)


class TestCluster_Rescale(QuippyTestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 3, 3, 3)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=2)

      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.t = create_cluster_info_from_mark(self.at, "terminate=F cluster_allow_modification=F")

      self.r_scale = 0.9

      self.args = {'randomise_buffer' : False}
      self.cluster = carve_cluster(self.at, args_str(self.args), self.t)

      self.args_r_scale =  {'randomise_buffer' : False,
                            'do_rescale_r': True,
                            'r_scale': self.r_scale }
      self.cluster_r_scale = carve_cluster(self.at, args_str(self.args_r_scale), self.t)

   def test_lattice(self):
      # Check non-zero lattice components
      mask = self.cluster.lattice != 0.0
      lattice_ratio = self.cluster_r_scale.lattice[mask]/self.cluster.lattice[mask]
      self.assertAlmostEqual(abs(lattice_ratio - self.r_scale).max(), 0.0)

   def test_pos(self):
      # Check non-zero position components
      mask = self.cluster.pos != 0.0
      pos_ratio = self.cluster_r_scale.pos[mask]/self.cluster.pos[mask]
      self.assertAlmostEqual(abs(pos_ratio - self.r_scale).max(), 0.0)


class TestCluster_RandomiseBuffer(QuippyTestCase):

   def setUp(self):
      system_reseed_rng(2065775975)
      
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 3, 3, 3)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=1)

      self.buffer = self.embed.copy()
      self.at.bfs_grow_list(self.buffer, n=1)

      self.buffer_outer_layer = self.buffer.copy()
      self.at.bfs_grow_list(self.buffer_outer_layer, n=1)

      self.at.hybrid_mark[self.buffer_outer_layer.int[1,:]]= HYBRID_BUFFER_OUTER_LAYER_MARK
      self.at.hybrid_mark[self.buffer.int[1,:]]= HYBRID_BUFFER_MARK
      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.all_marked = self.buffer_outer_layer.copy()

      self.args = {'terminate' :True,
                   'randomise_buffer' : False}
      t = create_cluster_info_from_mark(self.at, args_str(self.args))
      self.cluster = carve_cluster(self.at, args_str(self.args), t)

      self.args_randomise_buffer = {'terminate' :True,
                                    'randomise_buffer' : True}
      t = create_cluster_info_from_mark(self.at, args_str(self.args_randomise_buffer))
      self.cluster_randomise_buffer = carve_cluster(self.at, args_str(self.args_randomise_buffer), t)

   def test_cluster_hybrid_mark(self):
      self.assert_((self.at.hybrid_mark[self.embed.int[1,:]] == self.cluster.hybrid_mark[1:self.embed.n]).all())
      
   def test_randomised_atoms(self):
      # Atoms which are moved should be those marked as HYBRID_BUFFER_OUTER_LAYER_MARK
      pos_diff = abs((self.cluster.pos - self.cluster_randomise_buffer.pos).norm())
      self.assert_((self.cluster.hybrid_mark[pos_diff > 0.0] == HYBRID_BUFFER_OUTER_LAYER_MARK).all())

   def test_regenerate_buffer_outer_layer(self):
      # Make a copy of at, with BUFFER_OUTER_LAYER mark removed.
      # Should be automatically regenerated by carve_cluster()
      at2 = self.at.copy()
      at2.calc_connect()
      at2.hybrid_mark[at2.hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK] = HYBRID_BUFFER_MARK

      t = create_cluster_info_from_mark(at2, args_str(self.args_randomise_buffer))
      cluster_randomise_buffer_2 = carve_cluster(at2, args_str(self.args_randomise_buffer), t)

      self.assertEqual(list(self.cluster_randomise_buffer.hybrid_mark),
                       list(cluster_randomise_buffer_2.hybrid_mark))


class TestCluster_SplitQM(QuippyTestCase):

   def setUp(self):
      self.dia = diamond(5.44, 14)
      self.at = supercell(self.dia, 3, 3, 3)
      # self.at.set_cutoff(5.0)
      self.at.set_cutoff(2.5)
      self.at.calc_connect()
      
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.embed = self.at.bfs_grow_single(1, n=1)

      self.buffer = self.embed.copy()
      self.at.bfs_grow_list(self.buffer, n=1)

      self.buffer_outer_layer = self.buffer.copy()
      self.at.bfs_grow_list(self.buffer_outer_layer, n=1)

      self.at.hybrid_mark[self.buffer_outer_layer.int[1,:]]= HYBRID_BUFFER_OUTER_LAYER_MARK
      self.at.hybrid_mark[self.buffer.int[1,:]]= HYBRID_BUFFER_MARK
      self.at.hybrid_mark[self.embed.int[1,:]]= HYBRID_ACTIVE_MARK

      self.all_marked = self.buffer_outer_layer.copy()

      self.args = {'terminate' :True,
                   'cluster_allow_modification': False,
                   'randomise_buffer' : False}

   def tearDown(self):
      os.remove('create_cluster_abort.xyz')
      os.remove('create_cluster_abort.xyz.idx')

   def test_split_qm_raises_runtime_error(self):
      # Split QM region by marking another atom
      self.at.hybrid_mark[107] = HYBRID_ACTIVE_MARK
      verbosity_push(PRINT_SILENT)
      self.assertRaises(RuntimeError, create_cluster_info_from_mark, self.at, args_str(self.args))
      verbosity_pop()
  

class TestCluster_Surface_Dia(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 3, 3, 3)
      lat = self.at.lattice
      lat[3,3] *= 2
      self.at.set_lattice(lat, scale_positions=False)
      self.at.set_cutoff(1.2*bond_length(14, 14))
      self.at.calc_connect()

      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)

      embed = self.at.bfs_grow_single(1, n=3)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T")
      self.cluster3 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T", cluster_info=self.t)

      embed = self.at.bfs_grow_single(1, n=4)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T")
      self.cluster4 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T", cluster_info=self.t)

      #self.at.print_xyz("tmp_at.xyz", all_properties=True)
      #self.cluster3.print_xyz("tmp_cluster3.xyz", all_properties=True)
      #self.cluster4.print_xyz("tmp_cluster4.xyz", all_properties=True)

   def test_surface_3_hops(self):
      self.assertEqual(self.cluster3.n,53)

   def test_surface_4_hops(self):
      self.assertEqual(self.cluster4.n,94)


class TestCluster_Surface_FCC(QuippyTestCase):

   def setUp(self):
      self.at = supercell(fcc_z111(4.05), 5, 5, 3)
      self.at.set_atoms([ 13 for x in frange(1,self.at.n) ])

      lat = self.at.lattice
      lat[3,3] *= 2
      self.at.set_lattice(lat, scale_positions=False)
      self.at.set_cutoff(3.0)
      self.at.calc_connect()

      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)

      embed = self.at.bfs_grow_single(1, n=1, nneighb_only=False)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F cluster_hopping_nneighb_only=F cluster_heuristics_nneighb_only=F")
      self.cluster1 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F randomise_buffer=F", cluster_info=self.t)

      embed = self.at.bfs_grow_single(1, n=2, nneighb_only=False)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F cluster_hopping_nneighb_only=F cluster_heuristics_nneighb_only=F")
      self.cluster2 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F randomise_buffer=F", cluster_info=self.t)

      hollow_atom = 28
      embed.delete(hollow_atom)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_mark(self.at, "cluster_allow_modification=T reduce_n_cut_bonds=T cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F cluster_hopping_nneighb_only=F cluster_heuristics_nneighb_only=F")
      self.cluster2h = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F randomise_buffer=F", cluster_info=self.t)

      # self.at.write("tmp_at.xyz")
      # self.cluster1.write("tmp_cluster1.xyz")
      # self.cluster2.write("tmp_cluster2.xyz")
      # self.cluster2h.write("tmp_cluster2h.xyz")

   def test_surface_1_hops(self):
      self.assertEqual(self.cluster1.n,10)

   def test_surface_2_hops(self):
      self.assertEqual(self.cluster2.n,37)

   def test_surface_2_hops_hollow(self):
      self.assertEqual(self.cluster2h.n,37)


class TestCluster_CrackTip(QuippyTestCase):

   def setUp(self):
      slab = supercell(diamond(5.44, 14), 10, 10, 1)
      slab.map_into_cell()
      width = slab.lattice[1,1]
      notch_width = 0.5*width
      notch_height = notch_width/2.0

      slab.lattice[1,1] += 50.0
      slab.lattice[2,2] += 50.0
      slab.set_lattice(slab.lattice, scale_positions=False)

      mask = fzeros(slab.n, dtype=int)
      mask[:] = 1
      for i in frange(slab.n):
         if ((slab.pos[2,i] < -(0.5*notch_height/notch_width*(slab.pos[1,i]+width/2.0)) + notch_height/2.0) and
             (slab.pos[2,i] >  (0.5*notch_height/notch_width*(slab.pos[1,i]+width/2.0)) - notch_height/2.0)):
            mask[i] = 0

      self.at = slab.select(mask)
      self.at.set_cutoff(5.0)
      self.at.calc_connect()
      
      embed = self.at.bfs_grow_single(1, n=5, nneighb_only=True)
      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK

      self.t1 = create_cluster_info_from_mark(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T cluster_hopping_nneighb_only=F cluster_heuristics_nneighb_only=F cluster_allow_modification=F")
      self.cluster1 = carve_cluster(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T randomise_buffer=F", cluster_info=self.t1)

      self.t2 = create_cluster_info_from_mark(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T cluster_hopping_nneighb_only=F cluster_heuristics_nneighb_only=F cluster_allow_modification=F")
      self.cluster2 = carve_cluster(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T randomise_buffer=F", cluster_info=self.t2)

      #AtomsList([self.at,self.cluster1,self.cluster2]).show()
      #raw_input()

   def test_cluster_n(self):
      self.assertEqual(self.cluster1.n, 97)

   def test_cluster_n_term(self):
      self.assertEqual((self.cluster1.z == 1).sum(), 40)

   def test_cluster_lattice_x(self):
      self.assertAlmostEqual(self.cluster1.lattice[1,1], 28.0)

   def test_cluster_lattice_y(self):
      self.assertAlmostEqual(self.cluster1.lattice[2,2], 28.0)

   def test_cluster_lattice_z(self):
      self.assertAlmostEqual(self.cluster1.lattice[3,3], self.at.lattice[3,3])
      
   def test_cluster_nneighb_only_table(self):
      self.assertEqual(sorted(self.t1.int[1,:]), sorted(self.t2.int[1,:]))


def mark_atoms(at, nneighb_only=True, alt_connect=None):
   at.add_property('hybrid_mark', HYBRID_NO_MARK)
   embed = at.bfs_grow_single(1, n=1, nneighb_only=nneighb_only, alt_connect=alt_connect)

   buffer = embed.copy()
   at.bfs_grow_list(buffer, n=1, nneighb_only=nneighb_only, alt_connect=alt_connect)

   buffer_outer_layer = buffer.copy()
   at.bfs_grow_list(buffer_outer_layer, n=1, nneighb_only=nneighb_only, alt_connect=alt_connect)

   at.hybrid_mark[:] = HYBRID_NO_MARK
   at.hybrid_mark[buffer_outer_layer.int[1,:]]= HYBRID_BUFFER_OUTER_LAYER_MARK
   at.hybrid_mark[buffer.int[1,:]]= HYBRID_BUFFER_MARK
   at.hybrid_mark[embed.int[1,:]]= HYBRID_ACTIVE_MARK


class TestCluster_HystereticConnect_Equiv(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 3, 3, 3)

      self.at.set_cutoff(1.2*bond_length(14, 14))
      self.at.calc_connect()
      self.at.calc_connect_hysteretic(1.2, 1.4, self.at.hysteretic_connect)

      mark_atoms(self.at)

      self.args = {'terminate': True,
                   'randomise_buffer' : False,
                   'hysteretic_connect': False}
      self.t = create_cluster_info_from_mark(self.at, args_str(self.args))
      self.cluster = carve_cluster(self.at, args_str(self.args), cluster_info=self.t)

      self.args_hysteretic = {'terminate' :True,
                              'randomise_buffer' : False,
                              'hysteretic_connect': True,
                              'cluster_hopping_nneighb_only': False,
                              'cluster_heuristics_nneighb_only': False}
      self.t_hysteretic = create_cluster_info_from_mark(self.at, args_str(self.args_hysteretic))
      self.cluster_hysteretic = carve_cluster(self.at, args_str(self.args_hysteretic), cluster_info=self.t_hysteretic)

   def test_cluster_info(self):
      self.assertEqual(self.t, self.t_hysteretic)

   def test_cluster(self):
      self.assertEqual(self.cluster, self.cluster_hysteretic)


class TestCluster_HystereticConnect_MoveAtom(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 3, 3, 3)
      self.at.set_cutoff(1.2*bond_length(14, 14))      
      self.at.calc_connect()
      mark_atoms(self.at)

      self.move_atom = 66

      self.at_hysteretic = self.at.copy()
      self.at_hysteretic.calc_connect_hysteretic(1.2, 1.4, self.at_hysteretic.hysteretic_connect)

      self.at_move = self.at.copy()
      
      self.at_move.pos[2,self.move_atom] -= 0.5
      self.at_hysteretic.pos[2,self.move_atom] -= 0.5

      self.at_move.calc_connect()
      mark_atoms(self.at_move)

      self.args = {'terminate': True,
                   'randomise_buffer' : False,
                   'hysteretic_connect': False}

      self.t = create_cluster_info_from_mark(self.at, args_str(self.args))
      self.cluster = carve_cluster(self.at, args_str(self.args), cluster_info=self.t)

      self.t_move = create_cluster_info_from_mark(self.at_move, args_str(self.args))
      self.cluster_move = carve_cluster(self.at_move, args_str(self.args), cluster_info=self.t)

      self.at_hysteretic.calc_connect_hysteretic(1.2, 1.4, self.at_hysteretic.hysteretic_connect)
      mark_atoms(self.at_hysteretic, nneighb_only=False, alt_connect=self.at_hysteretic.hysteretic_connect)

      self.args_hysteretic = {'terminate' :True,
                              'randomise_buffer' : False,
                              'hysteretic_connect': True,
                              'cluster_hopping_nneighb_only': False,
                              'cluster_heuristics_nneighb_only': False}
      self.t_hysteretic = create_cluster_info_from_mark(self.at_hysteretic, args_str(self.args_hysteretic))
      self.cluster_hysteretic = carve_cluster(self.at_hysteretic, args_str(self.args_hysteretic), cluster_info=self.t_hysteretic)

   def test_mark(self):
      self.assertEqual(self.at.hybrid_mark[self.move_atom], HYBRID_BUFFER_OUTER_LAYER_MARK)
   
   def test_mark_move(self):
      self.assertEqual(self.at_move.hybrid_mark[self.move_atom], HYBRID_NO_MARK)

   def test_mark_hysteretic(self):
      self.assertEqual(self.at_hysteretic.hybrid_mark[self.move_atom], HYBRID_BUFFER_OUTER_LAYER_MARK)      

   def test_cluster_info_int(self):
      self.assertArrayAlmostEqual(self.t.int, self.t_hysteretic.int)
      
   def test_cluster_info_real(self):
      r_diff = abs(self.t.real - self.t_hysteretic.real) > 0.0
      r_diff_index = r_diff.nonzero()[1]

      # find the index of the atom we moved in the cluster table
      row1 = self.t.subtable(rows=frange(self.t.n), intcols=[1])
      move_atom_index = row1.find(self.move_atom)

      # Consider each entry in cluster_info which doesn't match original
      for i in r_diff_index:
         if self.t.int[6,i] == 0:
            # it's not a termination atom, so should be the one we moved
            self.assertEqual(self.t.int[1,i],self.move_atom)
         else:
            # it's a Hydrogen, so it should be attached to the atom we moved
            self.assertEqual(self.t.int[6,i],move_atom_index)


class TestCluster_SilicaTetrahedra(QuippyTestCase):
      def setUp(self):
            at = Atoms("""1622
Lattice="133.83785769       0.00000000       0.00000000       0.00000000     103.28524004       0.00000000       0.00000000       0.00000000       4.84038097" Properties=species:S:1:pos:R:3:Z:I:1:primitive_index:I:1:travel:I:3:orig_index:I:1:map_shift:I:3
Si             -1.96237006      1.87175075      1.02530679      14       2       0       0       0       1       0       0       0
Si              1.98983936     -0.04984537      1.07378369      14       3       0       0       0       2       0       0       0
O               0.45760858     -0.61800904      1.67220462       8       8       0       0       0       3       0       0       0
O               1.12886215     -2.27006565     -1.31033472       8       9       0       0       0       4       0       0       0
O              -1.79663016      1.25318665     -0.37839279       8       7       0       0       0       5       0       0       0
Si             -0.08153228     -1.69610667     -2.28824454      14       1       0       0       1       6       0       0       0
O              -1.23114359     -1.28912567     -1.37794458       8       4       0       0       1       7       0       0       0
O              -0.58847356      2.33685864      1.73112570       8       5       0       0       1       8       0       0       0
O               1.70533483      0.53183548     -0.38946669       8       6       0       0       1       9       0       0       0
O              -3.01573391     -2.40171874      1.10026603       8       9       0       0       0      10       0       0       0
Si             -4.20805398     -1.72857557      0.23566184      14       1       0       0       1      11       0       0       0
O              -2.42886765      0.62954805      2.09001100       8       6       0       0       1      12       0       0       0
Si              2.18700887      1.76853455     -1.31507538      14       2       0       1       0      13       0       0       0
O               3.70993066      2.33396064     -0.81374032       8       5       0       1       1      14       0       0       0
O               2.44937251      1.20155950      2.16565613       8       7       1       1       0      17       0       0       0
O               3.03523522     -1.16834817      1.05861145       8       4       1       1       1      18       0       0       0
Si             -2.01769707      7.08382833      1.05652135      14       2       0       0       0      19       0       0       0
Si              1.88539476      5.38286770      1.06750851      14       3       0       0       0      20       0       0       0
O               0.61773462      4.69403356      1.61336279       8       8       0       0       0      21       0       0       0
O               1.26204795      2.99877040     -1.40320851       8       9       0       0       0      22       0       0       0
O              -1.76329650      6.61204702     -0.37966973       8       7       0       0       0      23       0       0       0
Si              0.03846452      3.45306029     -2.28984001      14       1       0       0       1      24       0       0       0
O              -1.09311780      4.18175263     -1.24090029       8       4       0       0       1      25       0       0       0
O              -0.47866437      7.61122876      1.72712124       8       5       0       0       1      26       0       0       0
O               1.64299018      5.92201259     -0.43083828       8       6       0       0       1      27       0       0       0
O              -3.11568630      2.99300145      1.17187314       8       9       0       0       0      28       0       0       0
Si             -4.25792680      3.50093876      0.09532339      14       1       0       0       1      29       0       0       0
O              -2.39028305      5.99602386      2.11158076       8       6       0       0       1      30       0       0       0
Si              2.16578779      7.12916118     -1.28223222      14       2       0       1       0      31       0       0       0
O               3.61852870      7.63576869     -0.83175949       8       5       0       1       1      32       0       0       0
Si             -2.15454052      5.38255044     -1.21573147      14       3       0       0       0      33       0       0       0
O              -3.62767412      4.66489434     -0.67420786       8       8       0       0       0      34       0       0       0
O               2.50428697      6.53101265      2.13630956       8       7       1       1       0      35       0       0       0
O               2.98660038      4.05607163      1.07568591       8       4       1       1       1      36       0       0       0
Si             -1.90847428     12.44941949      1.06798369      14       2       0       0       0      37       0       0       0
Si              2.01269960     10.60764403      1.22281574      14       3       0       0       0      38       0       0       0
O               0.46124336     10.10221028      1.73453700       8       8       0       0       0      39       0       0       0
O               1.21842449      8.33930765     -1.40999242       8       9       0       0       0      40       0       0       0
O              -1.62820912     11.85332871     -0.37799570       8       7       0       0       0      41       0       0       0
Si              0.04963313      8.91237762     -2.18434942      14       1       0       0       1      42       0       0       0
O              -1.21230683      9.52826319     -1.31274111       8       4       0       0       1      43       0       0       0
O              -0.49972984     13.07945621      1.62049549       8       5       0       0       1      44       0       0       0
O               1.64662654     11.26571549     -0.29407369       8       6       0       0       1      45       0       0       0
O              -3.03391855      8.28636275      1.05412265       8       9       0       0       0      46       0       0       0
Si             -4.11017378      8.86246334      0.07348334      14       1       0       0       1      47       0       0       0
O              -2.50718872     11.31672862      2.01380300       8       6       0       0       1      48       0       0       0
Si              2.20367349     12.47380902     -1.19812302      14       2       0       1       0      49       0       0       0
O               3.56031289     12.95494076     -0.71788100       8       5       0       1       1      50       0       0       0
Si             -2.34224503     10.67019444     -1.35266345      14       3       0       0       0      51       0       0       0
O              -3.59854007     10.15407853     -0.71189322       8       8       0       0       0      52       0       0       0
O               2.40823717     11.77862959      2.10492659       8       7       1       1       0      53       0       0       0
O               2.97660659      9.53470796      1.12992099       8       4       1       1       1      54       0       0       0
Si             -1.94745983     17.82578717      1.03750311      14       2       0       0       0      55       0       0       0
Si              1.85644697     15.98544023      1.07198623      14       3       0       0       0      56       0       0       0
O               0.57240923     15.30835770      1.74523676       8       8       0       0       0      57       0       0       0
O               1.09774135     13.55416676     -1.27146988       8       9       0       0       0      58       0       0       0
O              -1.70812023     17.26154333     -0.41908807       8       7       0       0       0      59       0       0       0
Si              0.03995167     14.11822830     -2.28706893      14       1       0       0       1      60       0       0       0
O              -1.24792875     14.76862818     -1.30587364       8       4       0       0       1      61       0       0       0
O              -0.57477451     18.37028569      1.73635848       8       5       0       0       1      62       0       0       0
O               1.76486137     16.50969270     -0.41254801       8       6       0       0       1      63       0       0       0
O              -3.04066140     13.72088221      1.02210550       8       9       0       0       0      64       0       0       0
Si             -4.28190219     14.14059239      0.14977744      14       1       0       0       1      65       0       0       0
O              -2.45873226     16.57823893      2.15000991       8       6       0       0       1      66       0       0       0
Si              2.27887735     17.78534894     -1.34469437      14       2       0       1       0      67       0       0       0
O               3.63728500     18.29267575     -0.82746029       8       5       0       1       1      68       0       0       0
Si             -2.15719947     15.97883381     -1.38560848      14       3       0       0       0      69       0       0       0
O              -3.61699617     15.39904263     -0.77878705       8       8       0       0       0      70       0       0       0
O               2.54659193     17.23038807      1.97506460       8       7       1       1       0      71       0       0       0
O               2.95099219     14.82098458      1.11837922       8       4       1       1       1      72       0       0       0
Si             -2.01791919     23.08378069      1.15790477      14       2       0       0       0      73       0       0       0
Si              2.00725272     21.34748452      1.11549013      14       3       0       0       0      74       0       0       0
O               0.53712107     20.72004042      1.71929269       8       8       0       0       0      75       0       0       0
O               1.10146647     18.95833258     -1.31821306       8       9       0       0       0      76       0       0       0
O              -1.64047962     22.53725321     -0.31218075       8       7       0       0       0      77       0       0       0
Si             -0.03044683     19.47416042     -2.30065864      14       1       0       0       1      78       0       0       0
O              -1.11823321     20.18026234     -1.26949557       8       4       0       0       1      79       0       0       0
O              -0.51238611     23.69557483      1.67139676       8       5       0       0       1      80       0       0       0
O               1.80264256     21.99489368     -0.31157751       8       6       0       0       1      81       0       0       0
O              -2.99876281     18.89766690      1.09111199       8       9       0       0       0      82       0       0       0
Si             -4.14170440     19.53117464      0.17459179      14       1       0       0       1      83       0       0       0
O              -2.49613188     21.97795546      1.97584625       8       6       0       0       1      84       0       0       0
Si              2.20010400     23.11055736     -1.29893522      14       2       0       1       0      85       0       0       0
O               3.71348442     23.72019354     -0.82066366       8       5       0       1       1      86       0       0       0
Si             -2.26050018     21.29683005     -1.36829791      14       3       0       0       0      87       0       0       0
O              -3.54133917     20.76622685     -0.76245108       8       8       0       0       0      88       0       0       0
O               2.49859985     22.59475060      2.08278217       8       7       1       1       0      89       0       0       0
O               3.04783640     20.06268964      1.07038218       8       4       1       1       1      90       0       0       0
Si             -2.02282527    -24.80497230      1.05549173      14       2       0       1       0      91       0       0       0
Si              2.03567231    -26.54871242      1.17167549      14       3       0       1       0      92       0       0       0
O               0.54756813     25.96080268      1.76368764       8       8       0       0       0      93       0       0       0
O               1.18407355     24.22475618     -1.27168959       8       9       0       0       0      94       0       0       0
O              -1.74435398    -25.52770105     -0.43723843       8       7       0       1       0      95       0       0       0
Si             -0.05259460     24.91993005     -2.17860519      14       1       0       0       1      96       0       0       0
O              -1.08771973     25.36741975     -1.40491524       8       4       0       0       1      97       0       0       0
O              -0.54284981    -24.18203155      1.65833173       8       5       0       1       1      98       0       0       0
O               1.62655540    -26.07996653     -0.26687871       8       6       0       1       1      99       0       0       0
O              -3.00006827     24.36275482      1.05440846       8       9       0       0       0     100       0       0       0
Si             -4.22656390     24.80416087      0.09899430      14       1       0       0       1     101       0       0       0
O              -2.46049993    -26.06971449      2.11998697       8       6       0       1       1     102       0       0       0
Si              2.32330684    -24.79748299     -1.21107545      14       2       0       2       0     103       0       0       0
O               3.60576113    -24.26266031     -0.78396748       8       5       0       2       1     104       0       0       0
Si             -2.16441139     26.63505164     -1.35030546      14       3       0       0       0     105       0       0       0
O              -3.60963515     26.06041888     -0.81262407       8       8       0       0       0     106       0       0       0
O               2.36907579    -25.43639933      2.16018632       8       7       1       2       0     107       0       0       0
O               3.10938470     25.43924028      1.11709936       8       4       1       1       1     108       0       0       0
Si             -1.90759508    -19.59920956      1.19541096      14       2       0       1       0     109       0       0       0
Si              2.04523903    -21.39130286      1.08506241      14       3       0       1       0     110       0       0       0
O               0.49205381    -21.96858151      1.58235661       8       8       0       1       0     111       0       0       0
O               1.09582278    -23.76450357     -1.32765721       8       9       0       1       0     112       0       0       0
O              -1.63946225    -20.11091715     -0.41122925       8       7       0       1       0     113       0       0       0
Si              0.03936363    -23.02072648     -2.24299009      14       1       0       1       1     114       0       0       0
O              -1.25519560    -22.53371573     -1.31082052       8       4       0       1       1     115       0       0       0
O              -0.50649310    -19.03812809      1.62992697       8       5       0       1       1     116       0       0       0
O               1.77675658    -20.74422005     -0.39116050       8       6       0       1       1     117       0       0       0
O              -2.93889235    -23.65460020      1.00660717       8       9       0       1       0     118       0       0       0
Si             -4.14958793    -23.10518999      0.15748407      14       1       0       1       1     119       0       0       0
O              -2.54105735    -20.66722068      2.03227709       8       6       0       1       1     120       0       0       0
Si              2.15009445    -19.44954410     -1.27289905      14       2       0       2       0     121       0       0       0
O               3.53954937    -18.87920811     -0.67716855       8       5       0       2       1     122       0       0       0
Si             -2.34045712    -21.27329612     -1.29826899      14       3       0       1       0     123       0       0       0
O              -3.68506496    -21.90476731     -0.83257950       8       8       0       1       0     124       0       0       0
O               2.53791714    -20.08611282      2.15423530       8       7       1       2       0     125       0       0       0
O               3.03126347    -22.60285594      1.16343783       8       4       1       2       1     126       0       0       0
Si             -1.93521007    -14.11290172      1.09649717      14       2       0       1       0     127       0       0       0
Si              2.01418064    -16.05115654      1.08953919      14       3       0       1       0     128       0       0       0
O               0.48242037    -16.48905102      1.75865614       8       8       0       1       0     129       0       0       0
O               1.12423237    -18.28505578     -1.34343918       8       9       0       1       0     130       0       0       0
O              -1.70882948    -14.73410049     -0.32549083       8       7       0       1       0     131       0       0       0
Si              0.09391915    -17.80795947     -2.14971795      14       1       0       1       1     132       0       0       0
O              -1.15618959    -17.23318483     -1.41107010       8       4       0       1       1     133       0       0       0
O              -0.65385625    -13.63890612      1.63399552       8       5       0       1       1     134       0       0       0
O               1.78610787    -15.33096671     -0.40640413       8       6       0       1       1     135       0       0       0
O              -3.07690240    -18.35640184      1.18624552       8       9       0       1       0     136       0       0       0
Si             -4.13047185    -17.67466266      0.18363101      14       1       0       1       1     137       0       0       0
O              -2.48833942    -15.43030506      1.97828192       8       6       0       1       1     138       0       0       0
Si              2.31882679    -14.18551891     -1.24574743      14       2       0       2       0     139       0       0       0
O               3.68373151    -13.52852645     -0.71230440       8       5       0       2       1     140       0       0       0
Si             -2.23487038    -15.95687951     -1.26436225      14       3       0       1       0     141       0       0       0
O              -3.70701613    -16.66120345     -0.65797541       8       8       0       1       0     142       0       0       0
O               2.37789508    -14.80141701      2.14245953       8       7       1       2       0     143       0       0       0
O               2.94776594    -17.21006533      1.15242273       8       4       1       2       1     144       0       0       0
Si             -2.03686409     -8.87443703      1.06792660      14       2       0       1       0     145       0       0       0
Si              1.86727329    -10.70009877      1.05099640      14       3       0       1       0     146       0       0       0
O               0.50272758    -11.31732685      1.61524699       8       8       0       1       0     147       0       0       0
O               1.21099343    -12.98429596     -1.37731608       8       9       0       1       0     148       0       0       0
O              -1.71403794     -9.42933247     -0.33905836       8       7       0       1       0     149       0       0       0
Si             -0.09748880    -12.35979482     -2.21417805      14       1       0       1       1     150       0       0       0
O              -1.25400078    -11.92201529     -1.31293454       8       4       0       1       1     151       0       0       0
O              -0.46413044     -8.23641819      1.58233920       8       5       0       1       1     152       0       0       0
O               1.73960143    -10.00696318     -0.38001875       8       6       0       1       1     153       0       0       0
O              -3.08283861    -12.98864377      1.19531792       8       9       0       1       0     154       0       0       0
Si             -4.24142774    -12.49223825      0.21250598      14       1       0       1       1     155       0       0       0
O              -2.41714638    -10.12003278      2.04148717       8       6       0       1       1     156       0       0       0
Si              2.18520432     -8.89717398     -1.26247185      14       2       0       2       0     157       0       0       0
O               3.65057896     -8.28778244     -0.79609922       8       5       0       2       1     158       0       0       0
Si             -2.15211442    -10.66156042     -1.38130156      14       3       0       1       0     159       0       0       0
O              -3.56494142    -11.24828711     -0.74350773       8       8       0       1       0     160       0       0       0
O               2.41388110     -9.39763628      2.02549304       8       7       1       2       0     161       0       0       0
O               3.02225120    -11.77447112      1.09810282       8       4       1       2       1     162       0       0       0
Si             -1.92178332     -3.48899854      1.20379233      14       2       0       1       0     163       0       0       0
Si              2.01157064     -5.25296372      1.15374148      14       3       0       1       0     164       0       0       0
O               0.48894399     -5.89905334      1.71573315       8       8       0       1       0     165       0       0       0
O               1.07211087     -7.59662253     -1.39790074       8       9       0       1       0     166       0       0       0
O              -1.81197691     -4.11577155     -0.35671966       8       7       0       1       0     167       0       0       0
Si             -0.00273098     -7.08390912     -2.26591844      14       1       0       1       1     168       0       0       0
O              -1.14712974     -6.55610804     -1.28904660       8       4       0       1       1     169       0       0       0
O              -0.54381542     -2.93137741      1.62861635       8       5       0       1       1     170       0       0       0
O               1.68945382     -4.83339978     -0.37603833       8       6       0       1       1     171       0       0       0
O              -2.92577402     -7.76252047      1.19797579       8       9       0       1       0     172       0       0       0
Si             -4.17594384     -7.17865024      0.09636578      14       1       0       1       1     173       0       0       0
O              -2.48988331     -4.71275266      2.00752078       8       6       0       1       1     174       0       0       0
Si              2.25630807     -3.61728282     -1.34898784      14       2       0       2       0     175       0       0       0
O               3.68337406     -2.99938889     -0.70979032       8       5       0       2       1     176       0       0       0
Si             -2.21227184     -5.41966167     -1.25405946      14       3       0       1       0     177       0       0       0
O              -3.68132327     -5.96650822     -0.69217018       8       8       0       1       0     178       0       0       0
O               2.37917851     -4.04423556      2.04374809       8       7       1       2       0     179       0       0       0
O               2.99751409     -6.51813971      1.17225309       8       4       1       2       1     180       0       0       0
Si              6.46395481      1.70638957      1.16286272      14       2       0       0       0     181       0       0       0
Si             10.27159503      0.00460171      1.16604539      14       3       0       0       0     182       0       0       0
O               8.99192259     -0.49602921      1.73546705       8       8       0       0       0     183       0       0       0
O               9.46548870     -2.28776855     -1.29263316       8       9       0       0       0     184       0       0       0
O               6.60986745      1.16985676     -0.36820575       8       7       0       0       0     185       0       0       0
Si              8.44473668     -1.81940116     -2.34062706      14       1       0       0       1     186       0       0       0
O               7.15407699     -1.13909761     -1.29932692       8       4       0       0       1     187       0       0       0
O               7.85409224      2.30900262      1.68382204       8       5       0       0       1     188       0       0       0
O              10.11470117      0.53557756     -0.33248967       8       6       0       0       1     189       0       0       0
O               5.44525011     -2.44408848      1.06634313       8       9       0       0       0     190       0       0       0
Si              4.14455224     -1.70334585      0.18438006      14       1       0       0       1     191       0       0       0
O               5.83106644      0.57431087      2.00187047       8       6       0       0       1     192       0       0       0
Si             10.67691738      1.78853192     -1.26682519      14       2       0       1       0     193       0       0       0
O              12.07252855      2.40824794     -0.83478647       8       5       0       1       1     194       0       0       0
Si              6.21923052      0.05929777     -1.22178165      14       3       0       0       0     195       0       0       0
O               4.83587301     -0.59522946     -0.80589357       8       8       0       0       0     196       0       0       0
O              10.75248763      1.12826909      2.02630813       8       7       1       1       0     197       0       0       0
O              11.35415544     -1.12569701      1.17871982       8       4       1       1       1     198       0       0       0
Si              6.46136698      7.01089826      1.20011183      14       2       0       0       0     199       0       0       0
Si             10.41936350      5.41786966      1.19064685      14       3       0       0       0     200       0       0       0
O               9.03140670      4.75918849      1.63645645       8       8       0       0       0     201       0       0       0
O               9.56351204      3.04269600     -1.23794912       8       9       0       0       0     202       0       0       0
O               6.74092570      6.54168067     -0.29829081       8       7       0       0       0     203       0       0       0
Si              8.41209901      3.50268125     -2.21517022      14       1       0       0       1     204       0       0       0
O               7.19339168      4.09630641     -1.23845743       8       4       0       0       1     205       0       0       0
O               7.82339820      7.69601427      1.72769562       8       5       0       0       1     206       0       0       0
O              10.13062920      6.00567073     -0.36770137       8       6       0       0       1     207       0       0       0
O               5.37019810      2.90451737      1.06662348       8       9       0       0       0     208       0       0       0
Si              4.27712917      3.55432901      0.15217929      14       1       0       0       1     209       0       0       0
O               6.00162091      5.88201271      1.97982164       8       6       0       0       1     210       0       0       0
Si             10.60887931      7.15864374     -1.32411090      14       2       0       1       0     211       0       0       0
O              12.11074545      7.75000305     -0.77851243       8       5       0       1       1     212       0       0       0
Si              6.07245848      5.31343599     -1.22367218      14       3       0       0       0     213       0       0       0
O               4.67474450      4.68461742     -0.67138115       8       8       0       0       0     214       0       0       0
O              10.77671315      6.43338548      2.16014963       8       7       1       1       0     215       0       0       0
O              11.40841546      4.07681007      1.06741229       8       4       1       1       1     216       0       0       0
Si              6.34231627     12.38770285      1.16205744      14       2       0       0       0     217       0       0       0
Si             10.41915427     10.67680579      1.15000666      14       3       0       0       0     218       0       0       0
O               9.00679457     10.14693362      1.61955807       8       8       0       0       0     219       0       0       0
O               9.45518520      8.22178844     -1.24442052       8       9       0       0       0     220       0       0       0
O               6.59304368     11.94012934     -0.40983089       8       7       0       0       0     221       0       0       0
Si              8.35019848      8.91750751     -2.19724402      14       1       0       0       1     222       0       0       0
O               7.14158951      9.51399727     -1.38948526       8       4       0       0       1     223       0       0       0
O               7.88809820     13.11452992      1.59243325       8       5       0       0       1     224       0       0       0
O              10.07871444     11.22499925     -0.30858683       8       6       0       0       1     225       0       0       0
O               5.44831001      8.30245895      1.10703606       8       9       0       0       0     226       0       0       0
Si              4.20829572      8.91671593      0.16198683      14       1       0       0       1     227       0       0       0
O               5.99364062     11.24435412      2.12317533       8       6       0       0       1     228       0       0       0
Si             10.69885108     12.39304514     -1.32504841      14       2       0       1       0     229       0       0       0
O              12.09917141     12.93903069     -0.67845124       8       5       0       1       1     230       0       0       0
Si              6.12964340     10.61403444     -1.20436258      14       3       0       0       0     231       0       0       0
O               4.75874796     10.12608749     -0.65329606       8       8       0       0       0     232       0       0       0
O              10.90941803     11.87664813      2.02306481       8       7       1       1       0     233       0       0       0
O              11.36442679      9.52587056      1.08477384       8       4       1       1       1     234       0       0       0
Si              6.43606201     17.76517079      1.16839539      14       2       0       0       0     235       0       0       0
Si             10.34938229     15.93199202      1.20236660      14       3       0       0       0     236       0       0       0
O               8.84447813     15.42643331      1.71333154       8       8       0       0       0     237       0       0       0
O               9.57087768     13.58546736     -1.26679159       8       9       0       0       0     238       0       0       0
O               6.64453970     17.08292134     -0.27597054       8       7       0       0       0     239       0       0       0
Si              8.44525788     14.29006358     -2.27339717      14       1       0       0       1     240       0       0       0
O               7.17843910     14.77820234     -1.30990906       8       4       0       0       1     241       0       0       0
O               7.90477780     18.33170392      1.67802332       8       5       0       0       1     242       0       0       0
O              10.12533855     16.61855344     -0.32313476       8       6       0       0       1     243       0       0       0
O               5.29709056     13.71435794      1.18080921       8       9       0       0       0     244       0       0       0
Si              4.24947873     14.12390045      0.19356916      14       1       0       0       1     245       0       0       0
O               5.94028046     16.50484331      2.16152692       8       6       0       0       1     246       0       0       0
Si             10.54068605     17.72621818     -1.25036996      14       2       0       1       0     247       0       0       0
O              12.00213168     18.29616261     -0.76606286       8       5       0       1       1     248       0       0       0
Si              6.17348111     15.98467101     -1.28056938      14       3       0       0       0     249       0       0       0
O               4.83319895     15.46573285     -0.66465701       8       8       0       0       0     250       0       0       0
O              10.74977772     17.24781008      2.09123929       8       7       1       1       0     251       0       0       0
O              11.50263485     14.81329442      1.14375816       8       4       1       1       1     252       0       0       0
Si              6.41332737     23.03830543      1.15565182      14       2       0       0       0     253       0       0       0
Si             10.42431694     21.25139455      1.09062270      14       3       0       0       0     254       0       0       0
O               8.96231470     20.77804303      1.61173629       8       8       0       0       0     255       0       0       0
O               9.61795056     18.90063307     -1.33685546       8       9       0       0       0     256       0       0       0
O               6.67274673     22.51343778     -0.42249731       8       7       0       0       0     257       0       0       0
Si              8.40940101     19.58777750     -2.18152011      14       1       0       0       1     258       0       0       0
O               7.17250539     20.06830225     -1.27632593       8       4       0       0       1     259       0       0       0
O               7.84356673     23.60990200      1.65128056       8       5       0       0       1     260       0       0       0
O              10.15977146     21.80677424     -0.41837473       8       6       0       0       1     261       0       0       0
O               5.43195128     18.88270403      1.17693230       8       9       0       0       0     262       0       0       0
Si              4.24290932     19.60820913      0.20121645      14       1       0       0       1     263       0       0       0
O               5.95748567     21.81975428      2.16685978       8       6       0       0       1     264       0       0       0
Si             10.70280008     23.03499969     -1.37657750      14       2       0       1       0     265       0       0       0
O              11.98262915     23.58155321     -0.80015149       8       5       0       1       1     266       0       0       0
Si              6.17182511     21.34936734     -1.38006609      14       3       0       0       0     267       0       0       0
O               4.73965338     20.67570556     -0.70017281       8       8       0       0       0     268       0       0       0
O              10.82611787     22.48825794      1.99945427       8       7       1       1       0     269       0       0       0
O              11.38944529     20.05051267      1.17810434       8       4       1       1       1     270       0       0       0
Si              6.41830722    -24.89100189      1.12168932      14       2       0       1       0     271       0       0       0
Si             10.29303702     26.62750751      1.12862726      14       3       0       0       0     272       0       0       0
O               8.96575535     26.05248192      1.57957338       8       8       0       0       0     273       0       0       0
O               9.53823085     24.26628263     -1.29560732       8       9       0       0       0     274       0       0       0
O               6.75256485    -25.53062512     -0.40712832       8       7       0       1       0     275       0       0       0
Si              8.41934595     24.92096140     -2.30408600      14       1       0       0       1     276       0       0       0
O               7.15510412     25.37861764     -1.23935433       8       4       0       0       1     277       0       0       0
O               7.81154026    -24.29685579      1.74867199       8       5       0       1       1     278       0       0       0
O              10.07843636    -26.10349771     -0.45159824       8       6       0       1       1     279       0       0       0
O               5.41981628     24.24191585      1.16823323       8       9       0       0       0     280       0       0       0
Si              4.26789254     24.90092951      0.25386343      14       1       0       0       1     281       0       0       0
O               5.99864261    -26.05433159      2.06992874       8       6       0       1       1     282       0       0       0
Si             10.66244815    -24.92842794     -1.26556231      14       2       0       2       0     283       0       0       0
O              11.97213460    -24.28788738     -0.69556403       8       5       0       2       1     284       0       0       0
Si              6.08814069     26.61052819     -1.25970908      14       3       0       0       0     285       0       0       0
O               4.82045818     26.09878837     -0.75729788       8       8       0       0       0     286       0       0       0
O              10.85524024    -25.50595739      2.08222966       8       7       1       2       0     287       0       0       0
O              11.46771524     25.52017406      1.11658432       8       4       1       1       1     288       0       0       0
Si              6.45250901    -19.53267541      1.19181951      14       2       0       1       0     289       0       0       0
Si             10.37618874    -21.31397210      1.07148175      14       3       0       1       0     290       0       0       0
O               8.94517380    -21.94808994      1.72251462       8       8       0       1       0     291       0       0       0
O               9.60635818    -23.73777063     -1.25595063       8       9       0       1       0     292       0       0       0
O               6.67157603    -20.12828862     -0.29324480       8       7       0       1       0     293       0       0       0
Si              8.37692775    -22.99390580     -2.18524526      14       1       0       1       1     294       0       0       0
O               7.18392792    -22.57895835     -1.25256204       8       4       0       1       1     295       0       0       0
O               7.73984193    -18.99649778      1.65457247       8       5       0       1       1     296       0       0       0
O              10.08535194    -20.66923670     -0.32907393       8       6       0       1       1     297       0       0       0
O               5.27621373    -23.64204639      1.15282154       8       9       0       1       0     298       0       0       0
Si              4.28845631    -23.13451900      0.22735678      14       1       0       1       1     299       0       0       0
O               5.85419679    -20.76643541      2.15926014       8       6       0       1       1     300       0       0       0
Si             10.69227593    -19.55393170     -1.20178202      14       2       0       2       0     301       0       0       0
O              12.09308528    -18.93263291     -0.73617424       8       5       0       2       1     302       0       0       0
Si              6.17535000    -21.27022826     -1.30364926      14       3       0       1       0     303       0       0       0
O               4.76337747    -21.82736703     -0.68201610       8       8       0       1       0     304       0       0       0
O              10.78189160    -20.08643014      2.01153140       8       7       1       2       0     305       0       0       0
O              11.40065903    -22.55814319      1.15505205       8       4       1       2       1     306       0       0       0
Si              6.40704572    -14.12425324      1.06179573      14       2       0       1       0     307       0       0       0
Si             10.27456822    -15.92933293      1.14171892      14       3       0       1       0     308       0       0       0
O               8.95987385    -16.60882680      1.69243561       8       8       0       1       0     309       0       0       0
O               9.51463537    -18.38548464     -1.33004929       8       9       0       1       0     310       0       0       0
O               6.74106554    -14.71789476     -0.27974114       8       7       0       1       0     311       0       0       0
Si              8.46571287    -17.83653098     -2.23809038      14       1       0       1       1     312       0       0       0
O               7.16681380    -17.15032343     -1.31987211       8       4       0       1       1     313       0       0       0
O               7.74837832    -13.66912371      1.65469445       8       5       0       1       1     314       0       0       0
O              10.04835965    -15.33540086     -0.44020624       8       6       0       1       1     315       0       0       0
O               5.36020535    -18.37357307      1.07420233       8       9       0       1       0     316       0       0       0
Si              4.23761002    -17.71507346      0.09196674      14       1       0       1       1     317       0       0       0
O               5.85951311    -15.46235735      2.08066691       8       6       0       1       1     318       0       0       0
Si             10.57951090    -14.21310816     -1.26869675      14       2       0       2       0     319       0       0       0
O              12.08842732    -13.56394733     -0.69556471       8       5       0       2       1     320       0       0       0
Si              6.12083709    -16.00617437     -1.27029002      14       3       0       1       0     321       0       0       0
O               4.77281561    -16.63798829     -0.81170901       8       8       0       1       0     322       0       0       0
O              10.79643054    -14.76433198      2.11758800       8       7       1       2       0     323       0       0       0
O              11.32767647    -17.12635043      1.11154432       8       4       1       2       1     324       0       0       0
Si              6.44118271     -8.92483481      1.04382748      14       2       0       1       0     325       0       0       0
Si             10.32498843    -10.66677564      1.02835882      14       3       0       1       0     326       0       0       0
O               8.87190581    -11.16168931      1.63081856       8       8       0       1       0     327       0       0       0
O               9.57654265    -12.94355663     -1.28041504       8       9       0       1       0     328       0       0       0
O               6.66660645     -9.42416690     -0.26609360       8       7       0       1       0     329       0       0       0
Si              8.41461324    -12.41504256     -2.15348738      14       1       0       1       1     330       0       0       0
O               7.14429925    -11.78501514     -1.26550843       8       4       0       1       1     331       0       0       0
O               7.87650563     -8.21461333      1.62971916       8       5       0       1       1     332       0       0       0
O              10.05417694    -10.17181541     -0.32823745       8       6       0       1       1     333       0       0       0
O               5.36415888    -13.03850821      1.04799166       8       9       0       1       0     334       0       0       0
Si              4.24895209    -12.42094074      0.22354174      14       1       0       1       1     335       0       0       0
O               5.82234031    -10.02553230      2.00071262       8       6       0       1       1     336       0       0       0
Si             10.59477668     -8.96337879     -1.24550300      14       2       0       2       0     337       0       0       0
O              12.09267742     -8.31619126     -0.67540336       8       5       0       2       1     338       0       0       0
Si              6.20660927    -10.70112377     -1.29289580      14       3       0       1       0     339       0       0       0
O               4.75716671    -11.25707708     -0.83599420       8       8       0       1       0     340       0       0       0
O              10.92301897     -9.46600870      2.05731762       8       7       1       2       0     341       0       0       0
O              11.49724422    -11.81648855      1.16933816       8       4       1       2       1     342       0       0       0
Si              6.50329271     -3.55310211      1.05934394      14       2       0       1       0     343       0       0       0
Si             10.40038332     -5.30443036      1.03110522      14       3       0       1       0     344       0       0       0
O               8.86112016     -5.89676712      1.71760804       8       8       0       1       0     345       0       0       0
O               9.58244227     -7.70376328     -1.30580990       8       9       0       1       0     346       0       0       0
O               6.69850783     -4.23608515     -0.44423915       8       7       0       1       0     347       0       0       0
Si              8.40583393     -7.16579539     -2.14980977      14       1       0       1       1     348       0       0       0
O               7.27531319     -6.55592626     -1.30455157       8       4       0       1       1     349       0       0       0
O               7.73311431     -3.03822263      1.68810066       8       5       0       1       1     350       0       0       0
O              10.16381485     -4.72406014     -0.29466995       8       6       0       1       1     351       0       0       0
O               5.40663424     -7.66712986      1.18229937       8       9       0       1       0     352       0       0       0
Si              4.14383256     -7.10468002      0.08555381      14       1       0       1       1     353       0       0       0
O               5.86928808     -4.77962314      2.03724185       8       6       0       1       1     354       0       0       0
Si             10.65636485     -3.52169773     -1.35677870      14       2       0       2       0     355       0       0       0
O              11.98715244     -3.02956155     -0.75609086       8       5       0       2       1     356       0       0       0
Si              6.19698508     -5.42734781     -1.37850349      14       3       0       1       0     357       0       0       0
O               4.73870784     -5.86626668     -0.83550401       8       8       0       1       0     358       0       0       0
O              10.90624720     -4.14837965      2.02771710       8       7       1       2       0     359       0       0       0
O              11.50537738     -6.48039926      1.18264149       8       4       1       2       1     360       0       0       0
Si             14.78184674      1.71318334      1.21850181      14       2       0       0       0     361       0       0       0
Si             18.64194717     -0.08747604      1.17613742      14       3       0       0       0     362       0       0       0
O              17.40541893     -0.60910937      1.57121751       8       8       0       0       0     363       0       0       0
O              17.94427905     -2.32170645     -1.27137494       8       9       0       0       0     364       0       0       0
O              15.13653827      1.11200726     -0.44825277       8       7       0       0       0     365       0       0       0
Si             16.79489242     -1.83700782     -2.31307368      14       1       0       0       1     366       0       0       0
O              15.69795851     -1.20065305     -1.41631525       8       4       0       0       1     367       0       0       0
O              16.21919999      2.38513559      1.60089472       8       5       0       0       1     368       0       0       0
O              18.49883806      0.63475593     -0.31515017       8       6       0       0       1     369       0       0       0
O              13.67115987     -2.43895875      1.14217581       8       9       0       0       0     370       0       0       0
Si             12.54681969     -1.75124995      0.14761915      14       1       0       0       1     371       0       0       0
O              14.20690371      0.56220319      2.11377738       8       6       0       0       1     372       0       0       0
Si             18.94576900      1.81019274     -1.22733623      14       2       0       1       0     373       0       0       0
O              20.45100615      2.30161783     -0.82304966       8       5       0       1       1     374       0       0       0
Si             14.42753305      0.00558376     -1.37931561      14       3       0       0       0     375       0       0       0
O              13.10657298     -0.65029857     -0.77473687       8       8       0       0       0     376       0       0       0
O              19.30595234      1.09517978      2.04831594       8       7       1       1       0     377       0       0       0
O              19.82960434     -1.19405353      1.05705328       8       4       1       1       1     378       0       0       0
Si             14.85998550      7.04153423      1.17422705      14       2       0       0       0     379       0       0       0
Si             18.65107336      5.35619267      1.08775065      14       3       0       0       0     380       0       0       0
O              17.31353810      4.75110788      1.73109819       8       8       0       0       0     381       0       0       0
O              17.95228265      3.02260733     -1.23615638       8       9       0       0       0     382       0       0       0
O              15.13954054      6.49855957     -0.35002147       8       7       0       0       0     383       0       0       0
Si             16.70835278      3.65023976     -2.26090241      14       1       0       0       1     384       0       0       0
O              15.54872249      4.14012336     -1.40009909       8       4       0       0       1     385       0       0       0
O              16.16242045      7.68519444      1.62580478       8       5       0       0       1     386       0       0       0
O              18.54358736      5.96314358     -0.27532171       8       6       0       0       1     387       0       0       0
O              13.73546935      2.98215591      1.05096130       8       9       0       0       0     388       0       0       0
Si             12.48706454      3.52484240      0.19764154      14       1       0       0       1     389       0       0       0
O              14.30336324      5.91110405      2.06654069       8       6       0       0       1     390       0       0       0
Si             18.95608199      7.04252635     -1.37792257      14       2       0       1       0     391       0       0       0
O              20.46159956      7.66814644     -0.66924746       8       5       0       1       1     392       0       0       0
Si             14.48091939      5.23451635     -1.38324291      14       3       0       0       0     393       0       0       0
O              13.03878556      4.73629873     -0.71298070       8       8       0       0       0     394       0       0       0
O              19.25216957      6.47486241      2.16193884       8       7       1       1       0     395       0       0       0
O              19.85983688      4.14242359      1.14364719       8       4       1       1       1     396       0       0       0
Si             14.79800030     12.33907710      1.08857515      14       2       0       0       0     397       0       0       0
Si             18.76331766     10.66806536      1.03808144      14       3       0       0       0     398       0       0       0
O              17.32738554     10.03091650      1.68701673       8       8       0       0       0     399       0       0       0
O              17.94415753      8.31164039     -1.27171465       8       9       0       0       0     400       0       0       0
O              15.13347706     11.89504556     -0.38088207       8       7       0       0       0     401       0       0       0
Si             16.82836372      8.86894090     -2.21603653      14       1       0       0       1     402       0       0       0
O              15.50564052      9.49823579     -1.30509977       8       4       0       0       1     403       0       0       0
O              16.25186328     12.96573575      1.57789061       8       5       0       0       1     404       0       0       0
O              18.48655317     11.34149891     -0.26851238       8       6       0       0       1     405       0       0       0
O              13.83358407      8.24471796      1.07777360       8       9       0       0       0     406       0       0       0
Si             12.51996520      8.96903902      0.14928405      14       1       0       0       1     407       0       0       0
O              14.25507440     11.18522304      2.13681039       8       6       0       0       1     408       0       0       0
Si             18.92998475     12.43273555     -1.29994150      14       2       0       1       0     409       0       0       0
O              20.38760248     13.08218915     -0.75890740       8       5       0       1       1     410       0       0       0
Si             14.44325687     10.70036562     -1.26882573      14       3       0       0       0     411       0       0       0
O              13.16002171     10.11343301     -0.80676690       8       8       0       0       0     412       0       0       0
O              19.29407573     11.83657970      2.07319935       8       7       1       1       0     413       0       0       0
O              19.83042889      9.51643943      1.13548614       8       4       1       1       1     414       0       0       0
Si             14.87503935     17.84573739      1.13114875      14       2       0       0       0     415       0       0       0
Si             18.68183683     15.97384236      1.16514564      14       3       0       0       0     416       0       0       0
O              17.39432754     15.40643402      1.69077662       8       8       0       0       0     417       0       0       0
O              17.92172092     13.53787338     -1.25013694       8       9       0       0       0     418       0       0       0
O              15.10639217     17.15517021     -0.27681712       8       7       0       0       0     419       0       0       0
Si             16.79287077     14.20284506     -2.32721898      14       1       0       0       1     420       0       0       0
O              15.52786052     14.87315597     -1.36669925       8       4       0       0       1     421       0       0       0
O              16.12566321     18.26858925      1.58200055       8       5       0       0       1     422       0       0       0
O              18.39721941     16.61653062     -0.40687290       8       6       0       0       1     423       0       0       0
O              13.65741928     13.63376049      1.14595578       8       9       0       0       0     424       0       0       0
Si             12.56739542     14.23070443      0.08911271      14       1       0       0       1     425       0       0       0
O              14.22432838     16.47625968      1.99856311       8       6       0       0       1     426       0       0       0
Si             18.91924663     17.69109615     -1.24264341      14       2       0       1       0     427       0       0       0
O              20.45285101     18.30327441     -0.76025613       8       5       0       1       1     428       0       0       0
Si             14.50042321     15.88654642     -1.33877293      14       3       0       0       0     429       0       0       0
O              13.18265206     15.43467249     -0.72032422       8       8       0       0       0     430       0       0       0
O              19.22396405     17.25106718      2.14492689       8       7       1       1       0     431       0       0       0
O              19.78108614     14.80352209      1.00929352       8       4       1       1       1     432       0       0       0
Si             14.81741158     23.09890349      1.16201228      14       2       0       0       0     433       0       0       0
Si             18.70814679     21.33724120      1.15710355      14       3       0       0       0     434       0       0       0
O              17.38578991     20.68980123      1.74467174       8       8       0       0       0     435       0       0       0
O              17.91644567     18.92653666     -1.37914364       8       9       0       0       0     436       0       0       0
O              15.12252743     22.52052754     -0.31009275       8       7       0       0       0     437       0       0       0
Si             16.67498375     19.59976466     -2.32600660      14       1       0       0       1     438       0       0       0
O              15.65395078     20.06297176     -1.40709278       8       4       0       0       1     439       0       0       0
O              16.13511411     23.64050658      1.70722853       8       5       0       0       1     440       0       0       0
O              18.51725559     21.94741135     -0.44812433       8       6       0       0       1     441       0       0       0
O              13.64614714     18.92137431      1.01063225       8       9       0       0       0     442       0       0       0
Si             12.58340411     19.44744353      0.08079427      14       1       0       0       1     443       0       0       0
O              14.20226358     21.83009849      2.06062131       8       6       0       0       1     444       0       0       0
Si             19.10392872     23.02904365     -1.24522319      14       2       0       1       0     445       0       0       0
O              20.48295912     23.67402424     -0.65068132       8       5       0       1       1     446       0       0       0
Si             14.58064650     21.34112101     -1.35147554      14       3       0       0       0     447       0       0       0
O              13.21631923     20.80540893     -0.69064512       8       8       0       0       0     448       0       0       0
O              19.16259473     22.58604691      2.12658838       8       7       1       1       0     449       0       0       0
O              19.72261481     20.08408459      1.17622098       8       4       1       1       1     450       0       0       0
Si             14.75932341    -24.81576934      1.08624461      14       2       0       1       0     451       0       0       0
Si             18.73784518     26.59910562      1.04997842      14       3       0       0       0     452       0       0       0
O              17.40631915     26.03160431      1.73040964       8       8       0       0       0     453       0       0       0
O              18.03395821     24.24025958     -1.30492538       8       9       0       0       0     454       0       0       0
O              15.08381469    -25.46843829     -0.36972049       8       7       0       1       0     455       0       0       0
Si             16.81393162     24.80938849     -2.33505848      14       1       0       0       1     456       0       0       0
O              15.68103455     25.41497364     -1.41710667       8       4       0       0       1     457       0       0       0
O              16.19016058    -24.35733030      1.63427402       8       5       0       1       1     458       0       0       0
O              18.44352301    -25.96528529     -0.39835952       8       6       0       1       1     459       0       0       0
O              13.67876032     24.28251087      1.09691716       8       9       0       0       0     460       0       0       0
Si             12.49957006     24.95106549      0.25812604      14       1       0       0       1     461       0       0       0
O              14.36418211    -26.00865691      2.05087611       8       6       0       1       1     462       0       0       0
Si             18.94304844    -24.87490429     -1.33922030      14       2       0       2       0     463       0       0       0
O              20.35473821    -24.23078982     -0.68007186       8       5       0       2       1     464       0       0       0
Si             14.58527095    -26.61701540     -1.33540710      14       3       0       1       0     465       0       0       0
O              13.09538181     26.01864062     -0.66600436       8       8       0       0       0     466       0       0       0
O              19.32996657    -25.44684457      2.14022158       8       7       1       2       0     467       0       0       0
O              19.74439653     25.49363887      1.12953627       8       4       1       1       1     468       0       0       0
Si             14.82852774    -19.57685864      1.17292845      14       2       0       1       0     469       0       0       0
Si             18.67623940    -21.41159145      1.05485797      14       3       0       1       0     470       0       0       0
O              17.40731181    -21.94517453      1.69164418       8       8       0       1       0     471       0       0       0
O              17.90927850    -23.61565267     -1.35301059       8       9       0       1       0     472       0       0       0
O              15.02830867    -20.03665381     -0.35491566       8       7       0       1       0     473       0       0       0
Si             16.77003008    -23.06520445     -2.30170973      14       1       0       1       1     474       0       0       0
O              15.68695320    -22.57897810     -1.23802651       8       4       0       1       1     475       0       0       0
O              16.11617563    -18.90167619      1.65337608       8       5       0       1       1     476       0       0       0
O              18.53692861    -20.72398386     -0.35730446       8       6       0       1       1     477       0       0       0
O              13.82227808    -23.62422795      1.14605299       8       9       0       1       0     478       0       0       0
Si             12.48113610    -23.01854834      0.09784303      14       1       0       1       1     479       0       0       0
O              14.30607024    -20.63985662      2.06815222       8       6       0       1       1     480       0       0       0
Si             18.97888792    -19.50701384     -1.23416754      14       2       0       2       0     481       0       0       0
O              20.47377658    -18.88233862     -0.71365279       8       5       0       2       1     482       0       0       0
Si             14.53315643    -21.26557840     -1.33246693      14       3       0       1       0     483       0       0       0
O              13.17130994    -21.85340686     -0.70698808       8       8       0       1       0     484       0       0       0
O              19.22405378    -20.05089428      2.15442713       8       7       1       2       0     485       0       0       0
O              19.81673060    -22.42093899      1.18072768       8       4       1       2       1     486       0       0       0
Si             14.77981645    -14.12898214      1.17184037      14       2       0       1       0     487       0       0       0
Si             18.61521572    -16.00642196      1.07661364      14       3       0       1       0     488       0       0       0
O              17.34050843    -16.65349884      1.64556680       8       8       0       1       0     489       0       0       0
O              17.89825601    -18.28056302     -1.30673446       8       9       0       1       0     490       0       0       0
O              15.04504372    -14.78634212     -0.37644149       8       7       0       1       0     491       0       0       0
Si             16.85085271    -17.75333145     -2.25563745      14       1       0       1       1     492       0       0       0
O              15.58447087    -17.11198053     -1.25971286       8       4       0       1       1     493       0       0       0
O              16.15946504    -13.71826942      1.59185735       8       5       0       1       1     494       0       0       0
O              18.46781163    -15.32524884     -0.25805196       8       6       0       1       1     495       0       0       0
O              13.83044651    -18.43955096      1.14269625       8       9       0       1       0     496       0       0       0
Si             12.66305233    -17.72310777      0.24634924      14       1       0       1       1     497       0       0       0
O              14.37052598    -15.38334403      2.12387522       8       6       0       1       1     498       0       0       0
Si             19.09773126    -14.18104751     -1.31036033      14       2       0       2       0     499       0       0       0
O              20.34745877    -13.64985022     -0.66061893       8       5       0       2       1     500       0       0       0
Si             14.43894956    -16.03047961     -1.34617983      14       3       0       1       0     501       0       0       0
O              13.07872631    -16.56842138     -0.83987148       8       8       0       1       0     502       0       0       0
O              19.19635841    -14.71696380      1.98691486       8       7       1       2       0     503       0       0       0
O              19.72602197    -17.16395766      1.04191976       8       4       1       2       1     504       0       0       0
Si             14.84321591     -8.84824356      1.13120952      14       2       0       1       0     505       0       0       0
Si             18.69606328    -10.66877216      1.06811392      14       3       0       1       0     506       0       0       0
O              17.38842755    -11.14735127      1.73813349       8       8       0       1       0     507       0       0       0
O              17.87299325    -12.96944927     -1.41655031       8       9       0       1       0     508       0       0       0
O              15.08196551     -9.49462392     -0.30347109       8       7       0       1       0     509       0       0       0
Si             16.71564266    -12.40251647     -2.15298526      14       1       0       1       1     510       0       0       0
O              15.61401973    -11.83996831     -1.36498141       8       4       0       1       1     511       0       0       0
O              16.17978335     -8.21933101      1.69583242       8       5       0       1       1     512       0       0       0
O              18.40494576    -10.11730611     -0.30969535       8       6       0       1       1     513       0       0       0
O              13.68983702    -13.09634026      1.06716392       8       9       0       1       0     514       0       0       0
Si             12.60849917    -12.49496826      0.16002898      14       1       0       1       1     515       0       0       0
O              14.30111001    -10.09490768      2.11442098       8       6       0       1       1     516       0       0       0
Si             19.01749316     -8.81006138     -1.37056704      14       2       0       2       0     517       0       0       0
O              20.33682275     -8.34197327     -0.82038397       8       5       0       2       1     518       0       0       0
Si             14.53559629    -10.63177625     -1.32171581      14       3       0       1       0     519       0       0       0
O              13.23261505    -11.21706873     -0.80343935       8       8       0       1       0     520       0       0       0
O              19.29685011     -9.48399500      2.02893649       8       7       1       2       0     521       0       0       0
O              19.80994332    -11.93842456      1.03178557       8       4       1       2       1     522       0       0       0
Si             14.91699489     -3.50338256      1.04272062      14       2       0       1       0     523       0       0       0
Si             18.79365852     -5.36937890      1.20659022      14       3       0       1       0     524       0       0       0
O              17.34427883     -5.93356345      1.64024800       8       8       0       1       0     525       0       0       0
O              17.84245243     -7.71398184     -1.30522480       8       9       0       1       0     526       0       0       0
O              14.96030562     -4.10435989     -0.34321152       8       7       0       1       0     527       0       0       0
Si             16.68064588     -7.18743201     -2.25671531      14       1       0       1       1     528       0       0       0
O              15.58952254     -6.50998960     -1.26885393       8       4       0       1       1     529       0       0       0
O              16.18507280     -3.06453389      1.63738461       8       5       0       1       1     530       0       0       0
O              18.58862288     -4.64718512     -0.27415578       8       6       0       1       1     531       0       0       0
O              13.67767094     -7.66974431      1.01537908       8       9       0       1       0     532       0       0       0
Si             12.67023617     -7.12787583      0.18838158      14       1       0       1       1     533       0       0       0
O              14.35515497     -4.81444333      2.03534973       8       6       0       1       1     534       0       0       0
Si             19.07853753     -3.51297367     -1.32401701      14       2       0       2       0     535       0       0       0
O              20.42840188     -2.88645657     -0.78810158       8       5       0       2       1     536       0       0       0
Si             14.55519967     -5.36490225     -1.27607982      14       3       0       1       0     537       0       0       0
O              13.12864083     -5.85995666     -0.80223958       8       8       0       1       0     538       0       0       0
O              19.13611036     -4.17958080      2.02640045       8       7       1       2       0     539       0       0       0
O              19.76167001     -6.42852252      1.05887135       8       4       1       2       1     540       0       0       0
Si             23.20399721      1.81663782      1.02475250      14       2       0       0       0     541       0       0       0
Si             27.14315001      0.02500312      1.09771711      14       3       0       0       0     542       0       0       0
O              25.70266868     -0.49858414      1.58084879       8       8       0       0       0     543       0       0       0
O              26.31123392     -2.40709728     -1.24354281       8       9       0       0       0     544       0       0       0
O              23.35290852      1.27616673     -0.30828675       8       7       0       0       0     545       0       0       0
Si             25.05353505     -1.74055808     -2.22772540      14       1       0       0       1     546       0       0       0
O              24.03461733     -1.13717108     -1.26964238       8       4       0       0       1     547       0       0       0
O              24.55776217      2.37523282      1.75128814       8       5       0       0       1     548       0       0       0
O              26.97421534      0.65151202     -0.31383316       8       6       0       0       1     549       0       0       0
O              22.07190690     -2.32124939      1.11992370       8       9       0       0       0     550       0       0       0
Si             20.92720603     -1.84685293      0.09825551      14       1       0       0       1     551       0       0       0
O              22.69022499      0.60722331      2.04845574       8       6       0       0       1     552       0       0       0
Si             27.44845402      1.73203874     -1.28100683      14       2       0       1       0     553       0       0       0
O              28.77843604      2.44408150     -0.72855229       8       5       0       1       1     554       0       0       0
Si             22.83020339     -0.06035354     -1.20464151      14       3       0       0       0     555       0       0       0
O              21.47955041     -0.62758779     -0.80802110       8       8       0       0       0     556       0       0       0
O              27.61817269      1.10097452      2.10016948       8       7       1       1       0     557       0       0       0
O              28.20962436     -1.20192775      1.05305826       8       4       1       1       1     558       0       0       0
Si             23.25961284      7.11621623      1.10694566      14       2       0       0       0     559       0       0       0
Si             27.02319893      5.34323622      1.21778350      14       3       0       0       0     560       0       0       0
O              25.61352604      4.82486840      1.63728481       8       8       0       0       0     561       0       0       0
O              26.32356940      2.93710104     -1.34601342       8       9       0       0       0     562       0       0       0
O              23.48737371      6.50947326     -0.36473758       8       7       0       0       0     563       0       0       0
Si             25.05387806      3.62259687     -2.33865326      14       1       0       0       1     564       0       0       0
O              24.01745565      4.15244628     -1.33870400       8       4       0       0       1     565       0       0       0
O              24.63113325      7.72022834      1.69486634       8       5       0       0       1     566       0       0       0
O              26.86718800      5.86892695     -0.33703528       8       6       0       0       1     567       0       0       0
O              22.21081779      2.96511251      1.17645498       8       9       0       0       0     568       0       0       0
Si             20.98833664      3.57653981      0.20164098      14       1       0       0       1     569       0       0       0
O              22.58649373      5.83993389      2.07893040       8       6       0       0       1     570       0       0       0
Si             27.39436142      7.01128041     -1.26866764      14       2       0       1       0     571       0       0       0
O              28.85752805      7.61211196     -0.83959781       8       5       0       1       1     572       0       0       0
Si             22.97131880      5.32401622     -1.19757281      14       3       0       0       0     573       0       0       0
O              21.54688097      4.64859827     -0.81878008       8       8       0       0       0     574       0       0       0
O              27.65200854      6.54520937      2.05122859       8       7       1       1       0     575       0       0       0
O              28.19055448      4.05659694      1.06082273       8       4       1       1       1     576       0       0       0
Si             23.16040680     12.38341498      1.13883679      14       2       0       0       0     577       0       0       0
Si             27.10923777     10.64234660      1.13530561      14       3       0       0       0     578       0       0       0
O              25.74328897      9.99548227      1.64149219       8       8       0       0       0     579       0       0       0
O              26.29023887      8.25025265     -1.22369216       8       9       0       0       0     580       0       0       0
O              23.49654164     11.75270753     -0.30239695       8       7       0       0       0     581       0       0       0
Si             25.09036481      8.97788245     -2.16881855      14       1       0       0       1     582       0       0       0
O              24.03242584      9.48791915     -1.32733601       8       4       0       0       1     583       0       0       0
O              24.55169215     12.99362453      1.69095054       8       5       0       0       1     584       0       0       0
O              26.93209850     11.14708476     -0.36996062       8       6       0       0       1     585       0       0       0
O              22.11340752      8.27411000      1.06262347       8       9       0       0       0     586       0       0       0
Si             21.02100411      8.87138506      0.17849778      14       1       0       0       1     587       0       0       0
O              22.65615273     11.32797973      2.00573267       8       6       0       0       1     588       0       0       0
Si             27.34021438     12.38892351     -1.35397673      14       2       0       1       0     589       0       0       0
O              28.83912139     13.04091014     -0.74696368       8       5       0       1       1     590       0       0       0
Si             22.86696728     10.69008168     -1.39402666      14       3       0       0       0     591       0       0       0
O              21.56369914     10.00950912     -0.76002374       8       8       0       0       0     592       0       0       0
O              27.66616348     11.84964426      2.16687935       8       7       1       1       0     593       0       0       0
O              28.11913390      9.52290458      1.16630903       8       4       1       1       1     594       0       0       0
Si             23.14840527     17.68106249      1.04008552      14       2       0       0       0     595       0       0       0
Si             27.03692136     16.00677907      1.15457147      14       3       0       0       0     596       0       0       0
O              25.63744668     15.30605291      1.61614791       8       8       0       0       0     597       0       0       0
O              26.29516946     13.52607560     -1.39482688       8       9       0       0       0     598       0       0       0
O              23.33070155     17.08892132     -0.38291582       8       7       0       0       0     599       0       0       0
Si             25.07354138     14.20567744     -2.21135346      14       1       0       0       1     600       0       0       0
O              24.01787448     14.87926925     -1.30826514       8       4       0       0       1     601       0       0       0
O              24.65372691     18.32758924      1.57986186       8       5       0       0       1     602       0       0       0
O              26.79344455     16.48530962     -0.28130042       8       6       0       0       1     603       0       0       0
O              22.07420052     13.69693675      1.09729144       8       9       0       0       0     604       0       0       0
Si             20.90098956     14.21325686      0.24475891      14       1       0       0       1     605       0       0       0
O              22.76891247     16.53409027      2.00507001       8       6       0       0       1     606       0       0       0
Si             27.43973863     17.85435315     -1.38005042      14       2       0       1       0     607       0       0       0
O              28.79644141     18.43413012     -0.74592723       8       5       0       1       1     608       0       0       0
Si             22.86614248     15.93757581     -1.30588814      14       3       0       0       0     609       0       0       0
O              21.58402094     15.34111930     -0.84634228       8       8       0       0       0     610       0       0       0
O              27.60433005     17.19479733      1.97541013       8       7       1       1       0     611       0       0       0
O              28.25954580     14.86263207      1.02049900       8       4       1       1       1     612       0       0       0
Si             23.24679321     23.05761490      1.03152891      14       2       0       0       0     613       0       0       0
Si             27.09273803     21.33138190      1.05854122      14       3       0       0       0     614       0       0       0
O              25.63418504     20.66181393      1.68617655       8       8       0       0       0     615       0       0       0
O              26.41496165     18.89210901     -1.32306752       8       9       0       0       0     616       0       0       0
O              23.35734270     22.55780664     -0.25669733       8       7       0       0       0     617       0       0       0
Si             25.09187773     19.46794104     -2.33263939      14       1       0       0       1     618       0       0       0
O              23.96919803     20.07583384     -1.24493492       8       4       0       0       1     619       0       0       0
O              24.57461198     23.71290074      1.68677328       8       5       0       0       1     620       0       0       0
O              26.80668894     21.80184496     -0.27735616       8       6       0       0       1     621       0       0       0
O              22.17372384     18.85962914      1.01291516       8       9       0       0       0     622       0       0       0
Si             20.89721344     19.49542653      0.09323437      14       1       0       0       1     623       0       0       0
O              22.71422489     21.86615174      2.14015816       8       6       0       0       1     624       0       0       0
Si             27.44833692     23.06271390     -1.28272027      14       2       0       1       0     625       0       0       0
O              28.72781096     23.73643988     -0.73742240       8       5       0       1       1     626       0       0       0
Si             22.92542974     21.22129221     -1.39419454      14       3       0       0       0     627       0       0       0
O              21.56683581     20.70104398     -0.72067060       8       8       0       0       0     628       0       0       0
O              27.64491260     22.53319522      2.13707211       8       7       1       1       0     629       0       0       0
O              28.19428002     20.05125325      1.20099767       8       4       1       1       1     630       0       0       0
Si             23.16467930    -24.79489763      1.10998437      14       2       0       1       0     631       0       0       0
Si             27.14345306    -26.61350714      1.08567522      14       3       0       1       0     632       0       0       0
O              25.79800470     25.98094762      1.65202339       8       8       0       0       0     633       0       0       0
O              26.30486522     24.22333375     -1.37281303       8       9       0       0       0     634       0       0       0
O              23.48181233    -25.40207677     -0.40743156       8       7       0       1       0     635       0       0       0
Si             25.06258552     24.84712668     -2.26005810      14       1       0       0       1     636       0       0       0
O              24.05001600     25.44807045     -1.28151015       8       4       0       0       1     637       0       0       0
O              24.53222160    -24.20422325      1.70321678       8       5       0       1       1     638       0       0       0
O              26.85635802    -26.13040022     -0.29815324       8       6       0       1       1     639       0       0       0
O              22.08787626     24.37458883      1.04513537       8       9       0       0       0     640       0       0       0
Si             20.87548515     24.77497129      0.21914494      14       1       0       0       1     641       0       0       0
O              22.76637196    -26.13800088      2.12658998       8       6       0       1       1     642       0       0       0
Si             27.37150657    -24.92179766     -1.28443450      14       2       0       2       0     643       0       0       0
O              28.79180397    -24.22039725     -0.75598702       8       5       0       2       1     644       0       0       0
Si             22.90406214     26.63379458     -1.36515710      14       3       0       0       0     645       0       0       0
O              21.44348236     26.08853916     -0.83983059       8       8       0       0       0     646       0       0       0
O              27.55855967    -25.51631221      2.06030564       8       7       1       2       0     647       0       0       0
O              28.09009232     25.38452082      1.13263285       8       4       1       1       1     648       0       0       0
Si             23.20716719    -19.53375168      1.06590267      14       2       0       1       0     649       0       0       0
Si             27.04944168    -21.26603729      1.10338333      14       3       0       1       0     650       0       0       0
O              25.68035534    -21.90633221      1.67059690       8       8       0       1       0     651       0       0       0
O              26.24136600    -23.60107895     -1.22743726       8       9       0       1       0     652       0       0       0
O              23.48288634    -20.09057004     -0.27716225       8       7       0       1       0     653       0       0       0
Si             25.07267408    -23.14045267     -2.17197689      14       1       0       1       1     654       0       0       0
O              24.08010344    -22.55855915     -1.21909794       8       4       0       1       1     655       0       0       0
O              24.52920383    -18.97838092      1.72575998       8       5       0       1       1     656       0       0       0
O              26.89181967    -20.67804172     -0.41203930       8       6       0       1       1     657       0       0       0
O              22.03682912    -23.58741423      1.16974990       8       9       0       1       0     658       0       0       0
Si             21.01688768    -23.01788948      0.10665332      14       1       0       1       1     659       0       0       0
O              22.62336695    -20.77704277      2.14535663       8       6       0       1       1     660       0       0       0
Si             27.40305653    -19.55546556     -1.32908112      14       2       0       2       0     661       0       0       0
O              28.83367999    -18.93088144     -0.79153223       8       5       0       2       1     662       0       0       0
Si             22.83522518    -21.30472912     -1.27555248      14       3       0       1       0     663       0       0       0
O              21.52726757    -21.92097211     -0.83699885       8       8       0       1       0     664       0       0       0
O              27.56059230    -20.20236828      2.14895685       8       7       1       2       0     665       0       0       0
O              28.22669377    -22.46309064      1.07798658       8       4       1       2       1     666       0       0       0
Si             23.25016588    -14.29177031      1.18139004      14       2       0       1       0     667       0       0       0
Si             27.03244158    -15.91235653      1.16550290      14       3       0       1       0     668       0       0       0
O              25.78133633    -16.49761341      1.71189895       8       8       0       1       0     669       0       0       0
O              26.22370796    -18.30478761     -1.31189282       8       9       0       1       0     670       0       0       0
O              23.49742337    -14.79684351     -0.35931617       8       7       0       1       0     671       0       0       0
Si             25.16313807    -17.72191737     -2.21770723      14       1       0       1       1     672       0       0       0
O              23.96360577    -17.19801818     -1.28010261       8       4       0       1       1     673       0       0       0
O              24.64552099    -13.56516706      1.61823296       8       5       0       1       1     674       0       0       0
O              26.92072445    -15.34355751     -0.30188963       8       6       0       1       1     675       0       0       0
O              22.22532529    -18.41537205      1.18909587       8       9       0       1       0     676       0       0       0
Si             20.91508143    -17.79818317      0.23642320      14       1       0       1       1     677       0       0       0
O              22.58620450    -15.43204429      2.11366115       8       6       0       1       1     678       0       0       0
Si             27.48425625    -14.23656328     -1.23245489      14       2       0       2       0     679       0       0       0
O              28.85724391    -13.61133715     -0.84773018       8       5       0       2       1     680       0       0       0
Si             22.95100906    -16.07723452     -1.36183227      14       3       0       1       0     681       0       0       0
O              21.52038789    -16.53830915     -0.78072532       8       8       0       1       0     682       0       0       0
O              27.65866010    -14.85380208      2.11424284       8       7       1       2       0     683       0       0       0
O              28.23096047    -17.19372308      1.19744709       8       4       1       2       1     684       0       0       0
Si             23.18055803     -8.89069481      1.16996898      14       2       0       1       0     685       0       0       0
Si             27.15119628    -10.68730265      1.11919200      14       3       0       1       0     686       0       0       0
O              25.75906788    -11.26512303      1.75527969       8       8       0       1       0     687       0       0       0
O              26.34239518    -13.08189255     -1.34287822       8       9       0       1       0     688       0       0       0
O              23.40574915     -9.55197841     -0.32431315       8       7       0       1       0     689       0       0       0
Si             25.19081891    -12.42297135     -2.34648593      14       1       0       1       1     690       0       0       0
O              23.88358425    -11.82931820     -1.38254614       8       4       0       1       1     691       0       0       0
O              24.64875093     -8.36704230      1.63003708       8       5       0       1       1     692       0       0       0
O              26.95662197    -10.12612858     -0.37525331       8       6       0       1       1     693       0       0       0
O              22.07572320    -13.06067833      1.09373654       8       9       0       1       0     694       0       0       0
Si             20.88361786    -12.36707240      0.23274742      14       1       0       1       1     695       0       0       0
O              22.62950256    -10.14163273      2.02350833       8       6       0       1       1     696       0       0       0
Si             27.35789508     -8.87044558     -1.38026476      14       2       0       2       0     697       0       0       0
O              28.86227019     -8.28251431     -0.74484714       8       5       0       2       1     698       0       0       0
Si             22.91387217    -10.68836399     -1.20056668      14       3       0       1       0     699       0       0       0
O              21.47439246    -11.27165565     -0.76896947       8       8       0       1       0     700       0       0       0
O              27.54493979     -9.37378202      2.01897743       8       7       1       2       0     701       0       0       0
O              28.13288345    -11.83782280      1.15704074       8       4       1       2       1     702       0       0       0
Si             23.22835706     -3.52180310      1.14651012      14       2       0       1       0     703       0       0       0
Si             27.06195948     -5.33037598      1.02490308      14       3       0       1       0     704       0       0       0
O              25.70041505     -5.88392509      1.76831034       8       8       0       1       0     705       0       0       0
O              26.30857636     -7.66837960     -1.32063434       8       9       0       1       0     706       0       0       0
O              23.49100304     -4.20804576     -0.39184690       8       7       0       1       0     707       0       0       0
Si             25.13662694     -7.13776257     -2.21377141      14       1       0       1       1     708       0       0       0
O              23.95210144     -6.52788079     -1.36130589       8       4       0       1       1     709       0       0       0
O              24.62876592     -3.04205392      1.66377811       8       5       0       1       1     710       0       0       0
O              26.82163506     -4.77117825     -0.32228706       8       6       0       1       1     711       0       0       0
O              22.12130922     -7.74380143      1.02614884       8       9       0       1       0     712       0       0       0
Si             20.92241756     -7.03519815      0.16684829      14       1       0       1       1     713       0       0       0
O              22.71081688     -4.77549810      2.11458036       8       6       0       1       1     714       0       0       0
Si             27.47294782     -3.60699147     -1.20469347      14       2       0       2       0     715       0       0       0
O              28.77520229     -2.98598197     -0.84484326       8       5       0       2       1     716       0       0       0
Si             22.96213037     -5.29191119     -1.23030123      14       3       0       1       0     717       0       0       0
O              21.50844327     -6.00013069     -0.66129196       8       8       0       1       0     718       0       0       0
O              27.64851472     -4.20093458      2.16470278       8       7       1       2       0     719       0       0       0
O              28.10122276     -6.48110620      1.18933876       8       4       1       2       1     720       0       0       0
Si             31.57596492      1.82662599      1.04924369      14       2       0       0       0     721       0       0       0
Si             35.56277994      0.02355086      1.20059745      14       3       0       0       0     722       0       0       0
O              34.11382458     -0.50889137      1.60164175       8       8       0       0       0     723       0       0       0
O              34.75689363     -2.42714571     -1.34631478       8       9       0       0       0     724       0       0       0
O              31.81010186      1.23215313     -0.38041889       8       7       0       0       0     725       0       0       0
Si             33.53475420     -1.69272264     -2.29082974      14       1       0       0       1     726       0       0       0
O              32.43532741     -1.26129565     -1.34745720       8       4       0       0       1     727       0       0       0
O              32.89111876      2.30965855      1.67966964       8       5       0       0       1     728       0       0       0
O              35.27577887      0.59406800     -0.43694787       8       6       0       0       1     729       0       0       0
O              30.60019834     -2.39208189      1.07634637       8       9       0       0       0     730       0       0       0
Si             29.34682220     -1.73804685      0.19004744      14       1       0       0       1     731       0       0       0
O              31.00812469      0.49624461      1.99438909       8       6       0       0       1     732       0       0       0
Si             35.85084938      1.70832750     -1.34819949      14       2       0       1       0     733       0       0       0
O              37.11977808      2.40415341     -0.78656535       8       5       0       1       1     734       0       0       0
Si             31.32916986      0.09412374     -1.38189031      14       3       0       0       0     735       0       0       0
O              29.90888809     -0.49882101     -0.68758300       8       8       0       0       0     736       0       0       0
O              36.09082188      1.27165698      2.02021871       8       7       1       1       0     737       0       0       0
O              36.58350767     -1.26589492      1.08389264       8       4       1       1       1     738       0       0       0
Si             31.65696421      7.03388446      1.10447598      14       2       0       0       0     739       0       0       0
Si             35.38250412      5.24297696      1.08434054      14       3       0       0       0     740       0       0       0
O              34.19059680      4.75230097      1.59661693       8       8       0       0       0     741       0       0       0
O              34.61908162      3.04411549     -1.23970617       8       9       0       0       0     742       0       0       0
O              31.77358764      6.52632597     -0.31637230       8       7       0       0       0     743       0       0       0
Si             33.62437356      3.60869239     -2.29340853      14       1       0       0       1     744       0       0       0
O              32.26853527      4.07477118     -1.39993186       8       4       0       0       1     745       0       0       0
O              33.04210337      7.76205201      1.59576833       8       5       0       0       1     746       0       0       0
O              35.33336674      5.85074004     -0.39736399       8       6       0       0       1     747       0       0       0
O              30.49961199      2.87857108      1.01158113       8       9       0       0       0     748       0       0       0
Si             29.29439923      3.53664754      0.10243304      14       1       0       0       1     749       0       0       0
O              31.12976944      5.92427597      2.06254512       8       6       0       0       1     750       0       0       0
Si             35.77683390      7.19884707     -1.20960120      14       2       0       1       0     751       0       0       0
O              37.09631798      7.60778052     -0.76230737       8       5       0       1       1     752       0       0       0
Si             31.37459445      5.28330762     -1.25496898      14       3       0       0       0     753       0       0       0
O              29.98002989      4.66481433     -0.73506015       8       8       0       0       0     754       0       0       0
O              35.99936365      6.49805551      2.16278805       8       7       1       1       0     755       0       0       0
O              36.62998083      4.09488264      1.03741468       8       4       1       1       1     756       0       0       0
Si             31.57308422     12.43588275      1.17121443      14       2       0       0       0     757       0       0       0
Si             35.51761052     10.62669934      1.09347475      14       3       0       0       0     758       0       0       0
O              34.07232691     10.15796253      1.76931174       8       8       0       0       0     759       0       0       0
O              34.70468711      8.35413181     -1.34785570       8       9       0       0       0     760       0       0       0
O              31.74707208     11.92069386     -0.25689464       8       7       0       0       0     761       0       0       0
Si             33.62814086      8.80677353     -2.17471253      14       1       0       0       1     762       0       0       0
O              32.32559058      9.48061846     -1.29887072       8       4       0       0       1     763       0       0       0
O              32.89824928     12.96091027      1.57655891       8       5       0       0       1     764       0       0       0
O              35.18998750     11.31542681     -0.41554486       8       6       0       0       1     765       0       0       0
O              30.41446428      8.19914614      1.01307479       8       9       0       0       0     766       0       0       0
Si             29.38643933      8.80202218      0.23189602      14       1       0       0       1     767       0       0       0
O              31.08213120     11.21699051      1.98232943       8       6       0       0       1     768       0       0       0
Si             35.84017717     12.41364146     -1.23950519      14       2       0       1       0     769       0       0       0
O              37.17086624     12.97339778     -0.84353079       8       5       0       1       1     770       0       0       0
Si             31.34738886     10.57626061     -1.34575582      14       3       0       0       0     771       0       0       0
O              29.93998874      9.99606300     -0.78960449       8       8       0       0       0     772       0       0       0
O              36.05662432     11.81208054      2.08046051       8       7       1       1       0     773       0       0       0
O              36.53993690      9.55576299      1.11418469       8       4       1       1       1     774       0       0       0
Si             31.61892274     17.70322799      1.20087360      14       2       0       0       0     775       0       0       0
Si             35.41687640     16.03502808      1.15086890      14       3       0       0       0     776       0       0       0
O              34.05709278     15.35247866      1.60336173       8       8       0       0       0     777       0       0       0
O              34.63805963     13.53211431     -1.41301909       8       9       0       0       0     778       0       0       0
O              31.79088918     17.08022604     -0.30025477       8       7       0       0       0     779       0       0       0
Si             33.48342683     14.13993853     -2.23967306      14       1       0       0       1     780       0       0       0
O              32.34813184     14.78654527     -1.30269222       8       4       0       0       1     781       0       0       0
O              33.07260184     18.39765798      1.65389033       8       5       0       0       1     782       0       0       0
O              35.22518828     16.61606608     -0.33468826       8       6       0       0       1     783       0       0       0
O              30.60814401     13.59514415      1.08289871       8       9       0       0       0     784       0       0       0
Si             29.37051142     14.20522520      0.16254215      14       1       0       0       1     785       0       0       0
O              31.04424973     16.63311785      2.08608714       8       6       0       0       1     786       0       0       0
Si             35.85275840     17.82942726     -1.29832817      14       2       0       1       0     787       0       0       0
O              37.19336760     18.31046110     -0.83730048       8       5       0       1       1     788       0       0       0
Si             31.23943127     15.88642501     -1.20143913      14       3       0       0       0     789       0       0       0
O              29.95940704     15.31022710     -0.67809363       8       8       0       0       0     790       0       0       0
O              35.98764232     17.22457707      2.04335634       8       7       1       1       0     791       0       0       0
O              36.53014256     14.74865488      1.11107741       8       4       1       1       1     792       0       0       0
Si             31.55542400     23.08656785      1.20996828      14       2       0       0       0     793       0       0       0
Si             35.52039131     21.25250227      1.10298648      14       3       0       0       0     794       0       0       0
O              34.04645881     20.79878751      1.76002386       8       8       0       0       0     795       0       0       0
O              34.76519431     18.95386578     -1.41232079       8       9       0       0       0     796       0       0       0
O              31.76821950     22.51097163     -0.30349137       8       7       0       0       0     797       0       0       0
Si             33.54396735     19.49553726     -2.18738550      14       1       0       0       1     798       0       0       0
O              32.36433174     20.12505228     -1.37965103       8       4       0       0       1     799       0       0       0
O              33.04862865     23.57539450      1.60364268       8       5       0       0       1     800       0       0       0
O              35.33880170     21.85428611     -0.33443562       8       6       0       0       1     801       0       0       0
O              30.44852271     18.92294286      1.09771078       8       9       0       0       0     802       0       0       0
Si             29.36031988     19.49471152      0.19071219      14       1       0       0       1     803       0       0       0
O              31.14887254     21.81601754      1.99399343       8       6       0       0       1     804       0       0       0
Si             35.77048494     23.09462692     -1.34124966      14       2       0       1       0     805       0       0       0
O              37.25382392     23.73269759     -0.83518669       8       5       0       1       1     806       0       0       0
Si             31.24482283     21.40261570     -1.37519773      14       3       0       0       0     807       0       0       0
O              29.92260213     20.72882289     -0.70584563       8       8       0       0       0     808       0       0       0
O              36.03754537     22.52429037      1.97281960       8       7       1       1       0     809       0       0       0
O              36.50428148     20.08082442      1.10970250       8       4       1       1       1     810       0       0       0
Si             31.49212212    -24.93661345      1.08312601      14       2       0       1       0     811       0       0       0
Si             35.46426233     26.56899436      1.19709287      14       3       0       0       0     812       0       0       0
O              34.18042516     26.10436109      1.68007710       8       8       0       0       0     813       0       0       0
O              34.67214433     24.32860815     -1.21890850       8       9       0       0       0     814       0       0       0
O              31.86511060    -25.51530643     -0.34043610       8       7       0       1       0     815       0       0       0
Si             33.47009034     24.84915743     -2.23436726      14       1       0       0       1     816       0       0       0
O              32.33513801     25.47000530     -1.35139671       8       4       0       0       1     817       0       0       0
O              32.90647691    -24.35665788      1.62153431       8       5       0       1       1     818       0       0       0
O              35.30447858    -26.00081514     -0.36717734       8       6       0       1       1     819       0       0       0
O              30.55150735     24.30571271      1.04720751       8       9       0       0       0     820       0       0       0
Si             29.38094966     24.88948969      0.08757887      14       1       0       0       1     821       0       0       0
O              31.01991074    -26.01764732      2.05254946       8       6       0       1       1     822       0       0       0
Si             35.81973352    -24.76900873     -1.37551366      14       2       0       2       0     823       0       0       0
O              37.11490976    -24.26923537     -0.70449700       8       5       0       2       1     824       0       0       0
Si             31.27395448    -26.56973677     -1.33711827      14       3       0       1       0     825       0       0       0
O              29.87522032     26.04163067     -0.75104883       8       8       0       0       0     826       0       0       0
O              35.97211621    -25.42028257      2.00665587       8       7       1       2       0     827       0       0       0
O              36.47017129     25.53622867      1.13793774       8       4       1       1       1     828       0       0       0
Si             31.61569794    -19.50509118      1.18955864      14       2       0       1       0     829       0       0       0
Si             35.52238086    -21.35347988      1.06747477      14       3       0       1       0     830       0       0       0
O              34.11544652    -21.87483796      1.57263106       8       8       0       1       0     831       0       0       0
O              34.73320126    -23.70473810     -1.36823719       8       9       0       1       0     832       0       0       0
O              31.74314988    -20.15817910     -0.32712401       8       7       0       1       0     833       0       0       0
Si             33.46503244    -23.06822967     -2.28462345      14       1       0       1       1     834       0       0       0
O              32.36958245    -22.57347495     -1.25117682       8       4       0       1       1     835       0       0       0
O              33.05262392    -18.95584248      1.74161716       8       5       0       1       1     836       0       0       0
O              35.19791510    -20.64023807     -0.37830316       8       6       0       1       1     837       0       0       0
O              30.50421365    -23.69821892      1.07483934       8       9       0       1       0     838       0       0       0
Si             29.32903666    -23.15247816      0.07671618      14       1       0       1       1     839       0       0       0
O              31.06868937    -20.73559714      2.07747496       8       6       0       1       1     840       0       0       0
Si             35.71356365    -19.46973978     -1.36594296      14       2       0       2       0     841       0       0       0
O              37.13793838    -18.88465672     -0.82658526       8       5       0       2       1     842       0       0       0
Si             31.24140074    -21.24853537     -1.38612945      14       3       0       1       0     843       0       0       0
O              29.96974129    -21.99653200     -0.81326440       8       8       0       1       0     844       0       0       0
O              36.06112107    -20.16924437      1.99872537       8       7       1       2       0     845       0       0       0
O              36.54396816    -22.40556444      1.00887051       8       4       1       2       1     846       0       0       0
Si             31.61077171    -14.28546069      1.21987990      14       2       0       1       0     847       0       0       0
Si             35.41001066    -15.98755411      1.15505431      14       3       0       1       0     848       0       0       0
O              34.09465144    -16.51229787      1.64846814       8       8       0       1       0     849       0       0       0
O              34.66292813    -18.42203498     -1.40592095       8       9       0       1       0     850       0       0       0
O              31.88109193    -14.76164659     -0.44320221       8       7       0       1       0     851       0       0       0
Si             33.56946275    -17.76413059     -2.18579633      14       1       0       1       1     852       0       0       0
O              32.39344225    -17.23815760     -1.33386209       8       4       0       1       1     853       0       0       0
O              32.90250763    -13.54894025      1.59272137       8       5       0       1       1     854       0       0       0
O              35.25524199    -15.46454345     -0.41924445       8       6       0       1       1     855       0       0       0
O              30.60896038    -18.31866719      1.13781134       8       9       0       1       0     856       0       0       0
Si             29.43291003    -17.71358130      0.15885049      14       1       0       1       1     857       0       0       0
O              31.09824866    -15.38830635      1.98690087       8       6       0       1       1     858       0       0       0
Si             35.71160647    -14.18362903     -1.31839514      14       2       0       2       0     859       0       0       0
O              37.19727562    -13.70087318     -0.84017000       8       5       0       2       1     860       0       0       0
Si             31.30215708    -16.02396461     -1.34259185      14       3       0       1       0     861       0       0       0
O              29.89456907    -16.55291065     -0.74850522       8       8       0       1       0     862       0       0       0
O              35.98666212    -14.81503568      2.05903507       8       7       1       2       0     863       0       0       0
O              36.57849428    -17.14422638      1.15549936       8       4       1       2       1     864       0       0       0
Si             31.65936259     -8.81048241      1.05543948      14       2       0       1       0     865       0       0       0
Si             35.53114175    -10.71900359      1.12540226      14       3       0       1       0     866       0       0       0
O              34.10860089    -11.16112455      1.59937473       8       8       0       1       0     867       0       0       0
O              34.78847954    -13.04131859     -1.38054413       8       9       0       1       0     868       0       0       0
O              31.72210974     -9.37697956     -0.37067155       8       7       0       1       0     869       0       0       0
Si             33.63341238    -12.40974621     -2.31761430      14       1       0       1       1     870       0       0       0
O              32.27734635    -11.79047249     -1.30595254       8       4       0       1       1     871       0       0       0
O              32.91745648     -8.29318717      1.66211569       8       5       0       1       1     872       0       0       0
O              35.24567495    -10.06608798     -0.26423194       8       6       0       1       1     873       0       0       0
O              30.43314760    -13.05243498      1.13577142       8       9       0       1       0     874       0       0       0
Si             29.36798427    -12.50734973      0.09403211      14       1       0       1       1     875       0       0       0
O              30.97883004    -10.09555254      2.08896058       8       6       0       1       1     876       0       0       0
Si             35.73818628     -8.93787804     -1.31273123      14       2       0       2       0     877       0       0       0
O              37.08165759     -8.38956920     -0.74899477       8       5       0       2       1     878       0       0       0
Si             31.32972852    -10.70207064     -1.26563597      14       3       0       1       0     879       0       0       0
O              29.93452713    -11.30566388     -0.81197754       8       8       0       1       0     880       0       0       0
O              35.97677074     -9.39049836      2.03944611       8       7       1       2       0     881       0       0       0
O              36.49627918    -11.93224428      1.09784840       8       4       1       2       1     882       0       0       0
Si             31.63876958     -3.61557022      1.14635000      14       2       0       1       0     883       0       0       0
Si             35.44899340     -5.28771786      1.08050484      14       3       0       1       0     884       0       0       0
O              34.15178902     -5.82940991      1.67162401       8       8       0       1       0     885       0       0       0
O              34.70111890     -7.73913807     -1.36270883       8       9       0       1       0     886       0       0       0
O              31.81429563     -4.23282019     -0.33817330       8       7       0       1       0     887       0       0       0
Si             33.44084125     -7.01259778     -2.27832894      14       1       0       1       1     888       0       0       0
O              32.39559385     -6.48352346     -1.25058586       8       4       0       1       1     889       0       0       0
O              32.94768867     -3.00548571      1.70825188       8       5       0       1       1     890       0       0       0
O              35.20696542     -4.74931823     -0.28942571       8       6       0       1       1     891       0       0       0
O              30.49842250     -7.66957858      1.03262165       8       9       0       1       0     892       0       0       0
Si             29.31974368     -7.12992496      0.11044240      14       1       0       1       1     893       0       0       0
O              31.06318934     -4.66889536      2.07605817       8       6       0       1       1     894       0       0       0
Si             35.87715817     -3.50711482     -1.38451742      14       2       0       2       0     895       0       0       0
O              37.14232587     -2.87675011     -0.74150946       8       5       0       2       1     896       0       0       0
Si             31.34390668     -5.35831090     -1.32959768      14       3       0       1       0     897       0       0       0
O              29.91138942     -5.82582027     -0.78423030       8       8       0       1       0     898       0       0       0
O              36.07888084     -4.14126031      2.01350985       8       7       1       2       0     899       0       0       0
O              36.52963090     -6.49150306      1.00270978       8       4       1       2       1     900       0       0       0
Si             39.92225812      1.82350411      1.16080808      14       2       0       0       0     901       0       0       0
O              40.20422889      1.25489899     -0.29317482       8       7       0       0       0     905       0       0       0
Si             41.89420149     -1.73266323     -2.21729327      14       1       0       0       1     906       0       0       0
O              40.84954368     -1.10327374     -1.41382241       8       4       0       0       1     907       0       0       0
O              41.27301567      2.45682533      1.62786382       8       5       0       0       1     908       0       0       0
O              38.80165101     -2.41532187      1.17840337       8       9       0       0       0     910       0       0       0
Si             37.74127111     -1.73297738      0.13844726      14       1       0       0       1     911       0       0       0
O              39.38911100      0.63420763      2.14966599       8       6       0       0       1     912       0       0       0
Si             39.74779093     -0.04960191     -1.34250719      14       3       0       0       0     915       0       0       0
O              38.27498892     -0.61683361     -0.80092478       8       8       0       0       0     916       0       0       0
Si             39.97212300      7.11378649      1.20430999      14       2       0       0       0     919       0       0       0
O              40.12139176      6.61856386     -0.44337111       8       7       0       0       0     923       0       0       0
O              40.74593581      4.04415280     -1.37539291       8       4       0       0       1     925       0       0       0
O              41.35111891      7.75290010      1.74574801       8       5       0       0       1     926       0       0       0
O              38.90095535      2.87817502      1.00532569       8       9       0       0       0     928       0       0       0
Si             37.70885413      3.63041034      0.16206566      14       1       0       0       1     929       0       0       0
O              39.40679570      5.89415796      2.04137867       8       6       0       0       1     930       0       0       0
Si             39.57796112      5.24389046     -1.31759712      14       3       0       0       0     933       0       0       0
O              38.24285246      4.69285315     -0.81234835       8       8       0       0       0     934       0       0       0
Si             39.89503097     12.37017231      1.20218894      14       2       0       0       0     937       0       0       0
Si            -39.93854617     10.64324274      1.15227557      14       3       1       0       0     938       0       0       0
O              40.12166891     11.76912253     -0.42685641       8       7       0       0       0     941       0       0       0
Si             41.83556937      8.78577784     -2.32302265      14       1       0       0       1     942       0       0       0
O              40.73037694      9.50676638     -1.36530827       8       4       0       0       1     943       0       0       0
O              41.42073559     13.11377926      1.72365098       8       5       0       0       1     944       0       0       0
O             -40.25461236     11.29881214     -0.31300708       8       6       1       0       1     945       0       0       0
O              38.94898703      8.34908749      1.11354286       8       9       0       0       0     946       0       0       0
Si             37.68509556      8.91596284      0.12444899      14       1       0       0       1     947       0       0       0
O              39.39639560     11.30087963      2.07111296       8       6       0       0       1     948       0       0       0
Si            -39.76946731     12.34122623     -1.21084688      14       2       1       1       0     949       0       0       0
O             -38.30488986     12.97089009     -0.81205732       8       5       1       1       1     950       0       0       0
Si             39.64902672     10.69612751     -1.24995559      14       3       0       0       0     951       0       0       0
O              38.31649066     10.14571459     -0.67468458       8       8       0       0       0     952       0       0       0
O             -39.53883362     11.88341296      2.10143045       8       7       2       1       0     953       0       0       0
Si             39.87646956     17.84900837      1.04922661      14       2       0       0       0     955       0       0       0
Si            -40.05039511     16.06624775      1.12640001      14       3       1       0       0     956       0       0       0
O             -41.33098528     15.31817668      1.59713096       8       8       1       0       0     957       0       0       0
O             -40.77000184     13.56326627     -1.25327805       8       9       1       0       0     958       0       0       0
O              40.15715277     17.13205420     -0.34792237       8       7       0       0       0     959       0       0       0
Si            -41.90734645     14.25707972     -2.28850208      14       1       1       0       1     960       0       0       0
O              40.84922531     14.86600309     -1.30479270       8       4       0       0       1     961       0       0       0
O              41.42824044     18.36214224      1.62319680       8       5       0       0       1     962       0       0       0
O             -40.15302287     16.51023221     -0.33575723       8       6       1       0       1     963       0       0       0
O              38.97908331     13.60768424      1.04819692       8       9       0       0       0     964       0       0       0
Si             37.63498170     14.13696918      0.11721408      14       1       0       0       1     965       0       0       0
O              39.36708399     16.65378233      2.07124043       8       6       0       0       1     966       0       0       0
Si            -39.67800917     17.81465287     -1.22971392      14       2       1       1       0     967       0       0       0
O             -38.31478819     18.34852053     -0.72016701       8       5       1       1       1     968       0       0       0
Si             39.76026084     15.97628702     -1.37795893      14       3       0       0       0     969       0       0       0
O              38.28253466     15.32040934     -0.82230750       8       8       0       0       0     970       0       0       0
O             -39.44508748     17.14075067      2.04559933       8       7       2       1       0     971       0       0       0
O             -38.81637952     14.71825499      1.06568538       8       4       2       1       1     972       0       0       0
Si             39.96248139     23.14673830      1.20163753      14       2       0       0       0     973       0       0       0
Si            -39.89774751     21.32307276      1.20165171      14       3       1       0       0     974       0       0       0
O             -41.41062679     20.67482041      1.57703414       8       8       1       0       0     975       0       0       0
O             -40.79577127     18.89803083     -1.37631366       8       9       1       0       0     976       0       0       0
O              40.15861979     22.45761924     -0.37498073       8       7       0       0       0     977       0       0       0
Si            -41.90941599     19.49969590     -2.31151187      14       1       1       0       1     978       0       0       0
O              40.69455895     20.10019346     -1.23468078       8       4       0       0       1     979       0       0       0
O              41.39562897     23.63451155      1.60898903       8       5       0       0       1     980       0       0       0
O             -40.26438885     21.93473828     -0.30585408       8       6       1       0       1     981       0       0       0
O              38.81425165     18.86020499      1.19315556       8       9       0       0       0     982       0       0       0
Si             37.80107431     19.50354341      0.16894727      14       1       0       0       1     983       0       0       0
O              39.40039104     21.94788547      2.13038307       8       6       0       0       1     984       0       0       0
Si            -39.69560962     23.17883228     -1.31115320      14       2       1       1       0     985       0       0       0
O             -38.27031325     23.73554828     -0.68643131       8       5       1       1       1     986       0       0       0
Si             39.67797082     21.26226126     -1.24445738      14       3       0       0       0     987       0       0       0
O              38.19982645     20.67594802     -0.65720043       8       8       0       0       0     988       0       0       0
O             -39.37819211     22.44360696      2.03485921       8       7       2       1       0     989       0       0       0
O             -38.80520447     20.19915194      1.15509132       8       4       2       1       1     990       0       0       0
Si             40.03565584    -24.91887421      1.20771697      14       2       0       1       0     991       0       0       0
Si            -40.01639402    -26.59098389      1.18554775      14       3       1       1       0     992       0       0       0
O             -41.35577378     26.01009888      1.73393682       8       8       1       0       0     993       0       0       0
O             -40.78264987     24.18833063     -1.22875840       8       9       1       0       0     994       0       0       0
O              40.20388097    -25.40448690     -0.39451469       8       7       0       1       0     995       0       0       0
Si            -41.82404815     24.91168571     -2.25875751      14       1       1       0       1     996       0       0       0
O              40.71287380     25.36896826     -1.39720655       8       4       0       0       1     997       0       0       0
O              41.39510246    -24.31637401      1.74652575       8       5       0       1       1     998       0       0       0
O             -40.22775407    -25.96821588     -0.38391141       8       6       1       1       1     999       0       0       0
O              38.95032705     24.21430407      1.18282141       8       9       0       0       0    1000       0       0       0
Si             37.69219796     24.78377985      0.11417518      14       1       0       0       1    1001       0       0       0
O              39.53098277    -26.03896816      2.03746037       8       6       0       1       1    1002       0       0       0
Si            -39.62246888    -24.78024738     -1.25870120      14       2       1       2       0    1003       0       0       0
O             -38.23209525    -24.35495860     -0.66457986       8       5       1       2       1    1004       0       0       0
Si             39.66670623    -26.64001756     -1.25920787      14       3       0       1       0    1005       0       0       0
O              38.22677457     26.09468878     -0.78110246       8       8       0       0       0    1006       0       0       0
O             -39.54472228    -25.51681048      1.97554790       8       7       2       2       0    1007       0       0       0
O             -38.93106583     25.35342841      1.06812258       8       4       2       1       1    1008       0       0       0
Si             39.88057692    -19.46038472      1.11862562      14       2       0       1       0    1009       0       0       0
Si            -39.93438330    -21.40704976      1.16947344      14       3       1       1       0    1010       0       0       0
O             -41.32727313    -21.94733544      1.74979504       8       8       1       1       0    1011       0       0       0
O             -40.76443886    -23.64466496     -1.41715624       8       9       1       1       0    1012       0       0       0
O              40.24650530    -20.16549584     -0.43730438       8       7       0       1       0    1013       0       0       0
Si            -41.91172715    -23.18832916     -2.29042871      14       1       1       1       1    1014       0       0       0
O              40.65850520    -22.49044923     -1.31129426       8       4       0       1       1    1015       0       0       0
O              41.28830176    -18.99163181      1.67822644       8       5       0       1       1    1016       0       0       0
O             -40.25426423    -20.79906212     -0.36759858       8       6       1       1       1    1017       0       0       0
O              38.86256435    -23.58074511      1.17849992       8       9       0       1       0    1018       0       0       0
Si             37.63749855    -23.15727623      0.21887407      14       1       0       1       1    1019       0       0       0
O              39.47404066    -20.68391784      1.99473148       8       6       0       1       1    1020       0       0       0
Si            -39.67039822    -19.45977651     -1.34976981      14       2       1       2       0    1021       0       0       0
O             -38.26112947    -18.88144724     -0.74805999       8       5       1       2       1    1022       0       0       0
Si             39.65928154    -21.37724460     -1.32429336      14       3       0       1       0    1023       0       0       0
O              38.35416176    -21.95146539     -0.68739715       8       8       0       1       0    1024       0       0       0
O             -39.46220549    -20.21897990      2.05140632       8       7       2       2       0    1025       0       0       0
O             -38.99235089    -22.41925850      1.19000495       8       4       2       2       1    1026       0       0       0
Si             40.04776801    -14.13197311      1.05988471      14       2       0       1       0    1027       0       0       0
Si            -39.94191407    -15.96647878      1.11995091      14       3       1       1       0    1028       0       0       0
O             -41.29674253    -16.64842956      1.74991880       8       8       1       1       0    1029       0       0       0
O             -40.66017808    -18.41176923     -1.40527084       8       9       1       1       0    1030       0       0       0
O              40.10312998    -14.77873635     -0.39680494       8       7       0       1       0    1031       0       0       0
Si             41.86594322    -17.81141528     -2.30623223      14       1       0       1       1    1032       0       0       0
O              40.66320190    -17.21698123     -1.22145550       8       4       0       1       1    1033       0       0       0
O              41.34567472    -13.69046879      1.66476007       8       5       0       1       1    1034       0       0       0
O             -40.12018168    -15.46292258     -0.27501704       8       6       1       1       1    1035       0       0       0
O              38.81175639    -18.35379353      1.17255752       8       9       0       1       0    1036       0       0       0
Si             37.72165835    -17.75977188      0.09419738      14       1       0       1       1    1037       0       0       0
O              39.37812592    -15.39852495      2.03752792       8       6       0       1       1    1038       0       0       0
Si            -39.61823330    -14.10996958     -1.29745797      14       2       1       2       0    1039       0       0       0
O             -38.21633724    -13.69874710     -0.67334215       8       5       1       2       1    1040       0       0       0
Si             39.58146861    -15.93433803     -1.20042626      14       3       0       1       0    1041       0       0       0
O              38.28190982    -16.51155175     -0.78173138       8       8       0       1       0    1042       0       0       0
O             -39.50513509    -14.77715191      2.12392893       8       7       2       2       0    1043       0       0       0
O             -38.80493684    -17.18446796      1.17351590       8       4       2       2       1    1044       0       0       0
Si             40.06193291     -8.80218801      1.07964327      14       2       0       1       0    1045       0       0       0
Si            -39.95036758    -10.59141196      1.06317890      14       3       1       1       0    1046       0       0       0
O             -41.26664566    -11.33172090      1.70039206       8       8       1       1       0    1047       0       0       0
O             -40.78156095    -12.92535701     -1.30559018       8       9       1       1       0    1048       0       0       0
O              40.27251765     -9.56413645     -0.42173616       8       7       0       1       0    1049       0       0       0
Si             41.90584784    -12.41859623     -2.19835231      14       1       0       1       1    1050       0       0       0
O              40.69769614    -11.83400543     -1.34639951       8       4       0       1       1    1051       0       0       0
O              41.44950873     -8.22926632      1.69200413       8       5       0       1       1    1052       0       0       0
O             -40.15610390    -10.09537371     -0.45135268       8       6       1       1       1    1053       0       0       0
O              38.80841391    -13.09974559      1.11987551       8       9       0       1       0    1054       0       0       0
Si             37.80112093    -12.41688522      0.22008887      14       1       0       1       1    1055       0       0       0
O              39.49329336    -10.16734775      2.03667408       8       6       0       1       1    1056       0       0       0
Si             39.62349823    -10.58329223     -1.31215706      14       3       0       1       0    1059       0       0       0
O              38.30471946    -11.31127974     -0.78925064       8       8       0       1       0    1060       0       0       0
O             -38.98460886    -11.78629221      1.18518941       8       4       2       2       1    1062       0       0       0
Si             39.89925652     -3.46670789      1.11201366      14       2       0       1       0    1063       0       0       0
O              40.14260565     -4.10361366     -0.40603021       8       7       0       1       0    1067       0       0       0
Si             41.87794456     -7.08591144     -2.27499962      14       1       0       1       1    1068       0       0       0
O              40.82507045     -6.48957242     -1.33748419       8       4       0       1       1    1069       0       0       0
O              41.34567633     -3.01991505      1.67024898       8       5       0       1       1    1070       0       0       0
O              38.97482861     -7.64217830      1.16997210       8       9       0       1       0    1072       0       0       0
Si             37.66771213     -7.13934031      0.08521847      14       1       0       1       1    1073       0       0       0
O              39.51826281     -4.76844154      2.05815361       8       6       0       1       1    1074       0       0       0
Si             39.69250336     -5.38794422     -1.28406166      14       3       0       1       0    1077       0       0       0
O              38.29007516     -5.84804878     -0.66116080       8       8       0       1       0    1078       0       0       0
O             -29.97624495      7.58965272     -0.82250473       8       5       1       1       1    1112       0       0       0
Si            -35.57370388     12.44397304      1.12120672      14       2       1       0       0    1117       0       0       0
Si            -31.64933219     10.72744326      1.21544691      14       3       1       0       0    1118       0       0       0
O             -32.91989713     10.13131041      1.59822315       8       8       1       0       0    1119       0       0       0
O             -32.46432809      8.33309360     -1.27832936       8       9       1       0       0    1120       0       0       0
O             -35.30045528     11.83491301     -0.33607961       8       7       1       0       0    1121       0       0       0
Si            -33.63086788      8.91896014     -2.18844011      14       1       1       0       1    1122       0       0       0
O             -34.60404171      9.40880949     -1.31087793       8       4       1       0       1    1123       0       0       0
O             -34.07632279     12.97474245      1.76704788       8       5       1       0       1    1124       0       0       0
O             -31.85048648     11.29734321     -0.36472072       8       6       1       0       1    1125       0       0       0
O             -36.01011618     11.20099175      2.05231315       8       6       1       0       1    1128       0       0       0
Si            -31.27145035     12.33446473     -1.31469269      14       2       1       1       0    1129       0       0       0
O             -29.94800837     13.00426621     -0.84513175       8       5       1       1       1    1130       0       0       0
Si            -35.82020384     10.59910702     -1.30709825      14       3       1       0       0    1131       0       0       0
O             -37.19511934     10.10060188     -0.80923916       8       8       1       0       0    1132       0       0       0
O             -31.12034193     11.89040943      2.03472434       8       7       2       1       0    1133       0       0       0
O             -30.49473637      9.38935997      1.05827560       8       4       2       1       1    1134       0       0       0
Si            -35.44174831     17.67817548      1.15743300      14       2       1       0       0    1135       0       0       0
Si            -31.67497917     16.06683903      1.11142812      14       3       1       0       0    1136       0       0       0
O             -32.88334071     15.37868868      1.69111913       8       8       1       0       0    1137       0       0       0
O             -32.31613331     13.63134696     -1.24971436       8       9       1       0       0    1138       0       0       0
O             -35.19989178     17.26966391     -0.28609687       8       7       1       0       0    1139       0       0       0
Si            -33.48481500     14.26120815     -2.22967535      14       1       1       0       1    1140       0       0       0
O             -34.78752651     14.75923649     -1.33031076       8       4       1       0       1    1141       0       0       0
O             -34.09638953     18.26364762      1.65115137       8       5       1       0       1    1142       0       0       0
O             -31.82002715     16.57750790     -0.43842585       8       6       1       0       1    1143       0       0       0
O             -36.52414580     13.59189135      1.09524883       8       9       1       0       0    1144       0       0       0
Si            -37.73728851     14.11373442      0.15401223      14       1       1       0       1    1145       0       0       0
O             -35.96019194     16.62453356      2.00784535       8       6       1       0       1    1146       0       0       0
Si            -31.25269559     17.73872850     -1.23625208      14       2       1       1       0    1147       0       0       0
O             -29.86140619     18.36074199     -0.70949157       8       5       1       1       1    1148       0       0       0
Si            -35.75754135     16.03672058     -1.20171545      14       3       1       0       0    1149       0       0       0
O             -37.23925411     15.48901041     -0.77201800       8       8       1       0       0    1150       0       0       0
O             -31.01707104     17.19948006      2.05092473       8       7       2       1       0    1151       0       0       0
O             -30.45543666     14.75720766      1.10571942       8       4       2       1       1    1152       0       0       0
Si            -35.40009694     23.01228598      1.09405110      14       2       1       0       0    1153       0       0       0
Si            -31.64261488     21.22887979      1.13244276      14       3       1       0       0    1154       0       0       0
O             -33.04222159     20.69330308      1.74901983       8       8       1       0       0    1155       0       0       0
O             -32.38852762     18.98535809     -1.27231417       8       9       1       0       0    1156       0       0       0
O             -35.29701149     22.54457001     -0.41134888       8       7       1       0       0    1157       0       0       0
Si            -33.47204769     19.47946495     -2.18675624      14       1       1       0       1    1158       0       0       0
O             -34.68083353     20.03921162     -1.25935644       8       4       1       0       1    1159       0       0       0
O             -34.04072843     23.63610390      1.73088751       8       5       1       0       1    1160       0       0       0
O             -31.82362872     21.90834201     -0.38929074       8       6       1       0       1    1161       0       0       0
O             -36.63875222     19.00412420      1.12436274       8       9       1       0       0    1162       0       0       0
Si            -37.80402896     19.54385281      0.20386186      14       1       1       0       1    1163       0       0       0
O             -36.02637320     21.85622385      1.99158421       8       6       1       0       1    1164       0       0       0
Si            -31.35025510     23.11510663     -1.21306202      14       2       1       1       0    1165       0       0       0
O             -29.83668660     23.62796910     -0.77549401       8       5       1       1       1    1166       0       0       0
Si            -35.76133430     21.34463460     -1.25912173      14       3       1       0       0    1167       0       0       0
O             -37.13168120     20.79774453     -0.66345349       8       8       1       0       0    1168       0       0       0
O             -31.07100954     22.53332916      2.10251364       8       7       2       1       0    1169       0       0       0
O             -30.59401109     20.13577873      1.18946792       8       4       2       1       1    1170       0       0       0
Si            -35.53390201    -24.90501693      1.16321494      14       2       1       1       0    1171       0       0       0
Si            -31.59704820    -26.56726461      1.05756499      14       3       1       1       0    1172       0       0       0
O             -33.03783907     26.15565331      1.71685042       8       8       1       0       0    1173       0       0       0
O             -32.46165220     24.33721015     -1.33727402       8       9       1       0       0    1174       0       0       0
O             -35.21524488    -25.43644766     -0.38893662       8       7       1       1       0    1175       0       0       0
Si            -33.49015973     24.79470122     -2.19065757      14       1       1       0       1    1176       0       0       0
O             -34.73746036     25.41416224     -1.26032022       8       4       1       0       1    1177       0       0       0
O             -34.02402996    -24.30834240      1.76087168       8       5       1       1       1    1178       0       0       0
O             -31.82347324    -26.09950989     -0.41201421       8       6       1       1       1    1179       0       0       0
O             -36.48570310     24.27090483      1.10830530       8       9       1       0       0    1180       0       0       0
Si            -37.69562288     24.88857185      0.22782816      14       1       1       0       1    1181       0       0       0
O             -35.98484749    -26.11530602      2.07080899       8       6       1       1       1    1182       0       0       0
Si            -31.34887916    -24.91649340     -1.35103978      14       2       1       2       0    1183       0       0       0
O             -29.89492658    -24.35928417     -0.65006310       8       5       1       2       1    1184       0       0       0
Si            -35.70237502     26.60413774     -1.21893311      14       3       1       0       0    1185       0       0       0
O             -37.20845146     26.05737859     -0.73369250       8       8       1       0       0    1186       0       0       0
O             -30.99504241    -25.45321488      2.07514877       8       7       2       2       0    1187       0       0       0
O             -30.42989677     25.55029188      1.16976766       8       4       2       1       1    1188       0       0       0
Si            -35.51494863    -19.44169963      1.20632021      14       2       1       1       0    1189       0       0       0
Si            -31.55275374    -21.32678346      1.18731048      14       3       1       1       0    1190       0       0       0
O             -32.89771318    -21.80442298      1.65763128       8       8       1       1       0    1191       0       0       0
O             -32.45160807    -23.73212219     -1.24729342       8       9       1       1       0    1192       0       0       0
O             -35.18540684    -20.11442825     -0.40753021       8       7       1       1       0    1193       0       0       0
Si            -33.59697520    -23.09958362     -2.23868361      14       1       1       1       1    1194       0       0       0
O             -34.63358493    -22.41201782     -1.35538799       8       4       1       1       1    1195       0       0       0
O             -34.05930503    -18.88595502      1.59499622       8       5       1       1       1    1196       0       0       0
O             -31.84029390    -20.78653555     -0.42354888       8       6       1       1       1    1197       0       0       0
O             -36.62836522    -23.68709395      1.18909340       8       9       1       1       0    1198       0       0       0
Si            -37.69996177    -23.15657529      0.11240215      14       1       1       1       1    1199       0       0       0
O             -35.93149701    -20.72354749      2.09752124       8       6       1       1       1    1200       0       0       0
Si            -31.22465478    -19.47922441     -1.34897284      14       2       1       2       0    1201       0       0       0
O             -29.87518798    -19.04601443     -0.74402725       8       5       1       2       1    1202       0       0       0
Si            -35.69000120    -21.37582424     -1.37067023      14       3       1       1       0    1203       0       0       0
O             -37.15146563    -21.80333599     -0.79121309       8       8       1       1       0    1204       0       0       0
O             -31.12966497    -20.03049297      2.04764412       8       7       2       2       0    1205       0       0       0
O             -30.52577968    -22.49422127      1.14310819       8       4       2       2       1    1206       0       0       0
Si            -35.51356834    -14.18017920      1.20360740      14       2       1       1       0    1207       0       0       0
Si            -31.49097961    -16.03115341      1.10910884      14       3       1       1       0    1208       0       0       0
O             -33.00151007    -16.47988488      1.71383372       8       8       1       1       0    1209       0       0       0
O             -32.26820727    -18.31252280     -1.39060076       8       9       1       1       0    1210       0       0       0
O             -35.25736357    -14.76076683     -0.43253626       8       7       1       1       0    1211       0       0       0
Si            -33.60068714    -17.85952601     -2.25716490      14       1       1       1       1    1212       0       0       0
O             -34.78862389    -17.10642162     -1.26631680       8       4       1       1       1    1213       0       0       0
O             -34.11866765    -13.56900101      1.61739238       8       5       1       1       1    1214       0       0       0
O             -31.88447797    -15.36385613     -0.29082471       8       6       1       1       1    1215       0       0       0
O             -36.53173659    -18.29832140      1.05225398       8       9       1       1       0    1216       0       0       0
Si            -37.74425199    -17.83043219      0.09016391      14       1       1       1       1    1217       0       0       0
O             -35.90287150    -15.35708592      2.07013933       8       6       1       1       1    1218       0       0       0
Si            -31.20303353    -14.29003366     -1.39660151      14       2       1       2       0    1219       0       0       0
O             -29.97098869    -13.58327637     -0.70749855       8       5       1       2       1    1220       0       0       0
Si            -35.80228505    -16.08504565     -1.24898523      14       3       1       1       0    1221       0       0       0
O             -37.17713431    -16.51506194     -0.65445724       8       8       1       1       0    1222       0       0       0
O             -31.05704178    -14.79581699      2.10017213       8       7       2       2       0    1223       0       0       0
O             -30.45557312    -17.23613822      1.10061578       8       4       2       2       1    1224       0       0       0
Si            -31.60536990    -10.69391718      1.13025242      14       3       1       1       0    1226       0       0       0
O             -33.02515743    -11.16501557      1.59943540       8       8       1       1       0    1227       0       0       0
O             -32.28458622    -13.08190705     -1.25794030       8       9       1       1       0    1228       0       0       0
O             -35.35667011     -9.37789429     -0.33742150       8       7       1       1       0    1229       0       0       0
Si            -33.47091870    -12.33907265     -2.30420397      14       1       1       1       1    1230       0       0       0
O             -34.76584503    -11.88759454     -1.27700701       8       4       1       1       1    1231       0       0       0
O             -31.83811986    -10.00388110     -0.31883163       8       6       1       1       1    1233       0       0       0
O             -36.46379781    -13.03303122      1.15768265       8       9       1       1       0    1234       0       0       0
Si            -37.80248154    -12.46048579      0.09071889      14       1       1       1       1    1235       0       0       0
O             -36.06354064    -10.11582943      2.06967050       8       6       1       1       1    1236       0       0       0
Si            -31.21113616     -8.78222598     -1.37606220      14       2       1       2       0    1237       0       0       0
O             -29.88415464     -8.34499205     -0.71922886       8       5       1       2       1    1238       0       0       0
Si            -35.81941222    -10.57291289     -1.33369750      14       3       1       1       0    1239       0       0       0
O             -37.26929738    -11.25611205     -0.72042092       8       8       1       1       0    1240       0       0       0
O             -31.12371198     -9.37851902      1.97724997       8       7       2       2       0    1241       0       0       0
O             -30.58256510    -11.78072262      1.17517172       8       4       2       2       1    1242       0       0       0
Si            -27.18126345      7.15829607      1.09725614      14       2       1       0       0    1279       0       0       0
O             -25.70830688      7.64650646      1.62198042       8       5       1       0       1    1286       0       0       0
O             -23.47919428      5.93042077     -0.36827000       8       6       1       0       1    1287       0       0       0
Si            -22.80872093      7.06413623     -1.32212005      14       2       1       1       0    1291       0       0       0
O             -21.43430221      7.59275458     -0.73491624       8       5       1       1       1    1292       0       0       0
O             -22.67615484      6.58560764      2.03997988       8       7       2       1       0    1295       0       0       0
Si            -27.10038638     12.39695205      1.08318055      14       2       1       0       0    1297       0       0       0
Si            -23.18759259     10.60330087      1.09333355      14       3       1       0       0    1298       0       0       0
O             -24.53631891     10.16352084      1.77001087       8       8       1       0       0    1299       0       0       0
O             -24.06291475      8.30906427     -1.29756457       8       9       1       0       0    1300       0       0       0
O             -26.86993722     11.79201265     -0.28688911       8       7       1       0       0    1301       0       0       0
Si            -25.17816064      8.83076975     -2.34045368      14       1       1       0       1    1302       0       0       0
O             -26.29889157      9.55752366     -1.38932000       8       4       1       0       1    1303       0       0       0
O             -25.61718251     12.96571659      1.71661498       8       5       1       0       1    1304       0       0       0
O             -23.44681138     11.16312191     -0.34470248       8       6       1       0       1    1305       0       0       0
O             -28.19625996      8.23293464      1.20082534       8       9       1       0       0    1306       0       0       0
Si            -29.28512279      8.93411830      0.17252960      14       1       1       0       1    1307       0       0       0
O             -27.58825317     11.28353194      2.04435999       8       6       1       0       1    1308       0       0       0
Si            -22.92646208     12.42081106     -1.35792228      14       2       1       1       0    1309       0       0       0
O             -21.47925912     13.07671466     -0.73108787       8       5       1       1       1    1310       0       0       0
Si            -27.41667541     10.66055055     -1.24347204      14       3       1       0       0    1311       0       0       0
O             -28.85844185     10.05503655     -0.66423333       8       8       1       0       0    1312       0       0       0
O             -22.69531711     11.77197998      1.99506438       8       7       2       1       0    1313       0       0       0
O             -22.16097618      9.54238664      1.00339042       8       4       2       1       1    1314       0       0       0
Si            -27.15353798     17.71686532      1.11898551      14       2       1       0       0    1315       0       0       0
Si            -23.12042788     16.00882416      1.05998293      14       3       1       0       0    1316       0       0       0
O             -24.55584372     15.49775688      1.69313704       8       8       1       0       0    1317       0       0       0
O             -23.96319750     13.56254986     -1.29436419       8       9       1       0       0    1318       0       0       0
O             -26.82800715     17.27005319     -0.37072524       8       7       1       0       0    1319       0       0       0
Si            -25.14296310     14.23994042     -2.21134577      14       1       1       0       1    1320       0       0       0
O             -26.27964632     14.71987475     -1.37848393       8       4       1       0       1    1321       0       0       0
O             -25.62283105     18.40086319      1.73278572       8       5       1       0       1    1322       0       0       0
O             -23.35285263     16.66065825     -0.44258334       8       6       1       0       1    1323       0       0       0
O             -28.24136619     13.63722237      1.03614481       8       9       1       0       0    1324       0       0       0
Si            -29.39062239     14.24710838      0.26952413      14       1       1       0       1    1325       0       0       0
O             -27.70294680     16.47557267      2.14761465       8       6       1       0       1    1326       0       0       0
Si            -22.80950432     17.85956182     -1.20885127      14       2       1       1       0    1327       0       0       0
O             -21.56977331     18.31016919     -0.66743230       8       5       1       1       1    1328       0       0       0
Si            -27.36535606     15.89035277     -1.30851364      14       3       1       0       0    1329       0       0       0
O             -28.86773835     15.46255376     -0.66136034       8       8       1       0       0    1330       0       0       0
O             -22.74329157     17.27166610      2.01502190       8       7       2       1       0    1331       0       0       0
O             -22.12729518     14.84164613      1.01167530       8       4       2       1       1    1332       0       0       0
Si            -27.17614104     22.99081784      1.13484622      14       2       1       0       0    1333       0       0       0
Si            -23.22780816     21.22676812      1.08605746      14       3       1       0       0    1334       0       0       0
O             -24.58617406     20.74621881      1.67748479       8       8       1       0       0    1335       0       0       0
O             -24.01512084     18.92033328     -1.32409633       8       9       1       0       0    1336       0       0       0
O             -26.86615860     22.58190629     -0.35602479       8       7       1       0       0    1337       0       0       0
Si            -25.14400364     19.61265522     -2.23277131      14       1       1       0       1    1338       0       0       0
O             -26.34163084     20.22015621     -1.23372959       8       4       1       0       1    1339       0       0       0
O             -25.66968900     23.71816816      1.61830782       8       5       1       0       1    1340       0       0       0
O             -23.38967948     21.97420384     -0.27743524       8       6       1       0       1    1341       0       0       0
O             -28.13739393     18.91712028      1.14853952       8       9       1       0       0    1342       0       0       0
Si            -29.33095194     19.46015674      0.08924582      14       1       1       0       1    1343       0       0       0
O             -27.64565287     21.97189992      2.13461671       8       6       1       0       1    1344       0       0       0
Si            -22.88261642     23.05096520     -1.31726827      14       2       1       1       0    1345       0       0       0
O             -21.42990288     23.65886823     -0.78197354       8       5       1       1       1    1346       0       0       0
Si            -27.37809802     21.34330207     -1.39671847      14       3       1       0       0    1347       0       0       0
O             -28.68974645     20.75038476     -0.79996013       8       8       1       0       0    1348       0       0       0
O             -22.69724444     22.48972421      2.10512848       8       7       2       1       0    1349       0       0       0
O             -22.07993766     20.09988253      1.17200696       8       4       2       1       1    1350       0       0       0
Si            -27.18961099    -24.94190914      1.07160330      14       2       1       1       0    1351       0       0       0
Si            -23.15263752    -26.63584764      1.09758136      14       3       1       1       0    1352       0       0       0
O             -24.59011528     26.00445975      1.60115399       8       8       1       0       0    1353       0       0       0
O             -24.08189837     24.33267584     -1.25771802       8       9       1       0       0    1354       0       0       0
O             -26.92189863    -25.41620232     -0.30924653       8       7       1       1       0    1355       0       0       0
Si            -25.14439086     24.78083074     -2.20669118      14       1       1       0       1    1356       0       0       0
O             -26.26616036     25.54902937     -1.28264714       8       4       1       0       1    1357       0       0       0
O             -25.79832450    -24.22562288      1.74947832       8       5       1       1       1    1358       0       0       0
O             -23.41556528    -26.11232376     -0.27604508       8       6       1       1       1    1359       0       0       0
O             -28.11465565     24.18786283      1.12684601       8       9       1       0       0    1360       0       0       0
Si            -29.31333358     24.86211972      0.13115934      14       1       1       0       1    1361       0       0       0
O             -27.53097945    -26.08734045      2.04815977       8       6       1       1       1    1362       0       0       0
Si            -22.99029260    -24.85714305     -1.39100563      14       2       1       2       0    1363       0       0       0
O             -21.56168884    -24.24978147     -0.66582436       8       5       1       2       1    1364       0       0       0
Si            -27.31362823    -26.58317606     -1.25930753      14       3       1       1       0    1365       0       0       0
O             -28.73266346     25.96059780     -0.66748700       8       8       1       0       0    1366       0       0       0
O             -22.66549769    -25.40831733      2.02795664       8       7       2       2       0    1367       0       0       0
O             -22.06126643     25.41023923      1.19073773       8       4       2       1       1    1368       0       0       0
Si            -27.01792307    -19.44064083      1.19619964      14       2       1       1       0    1369       0       0       0
Si            -23.11103877    -21.24484767      1.20059596      14       3       1       1       0    1370       0       0       0
O             -24.49406750    -21.88674686      1.57109759       8       8       1       1       0    1371       0       0       0
O             -24.05503181    -23.65544863     -1.34297000       8       9       1       1       0    1372       0       0       0
O             -26.96323024    -20.18423730     -0.42818622       8       7       1       1       0    1373       0       0       0
Si            -25.23666398    -23.14502036     -2.16022524      14       1       1       1       1    1374       0       0       0
O             -26.41680418    -22.57473696     -1.31567660       8       4       1       1       1    1375       0       0       0
O             -25.61457900    -18.97762755      1.66893661       8       5       1       1       1    1376       0       0       0
O             -23.50642259    -20.80697384     -0.35526704       8       6       1       1       1    1377       0       0       0
O             -28.10650720    -23.67274719      1.11583082       8       9       1       1       0    1378       0       0       0
Si            -29.30621559    -23.06257509      0.18230517      14       1       1       1       1    1379       0       0       0
O             -27.61975742    -20.80492245      2.15096138       8       6       1       1       1    1380       0       0       0
Si            -22.93589991    -19.56908208     -1.21992017      14       2       1       2       0    1381       0       0       0
O             -21.54632057    -19.03912069     -0.69237731       8       5       1       2       1    1382       0       0       0
Si            -27.44670795    -21.30894378     -1.33225161      14       3       1       1       0    1383       0       0       0
O             -28.70346089    -21.86876915     -0.69223368       8       8       1       1       0    1384       0       0       0
O             -22.64031725    -20.09955800      1.99657631       8       7       2       2       0    1385       0       0       0
O             -22.22267490    -22.44276569      1.17381768       8       4       2       2       1    1386       0       0       0
Si            -27.10482752    -14.26228410      1.17799336      14       2       1       1       0    1387       0       0       0
Si            -23.24932259    -15.96193564      1.17858684      14       3       1       1       0    1388       0       0       0
O             -24.68797836    -16.60973739      1.67680906       8       8       1       1       0    1389       0       0       0
O             -23.90602530    -18.42872422     -1.32838098       8       9       1       1       0    1390       0       0       0
O             -26.82313757    -14.83307613     -0.40732498       8       7       1       1       0    1391       0       0       0
Si            -25.18541029    -17.84747538     -2.18505118      14       1       1       1       1    1392       0       0       0
O             -26.26044122    -17.08530075     -1.22523223       8       4       1       1       1    1393       0       0       0
O             -25.74401535    -13.57923317      1.67977572       8       5       1       1       1    1394       0       0       0
O             -23.46379948    -15.48054875     -0.27220978       8       6       1       1       1    1395       0       0       0
O             -28.08923275    -18.28922390      1.05824252       8       9       1       1       0    1396       0       0       0
Si            -29.24786696    -17.79131561      0.17695472      14       1       1       1       1    1397       0       0       0
O             -27.61255779    -15.31855390      2.03317173       8       6       1       1       1    1398       0       0       0
Si            -22.81485185    -14.21029786     -1.19843684      14       2       1       2       0    1399       0       0       0
O             -21.43945083    -13.55012771     -0.79905783       8       5       1       2       1    1400       0       0       0
Si            -27.39386922    -16.04146605     -1.26104598      14       3       1       1       0    1401       0       0       0
O             -28.81111472    -16.59686463     -0.69628075       8       8       1       1       0    1402       0       0       0
O             -22.64337485    -14.82076104      1.97398454       8       7       2       2       0    1403       0       0       0
O             -22.13673137    -17.25758906      1.19822668       8       4       2       2       1    1404       0       0       0
Si            -27.07608919     -8.91467045      1.04165796      14       2       1       1       0    1405       0       0       0
Si            -23.10750161    -10.62510605      1.14620838      14       3       1       1       0    1406       0       0       0
O             -24.53898835    -11.24847039      1.60093948       8       8       1       1       0    1407       0       0       0
O             -24.01610640    -12.94210745     -1.31425875       8       9       1       1       0    1408       0       0       0
O             -26.95953900     -9.51167037     -0.43398226       8       7       1       1       0    1409       0       0       0
Si            -25.21437619    -12.46225881     -2.19054673      14       1       1       1       1    1410       0       0       0
O             -26.31438256    -11.88442073     -1.26383185       8       4       1       1       1    1411       0       0       0
O             -25.72552253     -8.28579466      1.69069908       8       5       1       1       1    1412       0       0       0
O             -23.36524404    -10.13467304     -0.44974261       8       6       1       1       1    1413       0       0       0
O             -28.08821409    -13.08330491      1.01017728       8       9       1       1       0    1414       0       0       0
Si            -29.40767724    -12.34532629      0.10959780      14       1       1       1       1    1415       0       0       0
O             -27.54525550    -10.09698738      2.16194535       8       6       1       1       1    1416       0       0       0
Si            -22.81375990     -8.86345371     -1.34587704      14       2       1       2       0    1417       0       0       0
O             -21.42515601     -8.38511184     -0.69891940       8       5       1       2       1    1418       0       0       0
Si            -27.48990321    -10.72887641     -1.34055003      14       3       1       1       0    1419       0       0       0
O             -28.87414778    -11.32250738     -0.79801744       8       8       1       1       0    1420       0       0       0
O             -22.78428628     -9.40318131      2.01911196       8       7       2       2       0    1421       0       0       0
O             -22.16532382    -11.92786223      1.06118304       8       4       2       2       1    1422       0       0       0
O             -23.89222627     -7.61290936     -1.29804957       8       9       1       1       0    1426       0       0       0
Si            -25.19570962     -7.01196613     -2.33801931      14       1       1       1       1    1428       0       0       0
O             -28.12324223     -7.78494141      1.10860604       8       9       1       1       0    1432       0       0       0
O             -22.03359596     -6.58077418      1.02649501       8       4       2       2       1    1440       0       0       0
Si            -18.68742441      7.16759871      1.09719514      14       2       1       0       0    1459       0       0       0
Si            -14.91491082      5.41298850      1.21782954      14       3       1       0       0    1460       0       0       0
O             -16.20674614      4.82712294      1.67688914       8       8       1       0       0    1461       0       0       0
O             -18.39672757      6.47193438     -0.38357351       8       7       1       0       0    1463       0       0       0
O             -17.38038731      7.69454931      1.73846127       8       5       1       0       1    1466       0       0       0
O             -15.13208268      5.97643642     -0.38068950       8       6       1       0       1    1467       0       0       0
O             -19.32438666      5.88296301      1.99869699       8       6       1       0       1    1470       0       0       0
Si            -14.56578749      7.11225278     -1.23704684      14       2       1       1       0    1471       0       0       0
O             -13.10526945      7.59077670     -0.84606867       8       5       1       1       1    1472       0       0       0
Si            -19.05475155      5.40386117     -1.27967245      14       3       1       0       0    1473       0       0       0
O             -14.28355069      6.42482387      2.16703156       8       7       2       1       0    1475       0       0       0
O             -13.76796747      4.20071537      1.02864658       8       4       2       1       1    1476       0       0       0
Si            -18.67284139     12.36029355      1.15790415      14       2       1       0       0    1477       0       0       0
Si            -14.82257229     10.56317575      1.04238250      14       3       1       0       0    1478       0       0       0
O             -16.18334398     10.06478400      1.60342301       8       8       1       0       0    1479       0       0       0
O             -15.62849704      8.38827359     -1.22075179       8       9       1       0       0    1480       0       0       0
O             -18.53614409     11.76721444     -0.29131826       8       7       1       0       0    1481       0       0       0
Si            -16.84093012      8.79399304     -2.29785876      14       1       1       0       1    1482       0       0       0
O             -17.84855109      9.54159154     -1.24915796       8       4       1       0       1    1483       0       0       0
O             -17.26419008     13.08222817      1.57818937       8       5       1       0       1    1484       0       0       0
O             -14.95398449     11.30577298     -0.34552986       8       6       1       0       1    1485       0       0       0
O             -19.85399269      8.28970586      1.02006739       8       9       1       0       0    1486       0       0       0
Si            -20.87283327      8.86637176      0.26758480      14       1       1       0       1    1487       0       0       0
O             -19.14633244     11.22199275      2.08415940       8       6       1       0       1    1488       0       0       0
Si            -14.56474421     12.44966568     -1.28564882      14       2       1       1       0    1489       0       0       0
O             -13.17833751     13.06374789     -0.78552655       8       5       1       1       1    1490       0       0       0
Si            -19.00340365     10.71390641     -1.22681617      14       3       1       0       0    1491       0       0       0
O             -20.46482087     10.10327676     -0.74309026       8       8       1       0       0    1492       0       0       0
O             -14.31323977     11.86098482      1.97095886       8       7       2       1       0    1493       0       0       0
O             -13.80701076      9.44808295      1.03875141       8       4       2       1       1    1494       0       0       0
Si            -18.62387082     17.72572969      1.11008707      14       2       1       0       0    1495       0       0       0
Si            -14.82303214     15.94282157      1.03938869      14       3       1       0       0    1496       0       0       0
O             -16.11720228     15.36007627      1.62720083       8       8       1       0       0    1497       0       0       0
O             -15.57018757     13.66121683     -1.36833899       8       9       1       0       0    1498       0       0       0
O             -18.49668972     17.22950383     -0.45252955       8       7       1       0       0    1499       0       0       0
Si            -16.67403249     14.22892130     -2.20820976      14       1       1       0       1    1500       0       0       0
O             -18.00732921     14.73870121     -1.26364012       8       4       1       0       1    1501       0       0       0
O             -17.29869127     18.28535677      1.73513918       8       5       1       0       1    1502       0       0       0
O             -15.00921253     16.62344032     -0.28232535       8       6       1       0       1    1503       0       0       0
O             -19.80430083     13.72040502      1.03207741       8       9       1       0       0    1504       0       0       0
Si            -21.04759298     14.11257934      0.25141099      14       1       1       0       1    1505       0       0       0
O             -19.23493856     16.52008272      2.09918179       8       6       1       0       1    1506       0       0       0
Si            -14.52842633     17.82020925     -1.25757034      14       2       1       1       0    1507       0       0       0
O             -13.07191500     18.39031082     -0.73630683       8       5       1       1       1    1508       0       0       0
Si            -19.06235776     15.93476800     -1.20886068      14       3       1       0       0    1509       0       0       0
O             -20.42359264     15.41181940     -0.79875175       8       8       1       0       0    1510       0       0       0
O             -14.35362555     17.19413275      2.12491607       8       7       2       1       0    1511       0       0       0
O             -13.70678560     14.83840528      1.10082413       8       4       2       1       1    1512       0       0       0
Si            -18.61536931     23.08300542      1.20157345      14       2       1       0       0    1513       0       0       0
Si            -14.88631117     21.29506571      1.07863315      14       3       1       0       0    1514       0       0       0
O             -16.14569648     20.69096682      1.57411870       8       8       1       0       0    1515       0       0       0
O             -15.65297820     18.90085204     -1.40543985       8       9       1       0       0    1516       0       0       0
O             -18.58195771     22.42404027     -0.27272426       8       7       1       0       0    1517       0       0       0
Si            -16.83980823     19.56709802     -2.21855677      14       1       1       0       1    1518       0       0       0
O             -17.84893655     20.15656652     -1.27236430       8       4       1       0       1    1519       0       0       0
O             -17.24744883     23.71247046      1.67194166       8       5       1       0       1    1520       0       0       0
O             -15.08373810     21.79945624     -0.34113635       8       6       1       0       1    1521       0       0       0
O             -19.70151598     18.94084087      1.07455133       8       9       1       0       0    1522       0       0       0
Si            -20.88947832     19.46243202      0.22993080      14       1       1       0       1    1523       0       0       0
O             -19.24655213     21.90414897      2.11891226       8       6       1       0       1    1524       0       0       0
Si            -14.52784321     23.18621861     -1.25679585      14       2       1       1       0    1525       0       0       0
O             -13.10729376     23.67485717     -0.65674268       8       5       1       1       1    1526       0       0       0
Si            -19.08208359     21.28902601     -1.24302940      14       3       1       0       0    1527       0       0       0
O             -20.49892230     20.79213293     -0.65353664       8       8       1       0       0    1528       0       0       0
O             -14.34855759     22.59917989      2.08559450       8       7       2       1       0    1529       0       0       0
O             -13.84098556     20.08007711      1.12562425       8       4       2       1       1    1530       0       0       0
Si            -18.67844623    -24.80584534      1.15120405      14       2       1       1       0    1531       0       0       0
Si            -14.81150331    -26.54339173      1.08866769      14       3       1       1       0    1532       0       0       0
O             -16.24565620     25.98025378      1.57767513       8       8       1       0       0    1533       0       0       0
O             -15.54313544     24.28792374     -1.26078277       8       9       1       0       0    1534       0       0       0
O             -18.56439579    -25.38847112     -0.31061370       8       7       1       1       0    1535       0       0       0
Si            -16.81859541     24.84170419     -2.33098048      14       1       1       0       1    1536       0       0       0
O             -17.97750442     25.44033134     -1.38870168       8       4       1       0       1    1537       0       0       0
O             -17.38597086    -24.23153375      1.70833737       8       5       1       1       1    1538       0       0       0
O             -15.12400246    -26.13821966     -0.34564380       8       6       1       1       1    1539       0       0       0
O             -19.71707356     24.25381953      1.04001533       8       9       1       0       0    1540       0       0       0
Si            -21.04182669     24.87918660      0.08459434      14       1       1       0       1    1541       0       0       0
O             -19.19674102    -26.11017181      2.08982167       8       6       1       1       1    1542       0       0       0
Si            -14.59737294    -24.79845460     -1.38607475      14       2       1       2       0    1543       0       0       0
O             -13.23112456    -24.22854407     -0.80640610       8       5       1       2       1    1544       0       0       0
Si            -19.00223631    -26.54951150     -1.29626454      14       3       1       1       0    1545       0       0       0
O             -20.32045690     26.06248655     -0.83031494       8       8       1       0       0    1546       0       0       0
O             -14.22956866    -25.52987150      2.04932043       8       7       2       2       0    1547       0       0       0
O             -13.77638401     25.53849984      1.00579560       8       4       2       1       1    1548       0       0       0
Si            -18.80791424    -19.47107415      1.04728616      14       2       1       1       0    1549       0       0       0
Si            -14.76162520    -21.30521238      1.21140591      14       3       1       1       0    1550       0       0       0
O             -16.14878884    -21.83156664      1.68284154       8       8       1       1       0    1551       0       0       0
O             -15.52547828    -23.59866358     -1.28112577       8       9       1       1       0    1552       0       0       0
O             -18.45962356    -20.06743620     -0.37709723       8       7       1       1       0    1553       0       0       0
Si            -16.70601577    -23.07601989     -2.19416254      14       1       1       1       1    1554       0       0       0
O             -17.85646751    -22.48247985     -1.23914024       8       4       1       1       1    1555       0       0       0
O             -17.33439126    -18.88431213      1.75310776       8       5       1       1       1    1556       0       0       0
O             -15.03002514    -20.70492160     -0.35335404       8       6       1       1       1    1557       0       0       0
O             -19.78984257    -23.73370578      1.15429448       8       9       1       1       0    1558       0       0       0
Si            -20.92106764    -23.05197300      0.21519607      14       1       1       1       1    1559       0       0       0
O             -19.19690171    -20.72296218      2.15779858       8       6       1       1       1    1560       0       0       0
Si            -14.51739039    -19.58527194     -1.34848618      14       2       1       2       0    1561       0       0       0
O             -13.06757073    -19.05303295     -0.73685322       8       5       1       2       1    1562       0       0       0
Si            -19.01253408    -21.36699646     -1.38804242      14       3       1       1       0    1563       0       0       0
O             -20.33690700    -21.86411376     -0.66448557       8       8       1       1       0    1564       0       0       0
O             -14.34529097    -20.04698839      2.15659053       8       7       2       2       0    1565       0       0       0
O             -13.80837532    -22.43939908      1.10195371       8       4       2       2       1    1566       0       0       0
Si            -18.61937189    -14.28777587      1.14061690      14       2       1       1       0    1567       0       0       0
Si            -14.77340917    -16.07134702      1.17778178      14       3       1       1       0    1568       0       0       0
O             -16.25802749    -16.51798463      1.62093556       8       8       1       1       0    1569       0       0       0
O             -15.58292884    -18.38951109     -1.34261995       8       9       1       1       0    1570       0       0       0
O             -18.40250043    -14.83976782     -0.39493015       8       7       1       1       0    1571       0       0       0
Si            -16.85142076    -17.86134662     -2.23902603      14       1       1       1       1    1572       0       0       0
O             -18.01926033    -17.25730779     -1.22862668       8       4       1       1       1    1573       0       0       0
O             -17.30133028    -13.69754685      1.64814946       8       5       1       1       1    1574       0       0       0
O             -15.09669150    -15.46896889     -0.44364932       8       6       1       1       1    1575       0       0       0
O             -19.73803200    -18.35741504      1.13827616       8       9       1       1       0    1576       0       0       0
Si            -20.94583391    -17.69891952      0.17240980      14       1       1       1       1    1577       0       0       0
O             -19.20024547    -15.33838130      2.02684161       8       6       1       1       1    1578       0       0       0
Si            -14.59989034    -14.12739942     -1.31361935      14       2       1       2       0    1579       0       0       0
O             -13.21640775    -13.63594785     -0.67661910       8       5       1       2       1    1580       0       0       0
Si            -18.98866575    -15.94382991     -1.23259353      14       3       1       1       0    1581       0       0       0
O             -20.41993246    -16.57263627     -0.81495189       8       8       1       1       0    1582       0       0       0
O             -14.38420949    -14.88415992      2.10456587       8       7       2       2       0    1583       0       0       0
O             -13.76986150    -17.26092931      1.20130511       8       4       2       2       1    1584       0       0       0
Si            -18.67891712     -8.93168487      1.18026527      14       2       1       1       0    1585       0       0       0
Si            -14.88387960    -10.56072366      1.16840268      14       3       1       1       0    1586       0       0       0
O             -16.13084001    -11.20339565      1.63307332       8       8       1       1       0    1587       0       0       0
O             -15.50920663    -13.06435224     -1.38652389       8       9       1       1       0    1588       0       0       0
O             -18.40123837     -9.42435178     -0.27978684       8       7       1       1       0    1589       0       0       0
Si            -16.72155503    -12.40095345     -2.22471758      14       1       1       1       1    1590       0       0       0
O             -17.95470948    -11.93983356     -1.34039997       8       4       1       1       1    1591       0       0       0
O             -17.27640741     -8.32384474      1.72218747       8       5       1       1       1    1592       0       0       0
O             -15.08399282    -10.01757133     -0.30892570       8       6       1       1       1    1593       0       0       0
O             -19.82818535    -13.10224918      1.13123906       8       9       1       1       0    1594       0       0       0
Si            -20.93828515    -12.35634585      0.11485634      14       1       1       1       1    1595       0       0       0
O             -19.21451366     -9.98614326      2.14464478       8       6       1       1       1    1596       0       0       0
Si            -14.44572092     -8.81025624     -1.25299588      14       2       1       2       0    1597       0       0       0
O             -13.18234421     -8.37266673     -0.78859505       8       5       1       2       1    1598       0       0       0
Si            -18.96045229    -10.64401331     -1.38407738      14       3       1       1       0    1599       0       0       0
O             -20.31367606    -11.20026710     -0.66656188       8       8       1       1       0    1600       0       0       0
O             -14.40145341     -9.38800841      2.12630922       8       7       2       2       0    1601       0       0       0
O             -13.74885646    -11.92933638      1.10148917       8       4       2       2       1    1602       0       0       0
Si            -14.75313243     -5.37292074      1.14650590      14       3       1       1       0    1604       0       0       0
O             -16.27323795     -5.89479761      1.65077236       8       8       1       1       0    1605       0       0       0
O             -15.52794406     -7.64547061     -1.41063444       8       9       1       1       0    1606       0       0       0
Si            -16.68911675     -7.14576570     -2.29753111      14       1       1       1       1    1608       0       0       0
O             -18.03042482     -6.49849228     -1.22395753       8       4       1       1       1    1609       0       0       0
O             -14.94530811     -4.64412957     -0.37191195       8       6       1       1       1    1611       0       0       0
O             -19.76916841     -7.62652863      1.14909788       8       9       1       1       0    1612       0       0       0
Si            -20.97593795     -7.12190309      0.13393467      14       1       1       1       1    1613       0       0       0
Si            -19.05136718     -5.39420186     -1.30999246      14       3       1       1       0    1617       0       0       0
O             -20.43718288     -5.95428809     -0.83263599       8       8       1       1       0    1618       0       0       0
O             -14.28017161     -4.21833309      2.13819680       8       7       2       2       0    1619       0       0       0
O             -13.71390803     -6.55637138      1.18503828       8       4       2       2       1    1620       0       0       0
O              -7.19975535     -2.31737196     -1.36076225       8       9       1       0       0    1624       0       0       0
O              -8.86235412      2.34852204      1.72894417       8       5       1       0       1    1628       0       0       0
Si             -6.20412001      1.71921660     -1.34978090      14       2       1       1       0    1633       0       0       0
O              -4.74778219      2.33529035     -0.80373041       8       5       1       1       1    1634       0       0       0
Si            -10.33361272      7.14851095      1.18005238      14       2       1       0       0    1639       0       0       0
Si             -6.34867597      5.41055975      1.07605631      14       3       1       0       0    1640       0       0       0
O              -7.81551332      4.66821228      1.77074097       8       8       1       0       0    1641       0       0       0
O              -7.29842977      2.90407605     -1.31967796       8       9       1       0       0    1642       0       0       0
O             -10.16927450      6.44658952     -0.43586597       8       7       1       0       0    1643       0       0       0
Si             -8.31293817      3.58907087     -2.32438368      14       1       1       0       1    1644       0       0       0
O              -9.51099593      4.07988238     -1.30611810       8       4       1       0       1    1645       0       0       0
O              -8.86238357      7.71034408      1.67309111       8       5       1       0       1    1646       0       0       0
O              -6.69197523      5.85359084     -0.26568417       8       6       1       0       1    1647       0       0       0
O             -11.40168274      3.03116969      1.18193765       8       9       1       0       0    1648       0       0       0
Si            -12.54559903      3.51557076      0.16648440      14       1       1       0       1    1649       0       0       0
O             -10.83083914      5.94717844      2.00883464       8       6       1       0       1    1650       0       0       0
Si             -6.19137476      7.03444538     -1.28391355      14       2       1       1       0    1651       0       0       0
O              -4.76411432      7.65061885     -0.78737449       8       5       1       1       1    1652       0       0       0
Si            -10.68028190      5.28797229     -1.26667995      14       3       1       0       0    1653       0       0       0
O             -11.94896995      4.77080891     -0.74610036       8       8       1       0       0    1654       0       0       0
O              -5.82053873      6.55614369      2.13617523       8       7       2       1       0    1655       0       0       0
O              -5.44697436      4.20613015      1.02812052       8       4       2       1       1    1656       0       0       0
Si            -10.28677463     12.36989206      1.07359477      14       2       1       0       0    1657       0       0       0
Si             -6.47198820     10.70913396      1.14289970      14       3       1       0       0    1658       0       0       0
O              -7.80111299     10.11470480      1.69611228       8       8       1       0       0    1659       0       0       0
O              -7.26376421      8.30765104     -1.33589923       8       9       1       0       0    1660       0       0       0
O             -10.07956559     11.87852419     -0.27253284       8       7       1       0       0    1661       0       0       0
Si             -8.40955633      8.79887236     -2.18473181      14       1       1       0       1    1662       0       0       0
O              -9.59321722      9.39497959     -1.31024807       8       4       1       0       1    1663       0       0       0
O              -8.85901884     13.00106810      1.65207519       8       5       1       0       1    1664       0       0       0
O              -6.56930466     11.28405323     -0.37321565       8       6       1       0       1    1665       0       0       0
O             -11.38457302      8.22970649      1.02457021       8       9       1       0       0    1666       0       0       0
Si            -12.63278510      8.85866917      0.26085266      14       1       1       0       1    1667       0       0       0
O             -10.76181354     11.24136828      2.01444631       8       6       1       0       1    1668       0       0       0
Si             -6.08757664     12.52023182     -1.36623371      14       2       1       1       0    1669       0       0       0
O              -4.69926063     12.92899301     -0.74711800       8       5       1       1       1    1670       0       0       0
Si            -10.68755992     10.72284951     -1.30578817      14       3       1       0       0    1671       0       0       0
O             -11.95206164     10.16914889     -0.83431577       8       8       1       0       0    1672       0       0       0
O              -5.96985519     11.80194611      2.08265379       8       7       2       1       0    1673       0       0       0
O              -5.33396437      9.54361588      1.14147949       8       4       2       1       1    1674       0       0       0
Si            -10.28741232     17.71980671      1.12150316      14       2       1       0       0    1675       0       0       0
Si             -6.38953462     16.04537351      1.03797214      14       3       1       0       0    1676       0       0       0
O              -7.84814712     15.39281240      1.60700212       8       8       1       0       0    1677       0       0       0
O              -7.16342506     13.64418799     -1.39523096       8       9       1       0       0    1678       0       0       0
O             -10.02917993     17.12434907     -0.33595876       8       7       1       0       0    1679       0       0       0
Si             -8.36630177     14.24692034     -2.32976926      14       1       1       0       1    1680       0       0       0
O              -9.63409946     14.72386540     -1.34160707       8       4       1       0       1    1681       0       0       0
O              -9.01557628     18.35841506      1.60050060       8       5       1       0       1    1682       0       0       0
O              -6.73675193     16.56916888     -0.37035403       8       6       1       0       1    1683       0       0       0
O             -11.47273753     13.69210746      1.12682743       8       9       1       0       0    1684       0       0       0
Si            -12.65039248     14.20569105      0.14531160      14       1       1       0       1    1685       0       0       0
O             -10.79403714     16.62724929      2.14673608       8       6       1       0       1    1686       0       0       0
Si             -6.22978504     17.74180554     -1.31330101      14       2       1       1       0    1687       0       0       0
O              -4.80116691     18.30805895     -0.66239012       8       5       1       1       1    1688       0       0       0
Si            -10.55409091     16.04301744     -1.23270571      14       3       1       0       0    1689       0       0       0
O             -11.92338917     15.38453378     -0.82579267       8       8       1       0       0    1690       0       0       0
O              -5.89678982     17.17636881      1.98833234       8       7       2       1       0    1691       0       0       0
O              -5.26504086     14.80959293      1.02833446       8       4       2       1       1    1692       0       0       0
Si            -10.42086992     23.15420383      1.15108708      14       2       1       0       0    1693       0       0       0
Si             -6.34829675     21.33036001      1.16031184      14       3       1       0       0    1694       0       0       0
O              -7.86394085     20.71059351      1.72907134       8       8       1       0       0    1695       0       0       0
O              -7.15559557     18.88530763     -1.30995288       8       9       1       0       0    1696       0       0       0
O             -10.04764129     22.54596186     -0.31757970       8       7       1       0       0    1697       0       0       0
Si             -8.31178605     19.61456791     -2.20943649      14       1       1       0       1    1698       0       0       0
O              -9.56047124     20.04432821     -1.35077208       8       4       1       0       1    1699       0       0       0
O              -8.85524260     23.66824456      1.70117259       8       5       1       0       1    1700       0       0       0
O              -6.70622583     21.98573420     -0.37183207       8       6       1       0       1    1701       0       0       0
O             -11.32039392     18.96589210      1.12980534       8       9       1       0       0    1702       0       0       0
Si            -12.55185473     19.54697205      0.26032771      14       1       1       0       1    1703       0       0       0
O             -10.75260833     21.81152221      2.05353335       8       6       1       0       1    1704       0       0       0
Si             -6.23113857     23.09733072     -1.29556357      14       2       1       1       0    1705       0       0       0
O              -4.70444526     23.70591003     -0.78270820       8       5       1       1       1    1706       0       0       0
Si            -10.61679334     21.24341438     -1.28105297      14       3       1       0       0    1707       0       0       0
O             -12.01960814     20.82063788     -0.67643710       8       8       1       0       0    1708       0       0       0
O              -5.96603060     22.44476976      2.02426908       8       7       2       1       0    1709       0       0       0
O              -5.39492238     20.15743763      1.16865050       8       4       2       1       1    1710       0       0       0
Si            -10.37966779    -24.91157712      1.12728885      14       2       1       1       0    1711       0       0       0
Si             -6.49415772     26.56430619      1.06860866      14       3       1       0       0    1712       0       0       0
O              -7.76626789     26.01006909      1.65495877       8       8       1       0       0    1713       0       0       0
O              -7.21491114     24.24854463     -1.38501911       8       9       1       0       0    1714       0       0       0
O             -10.19312913    -25.45823006     -0.27017237       8       7       1       1       0    1715       0       0       0
Si             -8.32838893     24.81331664     -2.33519739      14       1       1       0       1    1716       0       0       0
O              -9.54545159     25.49667721     -1.29761736       8       4       1       0       1    1717       0       0       0
O              -8.89349392    -24.22027605      1.59185966       8       5       1       1       1    1718       0       0       0
O              -6.65472569    -25.99969381     -0.38995566       8       6       1       1       1    1719       0       0       0
O             -11.42729384     24.32486398      1.03036077       8       9       1       0       0    1720       0       0       0
Si            -12.47842930     24.82037400      0.17330576      14       1       1       0       1    1721       0       0       0
O             -10.94404879    -26.11579738      2.09475877       8       6       1       1       1    1722       0       0       0
Si             -6.05111457    -24.94461968     -1.23094598      14       2       1       2       0    1723       0       0       0
O              -4.66330760    -24.30237805     -0.75911858       8       5       1       2       1    1724       0       0       0
Si            -10.65172045     26.56007335     -1.24290725      14       3       1       0       0    1725       0       0       0
O             -12.07519997     26.08896687     -0.73005360       8       8       1       0       0    1726       0       0       0
O              -5.89403156    -25.36121065      2.16199208       8       7       2       2       0    1727       0       0       0
O              -5.42469405     25.50230309      1.02571560       8       4       2       1       1    1728       0       0       0
Si            -10.26891859    -19.61454300      1.18665695      14       2       1       1       0    1729       0       0       0
Si             -6.48722454    -21.40999977      1.19042465      14       3       1       1       0    1730       0       0       0
O              -7.86803883    -21.93184933      1.72938368       8       8       1       1       0    1731       0       0       0
O              -7.13220698    -23.57642654     -1.32745307       8       9       1       1       0    1732       0       0       0
O             -10.02918374    -20.16008608     -0.34720939       8       7       1       1       0    1733       0       0       0
Si             -8.37670959    -23.04276912     -2.14970990      14       1       1       1       1    1734       0       0       0
O              -9.47961360    -22.46223372     -1.22949089       8       4       1       1       1    1735       0       0       0
O              -8.97086496    -18.93867874      1.72366967       8       5       1       1       1    1736       0       0       0
O              -6.55822758    -20.75583572     -0.36059901       8       6       1       1       1    1737       0       0       0
O             -11.42054507    -23.73334863      1.10981955       8       9       1       1       0    1738       0       0       0
Si            -12.61671625    -22.99427084      0.23431573      14       1       1       1       1    1739       0       0       0
O             -10.94138247    -20.81618618      2.02291313       8       6       1       1       1    1740       0       0       0
Si             -6.03972596    -19.59067580     -1.34549045      14       2       1       2       0    1741       0       0       0
O              -4.75650107    -18.99141037     -0.80634889       8       5       1       2       1    1742       0       0       0
Si            -10.58298719    -21.23219534     -1.28655584      14       3       1       1       0    1743       0       0       0
O             -11.98993142    -21.85615325     -0.65884204       8       8       1       1       0    1744       0       0       0
O              -5.87325554    -20.16957806      2.01508655       8       7       2       2       0    1745       0       0       0
O              -5.45006634    -22.57263934      1.16209231       8       4       2       2       1    1746       0       0       0
Si            -10.24555115    -14.13329786      1.02947541      14       2       1       1       0    1747       0       0       0
Si             -6.35264306    -15.92681538      1.12463741      14       3       1       1       0    1748       0       0       0
O              -7.81580106    -16.66262561      1.62787439       8       8       1       1       0    1749       0       0       0
O              -7.12260696    -18.27061759     -1.36187370       8       9       1       1       0    1750       0       0       0
O             -10.15357600    -14.77654048     -0.37425650       8       7       1       1       0    1751       0       0       0
Si             -8.36335644    -17.77023924     -2.23969535      14       1       1       1       1    1752       0       0       0
O              -9.57288886    -17.14224522     -1.26384098       8       4       1       1       1    1753       0       0       0
O              -8.91576683    -13.71677770      1.76340097       8       5       1       1       1    1754       0       0       0
O              -6.69767797    -15.43477669     -0.42257415       8       6       1       1       1    1755       0       0       0
O             -11.34341135    -18.42886808      1.08073046       8       9       1       1       0    1756       0       0       0
Si            -12.66172146    -17.72501919      0.10595172      14       1       1       1       1    1757       0       0       0
O             -10.94132180    -15.42123858      2.13486764       8       6       1       1       1    1758       0       0       0
Si             -6.22778807    -14.21360652     -1.37983659      14       2       1       2       0    1759       0       0       0
O              -4.81938759    -13.66670559     -0.67786402       8       5       1       2       1    1760       0       0       0
Si            -10.65083410    -15.92791311     -1.24862804      14       3       1       1       0    1761       0       0       0
O             -11.92701495    -16.66100056     -0.67745386       8       8       1       1       0    1762       0       0       0
O              -5.99571514    -14.83790330      2.14552325       8       7       2       2       0    1763       0       0       0
O              -5.33077792    -17.15889642      1.03192777       8       4       2       2       1    1764       0       0       0
Si            -10.32403593     -8.82284664      1.13473062      14       2       1       1       0    1765       0       0       0
Si             -6.46736283    -10.63672814      1.08389693      14       3       1       1       0    1766       0       0       0
O              -7.88656022    -11.19568136      1.62531396       8       8       1       1       0    1767       0       0       0
O              -7.27127132    -13.11324818     -1.41513117       8       9       1       1       0    1768       0       0       0
O             -10.20100281     -9.51461971     -0.34961512       8       7       1       1       0    1769       0       0       0
Si             -8.34115297    -12.52480432     -2.34297694      14       1       1       1       1    1770       0       0       0
O              -9.52141336    -11.84848911     -1.23877233       8       4       1       1       1    1771       0       0       0
O              -8.95261752     -8.19831420      1.62446566       8       5       1       1       1    1772       0       0       0
O              -6.61401653    -10.16961082     -0.35397659       8       6       1       1       1    1773       0       0       0
O             -11.48528783    -13.02977795      1.12905594       8       9       1       1       0    1774       0       0       0
Si            -12.48793482    -12.53226539      0.07370362      14       1       1       1       1    1775       0       0       0
O             -10.78017197    -10.04973920      1.99018700       8       6       1       1       1    1776       0       0       0
Si             -6.10586022     -8.94420159     -1.30348879      14       2       1       2       0    1777       0       0       0
O              -4.77371236     -8.25359390     -0.72446560       8       5       1       2       1    1778       0       0       0
Si            -10.64659085    -10.72511831     -1.23286286      14       3       1       1       0    1779       0       0       0
O             -12.02166216    -11.17606396     -0.76616004       8       8       1       1       0    1780       0       0       0
O              -5.83081747     -9.52436879      2.15291004       8       7       2       2       0    1781       0       0       0
O              -5.28756992    -11.87804834      1.12124909       8       4       2       2       1    1782       0       0       0
Si            -10.33397546     -3.60623364      1.02375957      14       2       1       1       0    1783       0       0       0
Si             -6.42710761     -5.26698711      1.15940099      14       3       1       1       0    1784       0       0       0
O              -7.75111929     -5.91520110      1.61362405       8       8       1       1       0    1785       0       0       0
O              -7.12553725     -7.61466306     -1.41536688       8       9       1       1       0    1786       0       0       0
O             -10.08193016     -4.21106323     -0.30550163       8       7       1       1       0    1787       0       0       0
Si             -8.33408171     -7.10490705     -2.23640474      14       1       1       1       1    1788       0       0       0
O              -9.58185642     -6.46464486     -1.26169383       8       4       1       1       1    1789       0       0       0
O              -8.85776581     -2.87137622      1.58870314       8       5       1       1       1    1790       0       0       0
O              -6.59045402     -4.70911737     -0.38151134       8       6       1       1       1    1791       0       0       0
O             -11.39687185     -7.75282952      1.06171276       8       9       1       1       0    1792       0       0       0
Si            -12.50215321     -7.05274408      0.18720427      14       1       1       1       1    1793       0       0       0
O             -10.77972294     -4.82886588      2.06956985       8       6       1       1       1    1794       0       0       0
Si             -6.07598508     -3.49979484     -1.24037183      14       2       1       2       0    1795       0       0       0
O              -4.82334003     -2.90205520     -0.81170958       8       5       1       2       1    1796       0       0       0
Si            -10.62270368     -5.31668206     -1.28149797      14       3       1       1       0    1797       0       0       0
O             -11.93793632     -5.91471092     -0.72024568       8       8       1       1       0    1798       0       0       0
O              -5.86370786     -4.04977494      2.05951281       8       7       2       2       0    1799       0       0       0
O              -5.37614588     -6.61994820      1.11336788       8       4       2       2       1    1800       0       0       0""", format="string")

            at.set_cutoff(1.2*bond_length(14, 14))
            at.calc_connect()
            embed = at.bfs_grow_single(2, n=4, nneighb_only=True, min_images_only=True)
            at.add_property('hybrid_mark', HYBRID_NO_MARK)
            at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK

            cluster_options = {
               'cluster_periodic_x': False,
               'cluster_periodic_y': False,
               'cluster_periodic_z': True, 
               'terminate': True, 
               'cluster_allow_modification': True, 
               'randomise_buffer': False
            }

            self.t1 = create_cluster_info_from_mark(at, args_str(cluster_options))
            self.cluster1 = carve_cluster(at, args_str(cluster_options), cluster_info=self.t1)

            cluster_options['keep_whole_silica_tetrahedra'] = True

            self.t2 = create_cluster_info_from_mark(at, args_str(cluster_options))
            self.cluster2 = carve_cluster(at, args_str(cluster_options), cluster_info=self.t2)
      
      def test_cluster1_n(self):
            self.assertEqual(self.cluster1.n, 39)

      def test_cluster_index(self):
            self.assertEqual(sorted(self.cluster1.index),
                             [2, 3, 4, 6, 7, 8, 9, 10, 13, 14, 15, 16, 19, 20, 22, 23, 34, 161, 165, 168, 169, 172,
                              173, 174, 177, 183, 185, 188, 189, 190, 193, 194, 206, 207, 212, 341, 345, 348, 352])

      def test_cluster2_n(self):
            self.assertEqual(self.cluster2.n, 56)

      def test_cluster_2_index(self):
            self.assertEqual(sorted(self.cluster2.index),
                             [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 18, 19, 20, 22, 23, 31, 34, 161, 162,
                              162, 165, 168, 169, 172, 173, 174, 175, 175, 177, 179, 183, 184, 184, 185, 186, 188, 189, 190,
                              193, 194, 202, 206, 207, 211, 212, 341, 345, 348, 352, 355, 355])

      def test_cluster1_composition(self):
            descr = self.t1.str[1,:].stripstrings()
            self.assertEqual([(descr == 'h_active').sum(),
                              (descr == 'term').sum()],
                             [22,17])
            
      def test_cluster2_composition(self):
            descr = self.t2.str[1,:].stripstrings()
            self.assertEqual([(descr == 'h_active').sum(),
                              (descr == 'tetra').sum(),
                              (descr == 'n_cut_bond').sum(),
                              (descr == 'term').sum()],
                             [22,18,1,15])
            

if hasattr(quippy, 'Potential'):

   class TestCluster_ForceEnergy(QuippyTestCase):

      def setUp(self):
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

         system_reseed_rng(2065775975)
         self.pot = Potential('IP SW', param_str=self.xml)

         self.at = supercell(diamond(5.44, 14), 4, 4, 4)
         randomise(self.at.pos, 0.1)
         self.at.set_cutoff(1.2*bond_length(14, 14))
         self.at.calc_connect()

         self.embedlist = Table()
         self.embedlist.append((1,0,0,0))
         self.at.bfs_grow_list(self.embedlist, 2, nneighb_only=True)

         self.at.add_property('hybrid',0)
         self.at.hybrid[self.embedlist.int[1,:]] = 1

         self.bufferlist = self.embedlist.copy()
         self.at.bfs_grow_list(self.bufferlist, 2, nneighb_only=True)

         self.at.add_property('hybrid_mark', 0)
         self.at.hybrid_mark[self.at.hybrid == 1] = HYBRID_ACTIVE_MARK
         for i in self.bufferlist.int[1,:]:
            if not self.at.hybrid_mark[i]:
               self.at.hybrid_mark[i] = HYBRID_BUFFER_MARK

         self.args = {'cluster_calc_connect': True,
                      'randomise_buffer': False,
		      'reduce_n_cut_bonds': False}

         self.cluster_info = create_cluster_info_from_mark(self.at, args_str(self.args))
         self.cluster = carve_cluster(self.at, args_str(self.args), self.cluster_info)

         self.f_ref = FortranArray(
[[  4.49946168e-01,  7.69575734e-02, -2.02228557e-01, -1.40435257e-01,
    4.24383255e-01, -1.10172818e+00,  4.72538822e-01, -4.37414786e-01,
   -7.30100985e-02,  3.62099744e-01, -5.69033700e-01,  1.21494559e+00,
    7.75573582e-01,  2.35090293e-01,  4.18280814e-01, -5.23907757e-01,
   -6.75287389e-01,  1.19407061e-01,  9.00185579e-01, -5.59916356e-01,
   -8.46502946e-01, -6.02617918e-01,  2.27858864e-02,  3.24211513e-01,
   -7.13181479e-01, -2.19450323e-01,  4.40156920e-01, -2.55818217e-01,
   -7.01424357e-02, -1.36489641e+00,  6.97375756e-01,  1.41239494e-01,
    2.18391712e-01, -1.08687214e+00, -1.75728276e-01,  8.15331510e-01,
    4.64435176e-01, -1.01053774e+00, -3.30800450e-01,  2.40093383e-01,
    5.94888174e-01, -4.33133795e-01,  1.31644872e-01,  3.89500254e-01,
   -1.34791342e+00, -6.88306041e-01,  1.74442762e-01,  5.32850814e-01,
    6.92768502e-01,  5.88071977e-01, -7.88287370e-01,  7.86492096e-01,
   -7.24949286e-01,  9.55779993e-01,  6.93859128e-02, -4.12101406e-01,
    2.54612985e-01,  5.86817195e-01,  9.71595420e-01, -1.73161468e-01,
    1.84657043e-01,  1.50828956e-01, -3.12151206e-02,  1.09892018e+00,
    6.57705853e-01, -2.95372970e-01,  3.47177595e-01, -2.44153267e-02,
    8.36940652e-01, -1.16631200e+00,  4.72096037e-01,  9.25841651e-01,
    3.11920903e-01,  7.11262448e-02, -6.46862937e-02,  4.30115572e-01,
   -5.47042811e-01, -1.90791689e-01,  2.42259988e-01, -2.38187987e-01,
    5.12247408e-01,  3.61632622e-01, -4.87853480e-01, -2.30776954e-01,
   -3.13914343e-01, -2.09794901e-01, -1.10468369e-01, -7.78521426e-01,
   -2.38587975e-01, -5.43450639e-01,  3.52774477e-01, -1.66848529e-01,
    1.27154124e-01,  7.50659076e-01,  3.96113023e-01,  2.07323398e-01,
   -4.99075019e-01, -5.29301887e-01,  1.86151056e-01, -3.25059500e-01,
    1.84362024e-01,  5.51116248e-02,  1.85996241e-01, -3.40359972e-02,
   -2.74346584e-01,  1.33740126e-01,  4.96382211e-02,  3.34789740e-02,
   -3.52760459e-01,  6.88380597e-02, -1.33457928e-01, -4.70864915e-01,
   -2.88575581e-01,  9.70300878e-02, -6.88541909e-01, -1.22833217e-01,
   -8.93380864e-02,  3.75380265e-02, -2.97324278e-02, -2.16300046e-01,
   -1.00834572e+00,  4.66077531e-02,  1.81515552e-01,  1.36404906e-01,
   -1.27476063e-01, -1.44888523e-01,  8.24538942e-02, -2.08204693e-01,
    1.91580050e-01, -3.57169381e-01,  1.66993282e-01,  5.39958499e-01,
    1.66087725e-01,  6.04245290e-01,  1.50730337e-01,  5.68478542e-02,
    4.75074464e-01, -2.12528167e-01,  1.83083401e-01,  8.07952881e-01,
   -3.43285135e-02, -4.51512531e-02, -4.79364880e-01,  1.26510238e-01,
    1.54155175e-01,  2.27267819e-01,  1.05440967e+00, -3.24894061e-02,
   -1.30215099e-01,  2.40908005e-01, -1.38489693e-01, -5.99348143e-01,
    2.20770180e-01,  8.26173516e-02, -2.71166577e-01,  1.95148983e-01,
    4.92916653e-01,  3.33011849e-01,  4.76416295e-01, -1.67683048e-02,
    7.25267060e-02, -3.70163113e-01,  5.19814052e-02, -1.24078329e-01,
    9.32507431e-02, -1.73441022e-01,  2.68668772e-02,  1.78673671e-01,
   -1.80281830e-02,  7.08464543e-01,  3.86880290e-02, -1.10629667e-01,
    1.30180922e-01, -4.15479459e-01,  1.90083863e-01,  6.27924227e-01,
   -3.75574161e-01, -3.21365513e-02, -2.08900723e-01, -6.93098620e-01,
   -4.03635968e-01, -3.49420058e-01,  2.71468106e-01, -4.87996009e-01,
   -6.70855970e-01, -3.20492985e-01, -1.78380588e-01, -1.75744458e-01,
   -4.30193585e-01, -5.94052903e-02, -1.72020508e-01],               
 [ -6.09379142e-02,  9.27391444e-01,  4.45722576e-01,  4.06389548e-02,
   -1.04561975e-01,  1.50950893e-01, -6.36056757e-01, -3.31849259e-01,
   -5.90006826e-01,  2.08078307e-02,  2.45540196e-01, -1.72531066e-01,
   -8.20173466e-02,  5.94796741e-01, -6.46599161e-02, -1.74624753e+00,
   -3.55613378e-01,  2.32453017e-01,  1.04235858e-01,  5.94030244e-01,
    8.35712270e-01, -3.48814267e-01, -3.65511939e-01,  5.17756214e-01,
   -1.16385960e+00,  1.38701756e+00,  9.55549645e-01,  3.59684418e-01,
   -7.09091775e-01, -4.28883582e-01, -8.43329007e-01, -1.06212681e+00,
    7.22189726e-01,  1.94871136e-01, -8.42982813e-01, -2.57023848e-01,
    3.33400443e-01,  1.63707232e-01,  9.97798447e-01,  1.36888334e+00,
   -7.62729647e-02,  1.08849355e-01,  1.17890066e-01, -3.49449543e-01,
   -9.58633921e-01,  3.42144076e-01,  5.05608781e-01, -2.79823101e-01,
   -5.82034814e-01,  8.17926530e-01, -3.93984727e-01, -4.50121047e-01,
   -1.67681850e-01,  4.59230408e-02,  1.52193645e+00, -4.93984056e-01,
   -6.20819659e-01, -7.00565979e-02,  4.48323948e-01,  9.17532965e-02,
   -1.76184119e-02, -6.57757868e-02,  2.37541219e-01,  1.74459867e-01,
   -4.94929530e-01,  6.27806961e-01, -3.12370376e-02, -4.07852715e-01,
   -5.60423778e-01, -4.11525544e-01, -4.29649535e-01,  8.59980827e-01,
    9.45124387e-01,  2.52875218e-01, -6.25907546e-01,  5.67113389e-01,
   -3.42420053e-01, -8.12230865e-01,  1.30875125e-01, -9.94368591e-01,
   -3.28611123e-01, -4.86633755e-01,  2.10807254e-01,  1.40314751e-01,
    2.26348091e-01,  3.56708947e-01,  1.83979168e-01,  7.95962008e-01,
    2.40021486e-01, -6.83632654e-03, -3.54644226e-01,  2.77224614e-01,
   -1.10613081e-01, -7.35802969e-01, -3.83646527e-01,  1.11591053e-01,
    9.07626303e-01, -1.19809082e-01,  3.72836921e-01,  9.09212749e-02,
    2.21067501e-01,  4.64412651e-01,  5.57854840e-01,  1.55162271e-02,
   -1.27784544e-01,  2.01732299e-03,  6.18807651e-01,  4.36261061e-01,
    1.11945483e-01,  7.51857087e-01, -2.46538385e-01,  5.55052141e-02,
    2.85050506e-01, -1.85930792e-01,  6.64166636e-01, -8.24235017e-02,
   -9.41726328e-02, -1.76767443e-01,  2.87668514e-01, -1.19350314e-01,
    8.72125483e-01,  5.07211362e-01, -2.19490208e-01, -1.01981818e-01,
   -2.74790718e-01,  6.50053876e-02, -1.41376271e-01,  1.46660698e-01,
   -4.50164300e-01, -5.40399188e-04, -3.26385386e-01, -7.47982048e-01,
   -5.13529158e-01, -1.60388829e-01, -2.36913631e-01, -4.07750773e-01,
   -4.03467269e-01, -2.78933741e-02, -4.32696483e-01, -5.62112050e-01,
    3.91727427e-01, -6.48748518e-01,  8.71384786e-01, -2.47914307e-01,
    3.77508604e-01, -4.25021511e-02,  5.41081760e-01, -2.74312942e-02,
   -5.04326833e-01, -2.29213292e-02, -2.90958694e-01, -1.15830779e+00,
   -1.45976442e-01,  3.75887632e-01,  8.75816152e-01,  2.72256606e-01,
    1.36856994e-01,  2.20268835e-01,  1.84942843e-02,  1.18281776e-01,
    1.14928315e-01, -2.99458818e-01,  6.47391179e-02, -2.13266937e-01,
    1.48889706e-01, -1.02402903e-01,  3.73247416e-02, -8.43169662e-02,
   -5.26347968e-02,  6.21518630e-01,  2.03354977e-01, -9.48504141e-02,
    3.38416421e-02,  6.33399884e-01,  1.62767464e-01,  7.29588225e-01,
    1.52720531e-01, -1.94270111e-02, -3.63315534e-01, -6.34980449e-01,
   -3.96224963e-01, -2.78945887e-01,  3.58028315e-01, -1.52477847e-01,
    1.94015605e-02, -3.14307177e-01, -5.69134308e-01, -1.50456653e-01,
   -5.35582853e-01, -4.81011631e-01, -2.44996983e-01],               
 [  3.92952443e-01, -3.44793646e-01,  1.44150742e-01, -4.75283464e-01,
   -4.38075601e-01, -6.67373519e-01,  6.94790777e-01, -2.21969540e-03,
    7.68595065e-01, -6.70631405e-01, -9.66883026e-01, -3.33663179e-01,
    1.90347279e-02, -1.63420941e-01,  9.52184285e-02, -4.38125533e-01,
   -3.23407833e-01, -3.58151353e-01,  1.34675976e+00,  2.59479331e-01,
   -8.66493832e-01,  2.58248094e-01,  7.93117706e-01,  3.90202931e-02,
   -4.04624381e-01,  3.72107888e-01, -9.31221251e-01, -6.25683017e-01,
   -3.49323442e-01,  2.33099898e-01, -5.02532208e-01,  7.11213820e-01,
   -1.11157712e-01,  5.26227531e-02,  7.32510284e-01, -9.42931522e-01,
    1.19455949e+00, -4.04145547e-01,  3.19829128e-01,  6.66613152e-01,
    8.55262686e-01,  2.78329258e-01,  8.22779921e-01, -2.16964442e-02,
   -4.93437728e-01, -1.29623816e+00, -1.61900413e-01,  5.58001878e-02,
    9.85600448e-01,  2.97387788e-01,  7.90690269e-01, -4.85573210e-01,
    1.41760378e-01,  8.78492978e-01, -5.76106156e-01,  6.30911352e-02,
    6.38857585e-01, -8.59779818e-01,  5.37260041e-01,  8.19806741e-01,
   -5.73349676e-01,  2.98770513e-01,  1.00903802e+00, -1.78720507e-02,
    3.79308337e-01,  8.26671558e-01,  3.75771108e-01, -1.88886798e-01,
    8.43088897e-03, -7.91199079e-01, -2.51055587e-01, -8.38637628e-01,
   -9.28806782e-02,  1.63474509e+00, -5.39939790e-01, -3.64023988e-01,
   -5.29864857e-02, -1.78804225e-02, -7.55154696e-01, -6.14797969e-01,
   -3.97828996e-01,  9.96254933e-02, -5.31655053e-01,  1.64289867e-01,
    2.76841556e-01,  5.77101548e-02,  1.76690732e-02, -8.11655395e-02,
    2.44732031e-01,  5.45035353e-01,  2.65668201e-01,  4.83368807e-01,
    5.84391289e-01,  1.18517208e-01,  7.49150119e-01, -4.31513659e-01,
   -4.52571109e-01,  4.44815241e-01, -2.05236594e-01,  3.10695971e-01,
   -1.34111946e-01, -1.09625013e-01,  1.86855312e-01, -8.64448485e-02,
   -6.90240626e-02, -6.60781998e-02,  9.42334541e-03,  4.01393690e-01,
    1.80306247e-01,  7.62833289e-01, -3.08936611e-01,  1.36440749e-01,
    3.22479153e-01,  1.17990771e-01, -2.57455117e-01,  3.35264227e-01,
    4.81787735e-01, -4.79748535e-03,  3.46432105e-01, -2.46636933e-01,
   -4.84239006e-01,  4.61946815e-01,  3.73361580e-01,  6.37543488e-02,
   -9.37575963e-02,  3.46773658e-01,  2.26129932e-01, -2.54387049e-01,
    1.61006383e-01, -4.91439164e-01,  5.91764261e-02, -4.68954410e-02,
    2.20836776e-01,  5.60604091e-01,  1.71888518e-01,  3.54344853e-01,
    2.65693471e-02, -1.71038021e-01,  3.71009468e-01, -3.17587890e-01,
   -4.24280155e-01,  6.89488131e-01, -9.02908002e-01,  2.41301083e-01,
   -3.29190129e-01,  1.29659744e-01,  4.17330253e-01,  1.92729257e-02,
    4.84665457e-01,  2.37827815e-01, -9.40270185e-03,  6.95324594e-01,
    3.50825166e-01, -2.36605738e-01, -8.33539812e-01, -2.82641015e-01,
   -8.49122557e-02, -3.17346586e-01, -4.33322376e-01,  5.45668157e-02,
   -9.64683468e-02,  3.27406963e-01, -1.64836338e-01,  1.06741594e-01,
   -2.44606215e-01, -2.33641440e-01, -5.68962181e-01, -5.63979149e-01,
   -2.40652871e-01,  1.38102273e-01, -1.40593874e-01,  3.26682805e-02,
   -1.94269510e-01, -1.12266914e+00, -6.03174494e-01,  3.18843387e-02,
   -4.68457854e-01, -3.53587885e-02, -1.27924244e-01,  1.13189169e-01,
   -2.77125291e-01, -8.76382096e-03,  2.69710462e-01, -2.24548323e-01,
   -6.89565112e-01, -3.20368922e-01, -5.11605579e-01, -2.39169788e-01,
   -1.26308958e-01, -4.18427807e-01, -1.41271927e-01]])

         self.e_ref = -566.70771246811125
         self.cluster.set_cutoff(self.pot.cutoff())
         self.cluster.calc_connect()
         self.pot.calc(self.cluster, force=True, energy=True)


      def test_force(self):
         self.assertArrayAlmostEqual(self.cluster.force, self.f_ref)

      def test_energy(self):
         self.assertAlmostEqual(self.cluster.energy, self.e_ref)



if __name__ == '__main__':
   unittest.main()
