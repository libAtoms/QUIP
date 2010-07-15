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
      self.at.set_cutoff_factor(1.2)
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

      self.t = create_cluster_info_from_hybrid_mark(self.at, "terminate=F cluster_allow_modification=F")

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
      self.t = create_cluster_info_from_hybrid_mark(self.at, "terminate=T cluster_allow_modification=F", self.cut_bonds)

      sort_for_comparison(self.embed, self.t)

      self.n_term = self.t.n - self.embed.n
      self.term_atoms =  [122, 508, 510,  26, 124,  30,  98, 126, 100, 410, 508, 512,  26, 412,  32, 386, 416, 388,
                          482, 510, 512,  98, 486, 104, 386, 488, 390,  30,  32,   4, 100, 104,   6, 388, 390,   8]
      self.term_atoms.sort()

      self.cluster = carve_cluster(self.at, cluster_info=self.t, args_str="")
      self.cluster.set_cutoff_factor(1.2)
      self.cluster.calc_connect()

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

      self.t = create_cluster_info_from_hybrid_mark(self.at, "terminate=T cluster_allow_modification=F cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T")
      self.cluster = carve_cluster(self.at, cluster_info=self.t, cluster_periodic_z=True)
      self.cluster.set_cutoff_factor(1.2)
      self.cluster.calc_connect()
      
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
      cluster_info = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
      cluster = carve_cluster(self.at, args_str(self.args), cluster_info)
      self.assert_(cluster.z.sum() % 2 == 0)
      
   def test_15(self):
      # with z[1] = 15, need to remove an H atom
      self.at.z[1] = 15
      cluster_info = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
      cluster = carve_cluster(self.at, args_str(self.args), cluster_info)
      self.assert_(cluster.z.sum() % 2 == 0)

   def test_13(self):
      # with z[1] = 13, need to remove an H atom
      self.at.z[1] = 13
      cluster_info = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
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

      self.t = create_cluster_info_from_hybrid_mark(self.at, "terminate=F cluster_allow_modification=T")

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
      self.assertEqual(str(self.t.str[1,-1]), 'n_cut_bond')

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

      self.t = create_cluster_info_from_hybrid_mark(self.at, "terminate=T reduce_n_cut_bonds=F cluster_allow_modification=T")

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

      self.t = create_cluster_info_from_hybrid_mark(self.at, "terminate=F cluster_allow_modification=F")

      self.r_scale = 0.9

      self.args = {'randomise_buffer' : False}
      self.cluster = carve_cluster(self.at, args_str(self.args), self.t)

      self.args_r_scale =  {'randomise_buffer' : False,
                            'do_rescale_r': True,
                            'r_scale': self.r_scale }
      self.cluster_r_scale = carve_cluster(self.at, args_str(self.args_r_scale), self.t)

   def test_lattice(self):
      # Check non-zero lattice components
      lattice_ratio = self.cluster_r_scale.lattice/self.cluster.lattice
      self.assertAlmostEqual(abs(lattice_ratio[self.cluster.lattice != 0.0] - self.r_scale).max(), 0.0)

   def test_pos(self):
      # Check non-zero position components
      pos_ratio = self.cluster_r_scale.pos/self.cluster.pos
      self.assertAlmostEqual(abs(pos_ratio[self.cluster.pos != 0.0] - self.r_scale).max(), 0.0)


class TestCluster_RandomiseBuffer(QuippyTestCase):

   def setUp(self):
      system_reseed_rng(1)
      
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
      t = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
      self.cluster = carve_cluster(self.at, args_str(self.args), t)

      self.args_randomise_buffer = {'terminate' :True,
                                    'randomise_buffer' : True}
      t = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args_randomise_buffer))
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

      t = create_cluster_info_from_hybrid_mark(at2, args_str(self.args_randomise_buffer))
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

   def test_split_qm_raises_runtime_error(self):
      # Split QM region by marking another atom
      self.at.hybrid_mark[107] = HYBRID_ACTIVE_MARK
      verbosity_push(PRINT_SILENT)
      self.assertRaises(RuntimeError, create_cluster_info_from_hybrid_mark, self.at, args_str(self.args))
      verbosity_pop()
  

class TestCluster_Surface_Dia(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 3, 3, 3)
      lat = self.at.lattice
      lat[3,3] *= 2
      self.at.set_lattice(lat, scale_positions=False)
      self.at.calc_connect()

      self.at.add_property('hybrid_mark', HYBRID_NO_MARK)

      embed = self.at.bfs_grow_single(1, n=3)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_hybrid_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T")
      self.cluster3 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T", cluster_info=self.t)

      embed = self.at.bfs_grow_single(1, n=4)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_hybrid_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=T")
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
      self.t = create_cluster_info_from_hybrid_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F cluster_nneighb_only=F")
      self.cluster1 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F randomise_buffer=F", cluster_info=self.t)

      embed = self.at.bfs_grow_single(1, n=2, nneighb_only=False)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_hybrid_mark(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F cluster_nneighb_only=F")
      self.cluster2 = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F randomise_buffer=F", cluster_info=self.t)

      hollow_atom = 28
      embed.delete(hollow_atom)
      self.at.hybrid_mark[:] = HYBRID_NO_MARK
      self.at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK
      self.t = create_cluster_info_from_hybrid_mark(self.at, "cluster_allow_modification=T reduce_n_cut_bonds=T cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F cluster_nneighb_only=F")
      self.cluster2h = carve_cluster(self.at, "cluster_periodic_x=T cluster_periodic_y=T cluster_periodic_z=T terminate=F randomise_buffer=F", cluster_info=self.t)

      # self.at.print_xyz("tmp_at.xyz", all_properties=True)
      # self.cluster1.print_xyz("tmp_cluster1.xyz", all_properties=True)
      # self.cluster2.print_xyz("tmp_cluster2.xyz", all_properties=True)
      # self.cluster2h.print_xyz("tmp_cluster2h.xyz", all_properties=True)

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

      self.t1 = create_cluster_info_from_hybrid_mark(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T cluster_nneighb_only=F cluster_allow_modification=F")
      self.cluster1 = carve_cluster(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T randomise_buffer=F", cluster_info=self.t1)

      self.t2 = create_cluster_info_from_hybrid_mark(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T cluster_nneighb_only=T cluster_allow_modification=F")
      self.cluster2 = carve_cluster(self.at, "cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z=T terminate=T randomise_buffer=F", cluster_info=self.t2)

      #AtomsList([self.at,self.cluster1,self.cluster2]).show()
      #raw_input()

   def test_cluster_n(self):
      self.assertEqual(self.cluster1.n, 97)

   def test_cluster_n_term(self):
      self.assertEqual((self.cluster1.z == 1).count(), 40)

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

      self.at.set_cutoff_factor(1.2, 1.4)
      self.at.calc_connect()
      self.at.calc_connect_hysteretic(self.at.hysteretic_connect)

      mark_atoms(self.at)

      self.args = {'terminate': True,
                   'randomise_buffer' : False,
                   'hysteretic_connect': False}
      self.t = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
      self.cluster = carve_cluster(self.at, args_str(self.args), cluster_info=self.t)

      self.args_hysteretic = {'terminate' :True,
                              'randomise_buffer' : False,
                              'hysteretic_connect': True,
                              'cluster_nneighb_only': False}
      self.t_hysteretic = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args_hysteretic))
      self.cluster_hysteretic = carve_cluster(self.at, args_str(self.args_hysteretic), cluster_info=self.t_hysteretic)

   def test_cluster_info(self):
      self.assertEqual(self.t, self.t_hysteretic)

   def test_cluster(self):
      self.assertEqual(self.cluster, self.cluster_hysteretic)


class TestCluster_HystereticConnect_MoveAtom(QuippyTestCase):

   def setUp(self):
      self.at = supercell(diamond(5.44, 14), 3, 3, 3)
      self.at.set_cutoff_factor(1.2)
      self.at.calc_connect()
      mark_atoms(self.at)

      self.move_atom = 66

      self.at_hysteretic = self.at.copy()
      self.at_hysteretic.set_cutoff_factor(1.2, 1.4)
      self.at_hysteretic.calc_connect_hysteretic(self.at_hysteretic.hysteretic_connect)

      self.at_move = self.at.copy()
      
      self.at_move.pos[2,self.move_atom] -= 0.5
      self.at_hysteretic.pos[2,self.move_atom] -= 0.5

      self.at_move.calc_connect()
      mark_atoms(self.at_move)

      self.args = {'terminate': True,
                   'randomise_buffer' : False,
                   'hysteretic_connect': False}

      self.t = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
      self.cluster = carve_cluster(self.at, args_str(self.args), cluster_info=self.t)

      self.t_move = create_cluster_info_from_hybrid_mark(self.at_move, args_str(self.args))
      self.cluster_move = carve_cluster(self.at_move, args_str(self.args), cluster_info=self.t)

      self.at_hysteretic.calc_connect_hysteretic(self.at_hysteretic.hysteretic_connect)
      mark_atoms(self.at_hysteretic, nneighb_only=False, alt_connect=self.at_hysteretic.hysteretic_connect)

      self.args_hysteretic = {'terminate' :True,
                              'randomise_buffer' : False,
                              'hysteretic_connect': True,
                              'cluster_nneighb_only': False}
      self.t_hysteretic = create_cluster_info_from_hybrid_mark(self.at_hysteretic, args_str(self.args_hysteretic))
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
            from quippy.surface import crack_rotation_matrix, orthorhombic_slab

            aq = alpha_quartz(**sio2.quartz_params['ASAP_JRK'])
            unit_slab = orthorhombic_slab(aq, rot=crack_rotation_matrix(aq,(0,0,0,1),z=(-1,0,0)),verbose=False)
            slab = supercell(unit_slab, 10, 10, 1)

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

            at = slab.select(mask)
            at.calc_connect()

            embed = at.bfs_grow_single(2, n=4, nneighb_only=True)
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

            self.t1 = create_cluster_info_from_hybrid_mark(at, args_str(cluster_options))
            self.cluster1 = carve_cluster(at, args_str(cluster_options), cluster_info=self.t1)

            cluster_options['keep_whole_silica_tetrahedra'] = True

            self.t2 = create_cluster_info_from_hybrid_mark(at, args_str(cluster_options))
            self.cluster2 = carve_cluster(at, args_str(cluster_options), cluster_info=self.t2)

      
      def test_cluster1_n(self):
            self.assertEqual(self.cluster1.n, 39)

      def test_cluster_index(self):
            self.assertEqual(sorted(self.cluster1.index),
                             [2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 18, 19, 21, 22, 33, 160, 164, 167, 168, 171, 
                              172, 173, 176, 182, 184, 187, 188, 189, 192, 193, 205, 206, 211, 340, 344, 347, 351])

      def test_cluster2_n(self):
            self.assertEqual(self.cluster2.n, 56)

      def test_cluster_2_index(self):
            self.assertEqual(sorted(self.cluster2.index),
                             [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 17, 18, 19, 21, 22, 30, 33, 160, 161, 
                              161, 164, 167, 168, 171, 172, 173, 174, 174, 176, 178, 182, 183, 183, 184, 185, 187, 188, 
                              189, 192, 193, 201, 205, 206, 210, 211, 340, 344, 347, 351, 354, 354])

      def test_cluster1_composition(self):
            descr = self.t1.str[1,:].stripstrings()
            self.assertEqual([(descr == 'h_active').count(),
                              (descr == 'term').count()],
                             [22,17])
            
      def test_cluster2_composition(self):
            descr = self.t2.str[1,:].stripstrings()
            self.assertEqual([(descr == 'h_active').count(),
                              (descr == 'tetra').count(),
                              (descr == 'n_cut_bond').count(),
                              (descr == 'term').count()],
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

         system_reseed_rng(1)
         self.pot = Potential('IP SW', self.xml)

         self.at = supercell(diamond(5.44, 14), 4, 4, 4)
         matrix_randomise(self.at.pos, 0.1)
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

         self.cluster_info = create_cluster_info_from_hybrid_mark(self.at, args_str(self.args))
         self.cluster = carve_cluster(self.at, args_str(self.args), self.cluster_info)

         self.f_ref = FortranArray([[  4.49946168e-01,  -6.09379142e-02,   3.92952443e-01],
                                    [  7.69575734e-02,   9.27391444e-01,  -3.44793646e-01],
                                    [ -2.02228557e-01,   4.45722576e-01,   1.44150742e-01],
                                    [ -1.40435257e-01,   4.06389548e-02,  -4.75283464e-01],
                                    [  4.24383255e-01,  -1.04561975e-01,  -4.38075601e-01],
                                    [ -1.10172818e+00,   1.50950893e-01,  -6.67373519e-01],
                                    [  4.72538822e-01,  -6.36056757e-01,   6.94790777e-01],
                                    [ -4.37414786e-01,  -3.31849259e-01,  -2.21969540e-03],
                                    [ -7.30100985e-02,  -5.90006826e-01,   7.68595065e-01],
                                    [  3.62099744e-01,   2.08078307e-02,  -6.70631405e-01],
                                    [ -5.69033700e-01,   2.45540196e-01,  -9.66883026e-01],
                                    [  1.21494559e+00,  -1.72531066e-01,  -3.33663179e-01],
                                    [  7.75573582e-01,  -8.20173466e-02,   1.90347279e-02],
                                    [  2.35090293e-01,   5.94796741e-01,  -1.63420941e-01],
                                    [  4.18280814e-01,  -6.46599161e-02,   9.52184285e-02],
                                    [ -5.23907757e-01,  -1.74624753e+00,  -4.38125533e-01],
                                    [ -6.75287389e-01,  -3.55613378e-01,  -3.23407833e-01],
                                    [  1.16212994e-01,   2.33807426e-01,  -3.57987926e-01],
                                    [  9.00865980e-01,   1.02648958e-01,   1.34640116e+00],
                                    [ -5.61232766e-01,   5.94375079e-01,   2.58158081e-01],
                                    [ -8.45451051e-01,   8.35795673e-01,  -8.66315318e-01],
                                    [ -6.03214357e-01,  -3.48598437e-01,   2.59936659e-01],
                                    [  2.31133176e-02,  -3.65238052e-01,   7.94697994e-01],
                                    [  3.25241203e-01,   5.17339532e-01,   3.81696098e-02],
                                    [ -7.14249295e-01,  -1.16172581e+00,  -4.04560675e-01],
                                    [ -2.20298120e-01,   1.39035805e+00,   3.71767276e-01],
                                    [  4.41216139e-01,   9.55317131e-01,  -9.30998955e-01],
                                    [ -2.54206800e-01,   3.59321548e-01,  -6.26018673e-01],
                                    [ -6.96736644e-02,  -7.11004908e-01,  -3.48846847e-01],
                                    [ -1.36516934e+00,  -4.29532278e-01,   2.31083169e-01],
                                    [  6.97770531e-01,  -8.41390678e-01,  -5.02275486e-01],
                                    [  1.41386696e-01,  -1.06289729e+00,   7.13051235e-01],
                                    [  2.16032214e-01,   7.20638556e-01,  -1.09041574e-01],
                                    [ -1.08633723e+00,   1.95476000e-01,   5.03984300e-02],
                                    [ -1.74089762e-01,  -8.42900057e-01,   7.31890148e-01],
                                    [  8.16468165e-01,  -2.57161798e-01,  -9.42833365e-01],
                                    [  4.64668903e-01,   3.34698812e-01,   1.19421758e+00],
                                    [ -1.01092085e+00,   1.65574341e-01,  -4.04852770e-01],
                                    [ -3.33394643e-01,   9.96791104e-01,   3.19633559e-01],
                                    [  2.39140333e-01,   1.36827390e+00,   6.63554687e-01],
                                    [  5.94345388e-01,  -7.91705824e-02,   8.54480237e-01],
                                    [  1.47508485e-01,  -1.12399858e-01,   1.56800416e-01],
                                    [  6.53124727e-01,  -3.97324588e-01,   1.02858643e+00],
                                    [  8.24748414e-01,  -1.37549447e-01,  -4.62011429e-01],
                                    [ -1.06248544e+00,  -5.07628221e-01,  -6.69374435e-01],
                                    [ -1.07842657e+00,   7.35733209e-01,  -1.02314687e+00],
                                    [  7.65194291e-01,  -4.31216941e-03,   3.65084425e-01],
                                    [  5.29299880e-01,  -3.32875306e-02,   5.71925207e-02],
                                    [  9.97170843e-02,  -9.35899127e-01,   7.25971820e-01],
                                    [  5.11146119e-01,   6.45426932e-01,   1.27124290e-01],
                                    [ -7.16063471e-01,  -6.15202701e-01,   2.36629412e-01],
                                    [  4.99978360e-01,  -7.93199816e-01,  -8.18023604e-01],
                                    [ -4.59200973e-01,   1.47754149e-01,  -2.95829785e-01],
                                    [  8.77801932e-01,  -3.59533558e-01,   8.65156721e-01],
                                    [  1.11028081e-01,   1.01684690e+00,  -7.82088237e-01],
                                    [  1.80066667e-01,  -1.30602298e+00,   4.36386885e-01],
                                    [ -1.79257993e-01,  -4.60062135e-01,   5.53519253e-01],
                                    [  3.66863224e-01,   1.42357648e-01,  -7.55336911e-01],
                                    [  7.15913181e-01,   3.33915227e-01,   2.81834595e-01],
                                    [  2.04168319e-02,  -1.69610299e-01,   9.89817342e-01],
                                    [ -6.09305038e-01,   4.25196921e-01,  -3.70542517e-01],
                                    [ -1.21186196e-01,   6.54768293e-01,   1.61626847e-01],
                                    [  8.73419064e-02,   5.18958050e-01,   7.14043996e-01],
                                    [  6.61965457e-01,   7.61639411e-01,   4.25145129e-01],
                                    [  9.97457305e-01,  -6.16701687e-01,   5.14112990e-01],
                                    [ -8.48897570e-01,  -9.86977861e-03,   3.62741199e-01],
                                    [  2.22369986e-01,   2.08626067e-01,  -2.54993844e-01],
                                    [  4.38129418e-01,   2.86866878e-01,  -9.45716713e-01],
                                    [  7.91295349e-01,  -9.14754670e-01,   7.06984405e-01],
                                    [ -1.04577861e+00,  -6.54573933e-01,  -5.43965919e-01],
                                    [  1.39609425e-01,  -2.76304188e-01,   8.58081302e-02],
                                    [  4.53701905e-01,   7.54385234e-01,  -6.85468337e-01],
                                    [  6.69228527e-02,   7.00598876e-01,  -2.00854111e-01],
                                    [  1.16496459e-01,   2.99815846e-01,   1.42074550e+00],
                                    [ -7.26572069e-01,  -4.30238209e-01,  -2.05428845e-01],
                                    [  2.69605090e-01,   3.89939328e-02,  -2.31329989e-01],
                                    [ -2.46687266e-01,  -1.20256629e+00,   5.20644922e-01],
                                    [  1.56794373e-01,  -9.53040254e-01,   3.11878524e-01],
                                    [  9.33507508e-01,   4.26841978e-01,  -6.48622018e-01],
                                    [  1.20446227e-01,  -6.21457480e-01,  -7.69514188e-01],
                                    [  8.50756621e-01,  -2.55211979e-01,   3.47917616e-01],
                                    [  1.85877508e-01,  -7.69425227e-02,   5.07382039e-01],
                                    [ -3.17463710e-01,   9.47075403e-01,  -2.78355026e-01],
                                    [ -2.30241453e-01,   1.40891014e-01,   1.65651697e-01],
                                    [  4.18544278e-02,  -1.37782459e-01,   1.05387981e-01],
                                    [  1.47281023e-01,   1.94468019e-01,  -2.83733759e-01],
                                    [ -1.09366989e-01,   1.83049412e-01,   2.29903664e-02],
                                    [ -2.76739117e-01,   2.91226749e-01,  -3.10822919e-01],
                                    [ -2.37867238e-01,   2.42797751e-01,   2.43990716e-01],
                                    [ -2.32969857e-01,  -1.24974088e-01,   2.35939602e-01],
                                    [  1.98313547e-02,  -7.37231188e-03,   1.23468520e-01],
                                    [  2.56237205e-01,  -1.45255521e-01,   2.97490519e-01],
                                    [  2.99137533e-01,   1.99465504e-01,  -2.26007483e-01],
                                    [  2.06481977e-01,  -1.98168119e-01,  -9.46643277e-02],
                                    [  6.17503562e-01,  -5.97314712e-01,   6.44576251e-01],
                                    [ -5.62598722e-02,   9.76904376e-04,  -1.57937370e-01],
                                    [ -2.42280958e-01,   2.36697196e-01,  -3.18194607e-01],
                                    [ -2.81172967e-01,  -2.17291420e-01,   2.05360187e-01],
                                    [ -2.60423658e-01,   1.50854754e-01,   2.46604884e-01],
                                    [  1.53855965e-02,  -6.38126426e-02,  -3.42686077e-02],
                                    [ -4.20561134e-02,   1.20541235e-01,   8.30465868e-02],
                                    [  3.09456585e-01,   3.52639458e-01,  -3.61787415e-01],
                                    [  7.09162293e-02,   4.47336664e-02,  -4.43274025e-03],
                                    [ -3.44690378e-02,   1.53973910e-02,  -8.65738548e-02],
                                    [ -2.06508390e-01,  -2.78317036e-01,  -2.20405046e-01],
                                    [  1.32027425e-01,   2.22889847e-03,  -6.67310386e-02],
                                    [  3.87019643e-01,   4.61142807e-01,  -3.30992821e-01],
                                    [  1.49669705e-01,   1.48885108e-01,   1.15295482e-01],
                                    [ -2.78826083e-01,   2.85703136e-01,   3.51891737e-01],
                                    [  2.38709970e-01,   3.24305913e-01,   3.37567197e-01],
                                    [  1.58407720e-02,   8.83004670e-02,   4.18145655e-02],
                                    [  1.59750375e-01,   1.92604418e-01,  -7.12331242e-02],
                                    [ -9.96237904e-02,  -1.14903106e-01,  -8.87348041e-02],
                                    [ -8.90819367e-02,   1.04481854e-02,   3.36999138e-02],
                                    [ -4.34295935e-01,   4.06906064e-01,  -3.49432999e-01],
                                    [ -3.77794264e-01,  -2.43416678e-01,  -3.82285429e-01],
                                    [  2.00354110e-01,  -3.75600383e-01,   3.58856864e-01],
                                    [  3.63166344e-02,  -1.77376609e-01,  -4.66169670e-03],
                                    [  7.14002772e-02,   2.68308740e-02,   9.93218019e-02],
                                    [ -7.62631949e-02,   1.98447723e-01,   5.78624502e-02],
                                    [ -4.71095106e-01,   5.31929174e-01,  -4.92121110e-01],
                                    [  1.66861102e-01,   2.17884253e-01,   1.71237409e-01],
                                    [ -2.62254578e-01,   2.38394513e-01,   1.66049676e-01],
                                    [  1.36078350e-01,  -1.02311315e-01,   6.45373264e-02],
                                    [ -2.63617969e-01,  -2.11260197e-01,  -2.34455025e-01],
                                    [ -2.36017691e-01,   1.55667446e-01,   3.07520213e-01],
                                    [  8.22589645e-02,  -1.41214343e-01,   2.26725944e-01],
                                    [ -2.08491510e-01,   1.45810795e-01,  -2.54667149e-01],
                                    [ -2.07044776e-01,  -2.69556750e-01,  -2.19343055e-01],
                                    [ -9.23161020e-02,   1.07409489e-01,  -2.27460582e-01],
                                    [ -1.86869672e-01,  -1.81145047e-01,  -2.76345270e-01],
                                    [  1.09470747e-01,  -6.48010046e-02,  -5.84977647e-02],
                                    [  3.57836456e-01,  -4.29852517e-01,   4.12935392e-01],
                                    [  1.17227172e-01,   4.91319923e-02,   7.78943197e-02],
                                    [  1.51519075e-01,  -2.36003235e-01,   1.74184572e-01],
                                    [ -7.88325801e-02,  -7.16163337e-02,   1.07593165e-02],
                                    [  1.74897367e-01,  -1.17514273e-01,  -1.02430620e-01],
                                    [ -2.11085682e-01,  -2.76724048e-02,  -1.71422049e-01],
                                    [ -2.49433324e-02,   4.67375209e-02,  -1.15136473e-01],
                                    [  2.06245867e-01,  -2.62343149e-01,  -2.57770980e-01],
                                    [ -1.72070446e-01,   7.67104234e-02,  -9.94600006e-02],
                                    [ -2.42605989e-01,  -1.70987110e-01,   2.05681105e-01],
                                    [ -5.84595808e-01,   5.75890995e-01,  -5.93645878e-01],
                                    [  1.35091141e-02,   1.72851528e-03,  -1.22101483e-02],
                                    [  7.75629058e-03,   1.09443455e-02,   3.90964333e-02],
                                    [  1.19462948e-01,   2.11869614e-01,  -1.30654694e-01],
                                    [  2.12782287e-01,   2.27075281e-01,   2.02800224e-01],
                                    [ -3.06046096e-02,  -2.80149129e-02,   1.84235157e-02],
                                    [ -2.56722558e-01,  -2.11734508e-01,   1.88963287e-01],
                                    [  4.80053995e-02,   5.77646796e-02,   3.72979641e-02],
                                    [  1.35435559e-01,  -1.72599696e-01,   2.65403700e-01],
                                    [ -3.67900175e-01,  -2.75734814e-01,   4.14150857e-01],
                                    [ -5.92073489e-02,  -2.30247982e-02,   6.54465425e-02],
                                    [ -3.05333975e-01,   2.03842393e-01,   1.63977128e-01],
                                    [ -4.32347524e-01,   4.60589650e-01,  -4.14536959e-01],
                                    [  1.98780891e-01,   2.72206656e-01,  -2.81625547e-01],
                                    [  3.40750979e-01,  -2.14488739e-01,   2.59752188e-01],
                                    [  3.33736587e-01,   2.20047831e-01,  -3.17140113e-01],
                                    [  1.72732540e-01,  -1.00164610e-01,  -1.28916770e-01],
                                    [ -1.72341217e-02,   1.19697206e-01,   5.50582046e-02],
                                    [ -8.90958288e-02,   4.39367260e-02,   6.83712479e-02],
                                    [ -3.70647333e-01,  -2.98659286e-01,   3.26884454e-01],
                                    [ -1.54324275e-01,  -1.40420775e-01,  -7.41156879e-02],
                                    [ -1.24330781e-01,  -2.13525561e-01,   1.05875348e-01],
                                    [ -1.59637078e-01,  -1.07013259e-01,  -1.29860494e-01],
                                    [  7.65315561e-02,   1.47094214e-01,  -1.25918021e-01],
                                    [ -3.38247008e-01,  -3.42711549e-01,  -4.06718173e-01],
                                    [  6.42402133e-03,   8.25113283e-02,   4.42312453e-02],
                                    [  2.05676384e-01,   1.67520025e-01,  -1.47195804e-01],
                                    [  3.39862329e-01,   2.62006941e-01,   2.80312138e-01],
                                    [ -5.23524606e-02,  -8.14812687e-03,   6.80591065e-02],
                                    [ -1.09879514e-01,  -9.54292646e-02,   3.18533451e-02],
                                    [ -1.63164724e-02,  -1.11446566e-01,  -1.29839273e-01],
                                    [ -2.24431320e-01,   2.95477606e-01,  -2.74596811e-01],
                                    [  4.62479927e-01,   4.43113619e-01,  -4.79537978e-01],
                                    [  1.94521976e-01,   2.91250104e-01,   2.02284245e-01],
                                    [ -1.72461992e-01,   2.36474835e-01,  -2.59757355e-01],
                                    [ -3.16481383e-02,  -2.13198146e-02,  -3.49479982e-02],
                                    [  1.41198594e-01,  -2.08969901e-01,   2.07187910e-01],
                                    [ -2.62482851e-01,  -2.11509243e-01,   2.96765360e-01],
                                    [ -4.02854497e-01,  -3.96599759e-01,  -2.79887626e-01],
                                    [ -1.39483943e-01,  -6.72014510e-02,   8.46528938e-02],
                                    [  2.71970145e-01,   3.58695603e-01,   2.67464304e-01],
                                    [ -2.85913789e-01,   2.86365208e-01,   2.24989304e-01],
                                    [ -2.33590689e-01,   2.05218609e-01,  -2.47401387e-01],
                                    [ -3.20659877e-01,  -3.17468813e-01,  -3.19317587e-01],
                                    [  4.99770730e-02,  -5.97357086e-02,  -1.24110262e-02],
                                    [ -1.78347777e-01,  -1.49800869e-01,  -2.38577440e-01],
                                    [ -2.56104913e-04,  -9.50571990e-02,   6.44812326e-02],
                                    [  8.04232110e-02,  -1.42042393e-01,  -7.58670198e-02],
                                    [ -1.73947053e-01,  -2.44126270e-01,  -1.42540079e-01]])

         self.e_ref = -562.993702807
         self.pot.calc(self.cluster, calc_force=True, calc_energy=True)


      def test_force(self):
         self.assertArrayAlmostEqual(self.cluster.force, self.f_ref)

      def test_energy(self):
         self.assertAlmostEqual(self.cluster.energy, self.e_ref)



if __name__ == '__main__':
   unittest.main()
