from quippy import *
import unittest, quippy
from quippytest import *

test_electrostatic_embed = False

if test_electrostatic_embed:

  class TestDimer(QuippyTestCase):

    def setUp(self):
       self.pot = Potential('IP TS', param_str=quip_xml_parameters('ASAP', 'screened_LDA'))
       self.pot.print_()
       self.at = Atoms(n=2, lattice=numpy.diag([10.0, 10.0, 5.0]))
       half_cell = numpy.diag(self.at.lattice)/2.0
       self.at.pos[1] = [0.0,0.0,0.0] + half_cell
       self.at.pos[2] = [3.042*BOHR, 0.0, 0.0] + half_cell
       self.at.set_atoms([14, 8])
       self.at.set_cutoff(self.pot.cutoff())
       self.at.add_property('hybrid_mark', HYBRID_NO_MARK)
       self.at.hybrid_mark[1] = HYBRID_ACTIVE_MARK
       self.at.hybrid_mark[2] = HYBRID_ELECTROSTATIC_MARK
       self.at.calc_connect()
       self.pot.calc(self.at, force=True)
       self.ngrid = (50,50,10)

    def test_calc_esp(self):
       real_grid = fzeros((3,self.ngrid[0]*self.ngrid[1]*self.ngrid[2]))
       pot = fzeros(self.ngrid)
       self.pot.calc_electrostatic_potential(self.at, "hybrid_mark", self.ngrid, self.at.pos[:,1], numpy.diag(self.at.lattice),
                                             real_grid, pot, args_str="calc_short_range=F calc_sc_dipoles=F calc_dipoles=F pseudise=T")
       self.at.write('dimer.full.cube', data=pot/HARTREE)
       write_electrostatic_potential_cube(self.at, 'hybrid_mark', 'dimer.pos.cube', self.ngrid, (0., 0., 0.), numpy.diag(self.at.lattice), pot, write_efield=False)
       write_electrostatic_potential_cube(self.at, 'hybrid_mark', 'dimer.efield.cube', self.ngrid, (0., 0., 0.), numpy.diag(self.at.lattice), pot, write_efield=True)

       #make_periodic_potential(self.at, real_grid, self.ngrid, [False, False, True], 4.0, 1.0, 'hybrid_mark', pot)
       #write_electrostatic_potential(self.at, 'hybrid_mark', 'dimer.periodic.esp', self.ngrid, numpy.diag(self.at.lattice), pot)

    def test_calc_esp_no_pseudise(self):
       real_grid = fzeros((3,self.ngrid[0]*self.ngrid[1]*self.ngrid[2]))
       pot = fzeros(self.ngrid)
       self.pot.calc_electrostatic_potential(self.at, "hybrid_mark", self.ngrid, self.at.pos[:,1], numpy.diag(self.at.lattice),
                                             real_grid, pot, args_str="calc_short_range=F calc_sc_dipoles=F calc_dipoles=F pseudise=F")
       write_electrostatic_potential(self.at, 'hybrid_mark', 'dimer.no-pseudise.esp', self.ngrid, numpy.diag(self.at.lattice), pot)


  class TestCluster(QuippyTestCase):

    def setUp(self):
       self.pot = Potential('IP TS', param_str=quip_xml_parameters('ASAP', 'screened_LDA'))

    def test_cluster(self):
       from quippy.structure_tools import rotation_matrix, orthorhombic_slab

       aq = alpha_quartz(**sio2.quartz_params['ASAP_JRK'])
       unit_slab = orthorhombic_slab(aq, rot=rotation_matrix(aq,(0,0,0,1),z=(-1,0,0)),verbose=False)
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
       at.set_cutoff(self.pot.cutoff())
       at.calc_connect()

       embed = at.bfs_grow_single(2, n=4, nneighb_only=True)
       at.add_property('hybrid_mark', HYBRID_NO_MARK)
       buffer = embed.copy()
       at.bfs_grow_list(buffer, n=2, nneighb_only=True)
       at.hybrid_mark[buffer.int[1,:]] = HYBRID_BUFFER_MARK
       at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK

       cluster_options = {
          'cluster_periodic_x': False,
          'cluster_periodic_y': False,
          'cluster_periodic_z': True, 
          'terminate': True, 
          'cluster_allow_modification': True, 
          'randomise_buffer': False,
          'keep_whole_silica_tetrahedra': True,
          'map_into_cell': False
          }

       cluster_info = create_cluster_info_from_mark(at, args_str(cluster_options))
       cluster = carve_cluster(at, args_str(cluster_options), cluster_info=cluster_info)

       #electrostatic_mark = Table(4,0,0,0)
       #electrostatic_mark.append(cluster_info.int[1:4,:])  
       #at.bfs_grow_list(electrostatic_mark, n=2, nneighb_only=True)
       at.add_property('es_mark', HYBRID_ELECTROSTATIC_MARK)
       #at.es_mark[electrostatic_mark.int[1,:]] = HYBRID_ELECTROSTATIC_MARK
       at.es_mark[at.modified_hybrid_mark != HYBRID_NO_MARK] = HYBRID_NO_MARK

       # centre of cluster
       origin = (at.pos[:,at.modified_hybrid_mark != 0]).mean(axis=2)
       origin[3] = 0.0  # no shift in z-direction

       at.add_property('dummy_cluster_mark', HYBRID_NO_MARK)
       at.dummy_cluster_mark[at.modified_hybrid_mark != HYBRID_NO_MARK] = HYBRID_ACTIVE_MARK
       at.dummy_cluster_mark[at.es_mark != HYBRID_NO_MARK] = HYBRID_ACTIVE_MARK

       # bigger cluster, including the electrostatically embedded atoms
       cluster_options['terminate'] = False
       cluster_options['keep_whole_silica_tetrahedra'] = False
       cluster_info2 = create_cluster_info_from_mark(at, args_str(cluster_options), mark_name='dummy_cluster_mark')
       cluster2 = carve_cluster(at, args_str(cluster_options), cluster_info=cluster_info2)

       cluster2.pos[:] = cluster2.pos - numpy.tile(origin, [cluster2.n, 1]).T

       cluster2.set_cutoff_factor(1.2)
       cluster2.calc_connect()
       cluster2.coalesce_in_one_periodic_image(is_periodic=[False,False,True])

       # Put fractional coordinate (.5,.5,.5) at centre of smaller cell
       cluster2.pos[:] = cluster2.pos + numpy.tile(numpy.diag(cluster.lattice)/2.0, [cluster2.n, 1]).T
       cluster2.write('cluster2.xyz')

       self.pot.calc(at, force=True)
       ngrid = [50, 50, 10]
       real_grid = fzeros((3,ngrid[0]*ngrid[1]*ngrid[2]))
       pot = fzeros(ngrid)
       self.pot.calc_electrostatic_potential(at, "es_mark", ngrid, origin, numpy.diag(cluster.lattice),
                                                  real_grid, pot, args_str="calc_short_range=F calc_sc_dipoles=F calc_dipoles=F pseudise=T")

       write_electrostatic_potential(cluster, 'cluster.esp', ngrid, numpy.diag(cluster.lattice), pot)
       make_periodic_potential(cluster, real_grid, ngrid, [False, False, True], 4.0, 1.0, 'hybrid_mark', pot)
       write_electrostatic_potential(cluster, 'cluster.periodic.esp', ngrid, numpy.diag(cluster.lattice), pot)

       # Translate cluster to centre positions on origin of electrostatic potential
       cluster.pos[:] = cluster.pos - numpy.tile(origin, [cluster.n, 1]).T

       cluster.set_cutoff_factor(1.2)
       cluster.calc_connect()
       cluster.coalesce_in_one_periodic_image(is_periodic=[False,False,True])

       # Put fractional coordinate (.5,.5,.5) at centre of cell
       cluster.pos[:] = cluster.pos + numpy.tile(numpy.diag(cluster.lattice)/2.0, [cluster.n, 1]).T


  
if __name__ == '__main__':
   unittest.main()
