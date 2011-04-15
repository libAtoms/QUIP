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

import os
import numpy as np
from quippy import (Potential, Atoms, MPI_context, transform, print_title, verbosity_push,
                    HYBRID_NO_MARK)

from quippy import (crack_apply_load_increment,
                    crack_find_tip_local_energy, crack_k_to_g,
                    crack_setup_marks, crack_apply_strain_ramp,
                    crack_find_tip_percolation, crack_make_seed,
                    crack_strain_to_g, crack_calc_load_field,
                    crack_g_to_k, crack_make_slab, crack_uniform_load,
                    crack_check_coordination, crack_g_to_strain,
                    crack_measure_g, crack_update_connect,
                    crack_check_coordination_boundaries,
                    crack_is_edge_atom, crack_parse_name,
                    crack_update_selection, crack_find_tip,
                    crack_is_topbottom_edge_atom, crack_print_cio,
                    crack_update_selection_coordination,
                    crack_find_tip_coordination, crack_k_field,
                    crack_print_filename,
                    crack_update_selection_crack_front)

def makecrack(params, stem):
   """Given a CrackParams object `param`, construct and return a new crack slab Atoms object."""

   xmlfilename = stem+'.xml'

   print_title('Initialisation')

   verbosity_push(params.io_verbosity)
   params.print_()

   print("Initialising classical potential with args " + params.classical_args.strip() +
         " from file " + xmlfilename)
   classicalpot = Potential(params.classical_args, param_filename=xmlfilename)
   classicalpot.print_()

   mpi_glob = MPI_context()

   crack_slab, width, height, E, v, v2, bulk = crack_make_slab(params, classicalpot)

   # Save bulk cube (used for qm_rescale_r parameter in crack code)
   if params.qm_args.startswith('TB'):
      bigger_bulk = supercell(bulk, 2, 2, 2)
      bulk = bigger_bulk
   bulk.write(stem+'_bulk.xyz')

   crack_slab.write(stem+'_slab.xyz')

   crack_slab.params['OrigWidth'] = width
   crack_slab.params['OrigHeight'] = height

   crack_slab.params['YoungsModulus'] = E
   crack_slab.params['PoissonRatio_yx'] = v
   crack_slab.params['PoissonRatio_yz'] = v2

   # Open surfaces, remain periodic in z direction (normal to plane)
   # and optionally also in x direction if crack_double_ended is true
   if not params.crack_double_ended:
      crack_slab.lattice[1,1] = crack_slab.lattice[1,1] + params.crack_vacuum_size

   crack_slab.lattice[2,2] = crack_slab.lattice[2,2] + params.crack_vacuum_size
   crack_slab.set_lattice(crack_slab.lattice, False)

   # 3D crack with free surfaces at z = +/- depth/2
   if params.crack_free_surfaces:
      crack_slab.lattice[3,3] = crack_slab.lattice[3,3] + params.crack_vacuum_size

   crack_slab.set_lattice(crack_slab.lattice, False)

   # Map atoms into cell AFTER changing to the new lattice
   crack_slab.map_into_cell()

   miny, maxy = crack_slab.pos[2,:].min(), crack_slab.pos[2,:].max()
   assert abs((maxy-miny) - height) < 1e-5  # be sure that remapping didn't change height of slab

   # Add various properties to crack_slab
   crack_slab.add_property('hybrid', 0)
   crack_slab.add_property('hybrid_mark', HYBRID_NO_MARK)
   crack_slab.add_property('changed_nn', 0)
   crack_slab.add_property('move_mask', 0)
   crack_slab.add_property('nn', 0)
   crack_slab.add_property('old_nn', 0)
   crack_slab.add_property('md_old_changed_nn', 0)
   crack_slab.add_property('edge_mask', 0)
   crack_slab.add_property('crack_surface', False)
   crack_slab.add_property('crack_front', False)

   print_title('Fixing Atoms')

   # Fix top and bottom edges - anything within crack_edge_fix_tol of ymax or ymin is fixed

   miny, maxy = crack_slab.pos[2,:].min(), crack_slab.pos[2,:].max()

   crack_slab.move_mask[:] = 1
   crack_slab.move_mask[(abs(crack_slab.pos[2,:]-maxy) < params.crack_edge_fix_tol) |
                        (abs(crack_slab.pos[2,:]-miny) < params.crack_edge_fix_tol)] = 0

   print('%d atoms. %d fixed atoms' % (crack_slab.n, crack_slab.n - crack_slab.move_mask.count()))

   print_title('Setting edge mask')

   crack_slab.edge_mask[:] = 0

   minx, maxx = crack_slab.pos[1,:].min(), crack_slab.pos[1,:].max()
   crack_slab.edge_mask[(abs(crack_slab.pos[1,:]-minx) < params.selection_edge_tol) |
                        (abs(crack_slab.pos[1,:]-maxx) < params.selection_edge_tol)] = 1

   miny, maxy = crack_slab.pos[2,:].min(), crack_slab.pos[2,:].max()
   crack_slab.edge_mask[(abs(crack_slab.pos[2,:]-miny) < params.selection_edge_tol) |
                        (abs(crack_slab.pos[2,:]-maxy) < params.selection_edge_tol)] = 1

   if params.crack_free_surfaces:
      # Open surfaces at +/- z
      minz, maxz = crack_slab.pos[3,:].min(), crack_slab.pos[3,:].max()
      crack_slab.edge_mask[(abs(crack_slab.pos[3,:]-minz) < params.selection_edge_tol) |
                           (abs(crack_slab.pos[3,:]-maxz) < params.selection_edge_tol)] = 1

   crack_make_seed(crack_slab, params)

   if (params.crack_apply_initial_load):
      crack_calc_load_field(crack_slab, params, classicalpot, params.crack_loading, overwrite_pos=True, mpi=mpi_glob)

   crack_update_connect(crack_slab, params)

   if (not params.simulation_classical):
      if (params.selection_method.strip() == 'crack_front' or
          params.crack_tip_method.strip() == 'local_energy'):
         classicalpot.calc(crack_slab, local_energy=True)

      crack_setup_marks(crack_slab, params)
      crack_update_selection(crack_slab, params)
   
   return crack_slab


def crack_rescale_homogeneous_xy(at, params, new_strain):
    h0 = at.OrigHeight
    h = at.pos[2,:].max() - at.pos[2,:].min()

    eps0 = (h - h0)/h0
    eps1 = new_strain

    print 'Initial strain %.4f' % eps0

    t = np.diag([(1+eps1)/(1+eps0), (1+eps1)/(1+eps0), 1.0])

    b = transform(at, t)

    h1 = b.pos[2,:].max() - b.pos[2,:].min()
    print 'Final strain % .4f' % ((h1 - h0)/h0)

    b.params['G'] = crack_measure_g(b, b.YoungsModulus, b.PoissonRatio_yx, h0)

    crack_update_connect(b, params)

    return b
    
