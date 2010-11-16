#!/usr/bin/env python
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
import sys

def makecrack(params):
   """Given a CrackParams object `param`, construct and return a new crack slab Atoms object."""
   
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
   crack_slab.move_mask[logical_or(abs(crack_slab.pos[2,:]-maxy) < params.crack_edge_fix_tol,
                                   abs(crack_slab.pos[2,:]-miny) < params.crack_edge_fix_tol)] = 0

   print('%d atoms. %d fixed atoms' % (crack_slab.n, crack_slab.n - crack_slab.move_mask.count()))

   print_title('Setting edge mask')

   crack_slab.edge_mask[:] = 0

   minx, maxx = crack_slab.pos[1,:].min(), crack_slab.pos[1,:].max()
   crack_slab.edge_mask[logical_or(abs(crack_slab.pos[1,:]-minx) < params.selection_edge_tol,
                                   abs(crack_slab.pos[1,:]-maxx) < params.selection_edge_tol)] = 1

   miny, maxy = crack_slab.pos[2,:].min(), crack_slab.pos[2,:].max()
   crack_slab.edge_mask[logical_or(abs(crack_slab.pos[2,:]-miny) < params.selection_edge_tol,
                                   abs(crack_slab.pos[2,:]-maxy) < params.selection_edge_tol)] = 1

   if params.crack_free_surfaces:
      # Open surfaces at +/- z
      minz, maxz = crack_slab.pos[3,:].min(), crack_slab.pos[3,:].max()
      crack_slab.edge_mask[logical_or(abs(crack_slab.pos[3,:]-minz) < params.selection_edge_tol,
                                      abs(crack_slab.pos[3,:]-maxz) < params.selection_edge_tol)] = 1

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


         
if __name__ == '__main__':

   try:
      params = CrackParams()

      if len(sys.argv[1:]) < 1:
         print('Usage: makecrack [-V] <stem>')
         print('Reads parameter file <stem>.xml and writes NetCDF output file <stem>.nc')
         print('If -V option is given then XML file is validated against DTD.')
         print('')
         print('Available parameters and their default values are:')
         params.print_()

      validate = False
      if sys.argv[1] == '-V':
         validate = True
         del sys.argv[1]

      stem = sys.argv[1]
      xmlfilename = stem+'.xml'

      print('Reading parameters from file %s with XML validation %s.' %
            (xmlfilename, {True:'enabled', False: 'disabled'}[validate]))

      xmlfile = InOutput(xmlfilename,INPUT)
      params.read_xml(xmlfile, validate=validate)

      crack_slab = makecrack(params)
      if params.io_netcdf:
         crack_slab.write(stem+'.nc')
      else:
         crack_slab.write(stem+'.xyz')

   except RuntimeError, re:
      sys.stderr.write('error: %s\n' % str(re))
      sys.exit(1)
