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
import sys

def makecrack(params):
   """Given a CrackParams object `param`, construct and return a new crack slab Atoms object."""
   
   print_title('Initialisation')

   verbosity_push(params.io_verbosity)
   params.print_()

   print("Initialising classical potential with args " + params.classical_args.strip() +
         " from file " + xmlfilename)
   xmlfile.rewind()
   classicalpot = Potential(params.classical_args, xmlfile)
   classicalpot.print_()

   mpi_glob = MPI_context()

   print('Initialising metapotential')
   simple = MetaPotential('Simple', classicalpot, mpi_obj=mpi_glob)
   simple.print_()

   crack_slab, width, height, E, v, v2, bulk = crack_make_slab(params, classicalpot, simple)

   # Save bulk cube (used for qm_rescale_r parameter in crack code)
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

   # Add various properties to crack_slab
   crack_slab.add_property('hybrid', 0)
   crack_slab.add_property('hybrid_mark', HYBRID_NO_MARK)
   crack_slab.add_property('changed_nn', 0)
   crack_slab.add_property('move_mask', 0)
   crack_slab.add_property('nn', 0)
   crack_slab.add_property('old_nn', 0)
   crack_slab.add_property('md_old_changed_nn', 0)
   crack_slab.add_property('edge_mask', 0)

   print_title('Fixing Atoms')

   # Fix top and bottom edges - anything within crack_edge_fix_tol of ymax or ymin is fixed

   miny, maxy = crack_slab.pos[2,:].min(), crack_slab.pos[2,:].max()

   crack_slab.move_mask[:] = 1
   crack_slab.move_mask[logical_or(abs(crack_slab.pos[2,:]-maxy) < params.crack_edge_fix_tol,
                                   abs(crack_slab.pos[2,:]-miny) < params.crack_edge_fix_tol)] = 0

   print('%d atoms. %d fixed atoms' % (crack_slab.n, crack_slab.n - crack_slab.move_mask.count()))

   crack_make_seed(crack_slab, params)
   crack_setup_marks(crack_slab, params)

   if (params.crack_apply_initial_load):
      crack_calc_load_field(crack_slab, params, classicalpot, simple, params.crack_loading, overwrite_pos=True, mpi=mpi_glob)

   if (not params.simulation_classical):
      crack_update_selection(crack_slab, params)
   
   return crack_slab


if __name__ == '__main__':

   params = CrackParams()

   if len(sys.argv[1:]) != 1:
      print('Usage: makecrack <stem>')
      print('Reads parameter file <stem>.xml  and writes NetCDF output file <stem>.nc')
      print('')
      print('Available parameters and their default values are:')
      params.print_()

   stem = sys.argv[1]
   xmlfilename = stem+'.xml'
   ncfilename = stem+'.nc'

   print('Reading parameters from file '+xmlfilename)

   xmlfile = InOutput(xmlfilename,INPUT)
   params.read_xml(xmlfile)

   crack_slab = makecrack(params)
   crack_slab.write(ncfilename)



