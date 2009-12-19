#!/usr/bin/env python

from quippy import *
import sys

params = CrackParams()

if len(sys.argv[1:]) != 1:
   print('Usage: makecrack <stem>')
   print('Where <stem>.xml is parameter XML file and <stem>.xyz is the crack slab output XYZ file')
   print('')
   print('Available parameters and their default values are:')
   params.print_()


xmlfilename = stem+'.xml'
ncfilename = stem+'.nc'

print_title('Initialisation')

print('Reading parameters from file '+xmlfilename)

xmlfile = InOutput(xmlfilename,INPUT)
params.read_xml(xmlfile)

verbosity_push(params.io_verbosity)
params.print_()

print("Initialising classical potential with args " + params.classical_args +
      " from file " + xmlfilename)
xmlfile.rewind()
classicalpot = Potential(params.classical_args, xmlfile)
classicalpot.print_()

print('Initialising metapotential')
simple = MetaPotential('Simple', classicalpot, mpi_obj=mpi_glob)
simple.print_()

crack_slab, width, height, e, v, v2, bulk = crack_make_slab(params, classicalpot, simple)

# Save bulk cube (used for qm_rescale_r parameter in crack code)
bulk.write(stem+'_bulk.xyz')

crack_slab.params['OrigWidth'] = width
crack_slab.params['OrigHeight'] = height

crack_slab.params['YoungsModulus'] = E
crack_slab.params['PoissonRatio_yx'] = v
crack_slab.params['PoissonRatio_yz'] = v2

# Open x and y surfaces, remain periodic in z direction (normal to plane)
if params.crack_structure != 'graphene':
   lattice = crack_slab.lattice.copy()
     
   if not params.crack_double_ended:
        lattice[1,1] = lattice[1,1] + params.crack_vacuum_size
   
   lattice[2,2] = lattice[2,2] + params.crack_vacuum_size
   crack_slab.set_lattice(crack_slab, lattice)

# Add various properties to crack_slab
crack_slab.add_property('hybrid', 0)
crack_slab.add_property('hybrid_mark', HYBRID_NO_MARK)
crack_slab.add_property('changed_nn', 0)
crack_slab.add_property('move_mask', 0)
crack_slab.add_property('nn', 0)
crack_slab.add_property('old_nn', 0)
crack_slab.add_property('md_old_changed_nn', 0)
crack_slab.add_property('edge_mask', 0)

if params.crack_structure == 'alpha_quartz':
   crack_slab.add_property('efield', 0.0, n_cols=3)
   crack_slab.add_property('dipoles', 0.0, n_cols=3)
   crack_slab.add_property('efield_old1', 0.0, n_cols=3)
   crack_slab.add_property('efield_old2', 0.0, n_cols=3)
   crack_slab.add_property('efield_old3', 0.0, n_cols=3)

print_title('Fixing Atoms')

# Fix top and bottom edges - anything within crack_edge_fix_tol of ymax or ymin is fixed
maxy = 0.0; miny = 0.0; maxx = 0.0
for i in frange(crack_slab.n):
   if (crack_slab.pos[2,i] > maxy): maxy = crack_slab.pos[2,i]
   if (crack_slab.pos[2,i] < miny): miny = crack_slab.pos[2,i]

crack_slab.move_mask = 1
n_fixed = 0
for i in frange(crack_slab.n):
   if ((abs(crack_slab%pos(2,i)-maxy) < params%crack_edge_fix_tol) or
       (abs(crack_slab%pos(2,i)-miny) < params%crack_edge_fix_tol)):
      crack_slab.move_mask[i] = 0
      n_fixed = n_fixed + 1

print('%d atoms. %d fixed atoms' % (crack_slab.n, n_fixed))
  
crack_make_seed(crack_slab, params)

if (params.crack_apply_initial_load):
   call crack_calc_load_field(crack_slab, params, classicalpot, simple, params.crack_loading, overwrite_pos=True)

print_title('Initialising QM region')

crack_pos = crack_find_crack_pos(crack_slab, params)
crack_slab.params['CrackPosx'] = crack_pos[1]
crack_slab.params['CrackPosy'] = crack_pos[2]

crack_setup_marks(crack_slab, params)

params.io_print_all_properties = True

call crack_print(crack_slab, ncfilename, params)

