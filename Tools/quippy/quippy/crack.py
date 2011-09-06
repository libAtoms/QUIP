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

import os, warnings
import numpy as np
from quippy.mpi_context import MPI_context
from quippy.potential import Potential
from quippy.atoms import Atoms
from quippy.structures import transform, supercell
from quippy.system import print_title, verbosity_push
from quippy.clusters import HYBRID_NO_MARK
from quippy.units import GPA
from quippy.surface import J_PER_M2
from quippy.farray import fzeros

__all__ = []
 
try:
   import quippy.cracktools, quippy.crackparams
   from quippy.crackparams import *
   from quippy.cracktools import *
   __all__.extend(quippy.crackparams.__all__)
   __all__.extend(quippy.cracktools.__all__)
   
except ImportError:
   warnings.warn('crack utilities not available')

__all__.extend(['crack_rescale_homogeneous_xy',
                 'makecrack', 'crack_strain_energy_release_rate',
                 'crack_rotation_matrix'])


def makecrack_main(params, stem):
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

    if params.crack_fix_sides:
        maxx, minx = crack_slab.pos[1,:].min(), crack_slab.pos[1,:].max()
        crack_slab.move_mask[(abs(crack_slab.pos[1,:]-maxx) < params.crack_edge_fix_tol) |
                             (abs(crack_slab.pos[1,:]-minx) < params.crack_edge_fix_tol)] = 0


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

    if params.any_per_atom_tau():
        # Set up per_atom_tau property for ramped Langevin thermostat:
        #
        #    tau
        #    ^
        #    |\        /|                     |\        /|  max_tau
        #    | \      / |                     | \      / |
        #    |  \    /  |     constant E      |  \    /  |
        #    |   \  /   |      (tau = 0)      |   \  /   |
        #    |    \/    |                     |    \/    |
        #    +----------+---------------------+----------+---> x
        #   -w/2      -w/2+r                 w/2-r      w/2

        w_by_2   = crack_slab.OrigWidth/2.
        ramp_len = params.crack_thermostat_ramp_length
        max_tau  = params.crack_thermostat_ramp_max_tau
        print 'Adding thermostat ramp with length', ramp_len, 'max_tau', max_tau

        @np.vectorize
        def tau(x):
            if x < -w_by_2 + ramp_len/2:
                q = (x+w_by_2)/(ramp_len/2.)
                return max_tau*(1.- q)
            elif (x > -w_by_2 + ramp_len/2 and
                  x < -w_by_2 + ramp_len):
                q = (x+w_by_2-ramp_len/2.)/(ramp_len/2.)
                return max_tau*q
            elif (x > -w_by_2 + ramp_len and
                  x < w_by_2 - ramp_len):
                return 0.
            elif (x > w_by_2 - ramp_len and
                  x < w_by_2 - ramp_len/2):
                q = (x-w_by_2+ramp_len)/(ramp_len/2.)
                return max_tau*(1.- q)
            else:
                q = (x-w_by_2+ramp_len/2.)/(ramp_len/2.)
                return max_tau*q
        crack_slab.add_property('per_atom_tau', tau(crack_slab.pos[1,:]))

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


def makecrack_initial_velocity_field(params, stem, advance_step, advance_time):
    crack_slab_1 = makecrack_main(params, stem)

    # advance by one bond
    params.crack_seed_length = params.crack_seed_length + advance_step
    crack_slab_2 = makecrack_main(params, stem)

    crack_slab_1.add_property('velo', 0., n_cols=3)
    crack_slab_1.velo[...] = (crack_slab_2.pos - crack_slab_1.pos)/(advance_step*advance_time)

    return crack_slab_1

def makecrack(params, stem):
    if params.crack_initial_velocity_field:
        slab = makecrack_initial_velocity_field(params, stem, params.crack_initial_velocity_field_dx,
                                                params.crack_initial_velocity_field_dt)
    else:
        slab = makecrack_main(params, stem)

    return slab


def crack_strain_energy_release_rate(at, bulk=None, f_min=.8, f_max=.9, stem=None, avg_pos=False):

    print 'Analytical effective elastic modulus E\' = ', at.YoungsModulus/(1-at.PoissonRatio_yx**2), 'GPa'
    print 'Analytical energy release rate G = ', crack_measure_g(at, at.YoungsModulus, at.PoissonRatio_yx, at.OrigHeight), 'J/m^2'

    if bulk is None:
        if stem is None: raise ValueError('Either "bulk" or "stem" must be present')
        bulk = Atoms(stem+'_bulk.xyz')

    if not hasattr(at, 'local_energy') or not hasattr(bulk, 'energy'):
        if stem is None: raise ValueError('local_energy property not found in Atoms and "stem" is missing')
        xmlfile = stem+'.xml'
        params = CrackParams(xmlfile)
        pot = Potential(params.classical_args, param_filename=stem+'.xml')
        pot.print_()

        if not hasattr(at, 'local_energy'):
            if avg_pos:
                tmp_pos = at.pos.copy()
                at.pos[...] = at.avgpos
            at.set_cutoff(pot.cutoff()+1.)
            at.calc_connect()
            pot.calc(at, args_str="local_energy")
            if avg_pos:
                at.pos[...] = tmp_pos

        if not hasattr(bulk, 'energy'):
            bulk.set_cutoff(pot.cutoff()+1.)
            bulk.calc_connect()
            pot.calc(bulk, args_str='energy')

    h = at.pos[2,:].max() - at.pos[2,:].min()
    h0 = at.OrigHeight
    strain = (h - h0)/h0
    print 'Applied strain', strain

    x_min = f_min*at.OrigWidth - at.OrigWidth/2.
    x_max = f_max*at.OrigWidth - at.OrigWidth/2.
    strip = np.logical_and(at.move_mask == 1, np.logical_and(at.pos[1,:] > x_min, at.pos[1,:] < x_max))
    at.add_property('strip', strip, overwrite=True)

    strip_depth = at.lattice[3,3]
    strip_width = at.pos[1,strip].max() - at.pos[1,strip].min()
    strip_height = at.pos[2,strip].max() - at.pos[2,strip].min()
    strip_volume = strip_width*strip_height*strip_depth
    print 'Strip contains', strip.count(), 'atoms', 'width', strip_width, 'height', strip_height, 'volume', strip_volume

    strain_energy_density = (at.local_energy[strip].sum() - bulk.energy/bulk.n*strip.count())/strip_volume

    print 'Strain energy density in strip', strain_energy_density, 'eV/A**3'

    E_effective = 2*strain_energy_density/strain**2*GPA
    print 'Effective elastic modulus E =', E_effective, 'GPa'

    G_effective = strain_energy_density*strip_height*J_PER_M2
    print 'Effective energy release rate G =', G_effective, 'J/m^2'

    return G_effective


def crack_rotation_matrix(unit, y, z=None, x=None, tol=1e-5):
    """Return 3x3 matrix rotation matrix defining a crack with open
    surface defined by the plane `y`=(l,m.n) or (h,k,i,l), and either
    crack tip line `z` or crack propagation direction `x`."""

    axes = fzeros((3,3))
    if len(y) == 4:
        h, k, i, l = y
        y = [h, k, l]

    if (x is None and z is None) or (x is not None and z is not None):
        raise ValueError('exactly one of x and z must be non-null')

    axes[:,2] = np.dot(unit.g.T, y)     # plane defined by y=(lmn)

    if z is not None:
        axes[:,3] = np.dot(unit.lattice, z) # line defined by z=[pqr]

        axes[:,2] = axes[:,2]/axes[:,2].norm()
        axes[:,3] = axes[:,3]/axes[:,3].norm()

        if abs(np.dot(axes[:,2], axes[:,3])) > tol:
            raise ValueError('y (%s) and z (%s) directions are not perpendicular' % (y,z))

        axes[:,1] = np.cross(axes[:,2], axes[:,3])
    else:
        axes[:,1] = np.dot(unit.lattice, x)

        axes[:,2] = axes[:,2]/axes[:,2].norm()
        axes[:,1] = axes[:,1]/axes[:,1].norm()

        if abs(np.dot(axes[:,2], axes[:,3])) > tol:
            raise ValueError('y (%s) and x (%s) directions are not perpendicular' % (y,x))

        axes[:,3] = np.cross(axes[:,1], axes[:,2])

    # Rotation matrix is transpose of axes matrix
    return axes.T

