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
from math import sqrt, pi
import numpy as np
import numpy.ma as ma
from quippy.mpi_context import MPI_context
from quippy.potential import Potential
from quippy.atoms import Atoms
from quippy.structures import transform, supercell, MillerPlane, MillerDirection
from quippy.system import print_title, verbosity_push
from quippy.clusters import HYBRID_NO_MARK
from quippy.units import GPA
from quippy.surface import J_PER_M2
from quippy.farray import fzeros, farray

MPA_SQRT_M = 1e-3/GPA*sqrt(1.0e10)

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
                'crack_rescale_uniaxial',
                'makecrack',
                'crack_strain_energy_release_rate',
                'crack_strain',
                'crack_find_griffith_load',
                'stress_intensity_factor',
                'make_crack_advance_map',
                'find_crack_tip_coordination',
                'irwin_modeI_crack_tip_stress_field',
                'strain_to_G',
                'G_to_strain',
                'get_strain',
                'get_energy_release_rate',
                'get_stress_intensity_factor',
                'fit_crack_stress_field',
                'find_crack_tip_stress_field',
                'plot_stress_fields',
                'thin_strip_displacement_y',
                'print_crack_system',
                'ConstantStrainRate'])

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
    if params.crack_free_surfaces:
       depth = crack_slab.pos[3,:].max() - crack_slab.pos[3,:].min()
    else:
       depth = crack_slab.lattice[3,3]

    # Save bulk cube (used for qm_rescale_r parameter in crack code)
    if params.qm_args.startswith('TB'):
        bigger_bulk = supercell(bulk, 2, 2, 2)
        bulk = bigger_bulk
    bulk.write(stem+'_bulk.xyz')

    crack_slab.write(stem+'_slab.xyz')

    crack_slab.params['OrigWidth'] = width
    crack_slab.params['OrigHeight'] = height
    crack_slab.params['OrigDepth'] = depth

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
        crack_slab.pos[3,:] -= crack_slab.pos[3,:].mean() # center on z=0
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
    if params.crack_fix_dipoles:
	crack_slab.add_property('fixdip', False)

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


    print('%d atoms. %d fixed atoms' % (crack_slab.n, crack_slab.n - crack_slab.move_mask.sum()))

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

    if params.crack_fix_dipoles:
        print_title('Fixing dipoles')
        crack_slab.fixdip[(abs(crack_slab.pos[2,:]-maxy) < params.crack_fix_dipoles_tol) |
                          (abs(crack_slab.pos[2,:]-miny) < params.crack_fix_dipoles_tol)] = 1
        if params.crack_fix_sides:
                maxx, minx = crack_slab.pos[1,:].min(), crack_slab.pos[1,:].max()
                crack_slab.fixdip[(abs(crack_slab.pos[1,:]-maxx) < params.crack_fix_dipoles_tol) |
                                  (abs(crack_slab.pos[1,:]-minx) < params.crack_fix_dipoles_tol)] = 1


    if params.crack_curved_front:
       crack_make_seed_curved_front(crack_slab, params)
    else:
       crack_make_seed(crack_slab, params)
       if params.crack_apply_initial_load:
          crack_calc_load_field(crack_slab, params, classicalpot,
                              params.crack_loading, overwrite_pos=True,
                              mpi=mpi_glob)

    crack_slab.write('dump.xyz')
    crack_update_connect(crack_slab, params)

    if not params.simulation_classical:
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


try:
   from quippy.cracktools import crack_measure_g as _crack_measure_g
   
   def crack_measure_g(crack, YoungsModulus=None, PoissonRatio_yx=None, OrigHeight=None):
       if YoungsModulus is None:
           YoungsModulus = crack.YoungsModulus
       if PoissonRatio_yx is None:
           PoissonRatio_yx = crack.PoissonRatio_yx
       if OrigHeight is None:
           OrigHeight = crack.OrigHeight
       return _crack_measure_g(crack, YoungsModulus, PoissonRatio_yx, OrigHeight)

   import functools
   crack_measure_g = functools.update_wrapper(crack_measure_g, quippy.cracktools.crack_measure_g)

except ImportError:
   pass

def stress_intensity_factor(at):
    """
    Returns stress instensity factor for mode I loading (K_I) in MPa \sqrt(m)
    """

    G = crack_measure_g(at)
    return crack_g_to_k(G, at.YoungsModulus, at.PoissonRatio_yx, 'plane strain')/1.0e6

def crack_rescale_homogeneous_xy(at, params, strain=None, G=None):
    if strain is None and G is None:
       raise ValueError('either strain or G must be present')

    h0 = at.OrigHeight
    h = at.pos[2,:].max() - at.pos[2,:].min()

    print 'Original height %.4f A' % h0
    print 'Initial measured height %.4f A' % h

    if strain is None:
       strain = crack_g_to_strain(G, at.YoungsModulus, at.PoissonRatio_yx, h0)
    
    eps0 = (h - h0)/h0
    eps1 = strain

    print 'Initial strain_yy %.4f' % eps0
    print 'Initial G %.4f J/m^2' % crack_strain_to_g(eps0, at.YoungsModulus, at.PoissonRatio_yx, h0)

    t = np.diag([(1+eps1)/(1+eps0), (1+eps1)/(1+eps0), 1.0])

    b = transform(at, t)

    h1 = b.pos[2,:].max() - b.pos[2,:].min()
    print 'Final strain_yy % .4f' % ((h1 - h0)/h0)

    b.params['G'] = crack_measure_g(b, b.YoungsModulus, b.PoissonRatio_yx, h0)

    print 'Final G %.4f J/m^2' % b.params['G']

    crack_update_connect(b, params)

    return b

 
def crack_rescale_uniaxial(at, params, strain=None, G=None):
    if strain is None and G is None:
       raise ValueError('either strain or G must be present')

    h0 = at.OrigHeight
    h = at.pos[2,:].max() - at.pos[2,:].min()

    print 'Original height %.4f A' % h0
    print 'Initial measured height %.4f A' % h

    if strain is None:
       strain = crack_g_to_strain(G, at.YoungsModulus, at.PoissonRatio_yx, h0)
    
    eps0 = (h - h0)/h0
    eps1 = strain

    print 'Initial strain %.4f' % eps0
    print 'Initial G %.4f J/m^2' % crack_strain_to_g(eps0, at.YoungsModulus, at.PoissonRatio_yx, h0)

    t = np.diag([1.0, (1+eps1)/(1+eps0), 1.0])

    b = transform(at, t)

    h1 = b.pos[2,:].max() - b.pos[2,:].min()
    print 'Final strain_yy % .4f' % ((h1 - h0)/h0)

    b.params['G'] = crack_measure_g(b, b.YoungsModulus, b.PoissonRatio_yx, h0)

    print 'Final G %.4f J/m^2' % b.params['G']

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
    """
    Compute strain energy release rate G from elastic potential energy in a strip
    """

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
    print 'Strip contains', strip.sum(), 'atoms', 'width', strip_width, 'height', strip_height, 'volume', strip_volume

    strain_energy_density = (at.local_energy[strip].sum() - bulk.energy/bulk.n*strip.sum())/strip_volume

    print 'Strain energy density in strip', strain_energy_density, 'eV/A**3'

    E_effective = 2*strain_energy_density/strain**2*GPA
    print 'Effective elastic modulus E =', E_effective, 'GPa'

    G_effective = strain_energy_density*strip_height*J_PER_M2
    print 'Effective energy release rate G =', G_effective, 'J/m^2'

    return G_effective


def crack_make_seed_curved_front(slab, params):
   """Given a slab, introduce a 3D curved crack front."""

   orig_width = slab.params['OrigWidth']
   orig_height = slab.params['OrigHeight']
   orig_depth  = slab.params['OrigDepth']

   crack_length = params.crack_seed_length/orig_width
   crack_curvature = params.crack_curvature
   crack_depth = 1.0

   if params.crack_g > 0.:
      epsilon_max = crack_g_to_strain(params.crack_g,
                                     slab.YoungsModulus,
                                     slab.PoissonRatio_yx,
                                     orig_height)
   else:
      epsilon_max = params.crack_strain

   # ... then shift so coordinates are in ranges x=[0..a], y=[0..b], z=[0..c]
   slab.pos += np.tile(np.array([orig_width, orig_height, orig_depth])[:,np.newaxis]/2., slab.n)

   xlin = 0.5*crack_length*orig_width
   crack_center = farray([crack_length*orig_width, 0.5*orig_height, 0.])

   epsilon_min = 0.
   epsilon_local = 0.
   dy = 0.

   a = crack_curvature
   x1 = crack_length*orig_width
   x2 = x1 + crack_depth
   z1 = 0.
   z2 = orig_depth
   slope = (x2-x1)/(z2-z1)
   
   dymax = 0.5*epsilon_max*orig_height

   z = slab.pos[3,:]
   crack_x_center = crack_curvature*(z*z - z1*z1) + (slope-crack_curvature*(z2+z1))*(z-z1) + crack_center[1]
   ##slab.add_property('xc', crack_x_center, overwrite=True)

   dymin = abs(slab.pos[2,:]-crack_center[2])*epsilon_max

   dy = dymin
   dy[slab.pos[1,:] < xlin] = dymax
   mask = (slab.pos[1,:] > xlin) & (slab.pos[1,:] <= crack_x_center)
   ##slab.add_property('mask', mask, overwrite=True)
   dy[mask] = (dymax - dymin[mask])/(xlin - crack_x_center[mask])*abs(slab.pos[1,mask] - xlin) + dymax

   dy *= np.sign(slab.pos[2,:] - crack_center[2])
   ##slab.add_property('dy', dy, overwrite=True)
   slab.pos[2,:] += dy

   # Centre slab in cell
   slab.pos += np.tile((np.array([params.crack_vacuum_size/2.]*3) +
                        np.diag(slab.lattice)/2.)[:,np.newaxis], slab.n)
   slab.map_into_cell()

   slab.params['CrackPosx'] = crack_center[1] + params.crack_vacuum_size/2.
   slab.params['CrackPosy'] = crack_center[2] + params.crack_vacuum_size/2.


def crack_strain(at):
   """
   Returns strain of crack specimen
   """
   h = at.pos[2,:].max() - at.pos[2,:].min()
   h0 = at.OrigHeight
   eps = (h - h0)/h0
   return eps

def energy_difference(a,b,pot,eps,relax=False):
   """
   Compute energy difference between configurations `a` and `b` using
   `pot`, after the application of a homogenous-xy rescaling to strain
   `eps`. If `relax`=True, the rescaled confiugrations are
   structurally optimised.
   """
   
   ap = crack_rescale_homogeneous_xy(a, cp, eps)
   bp = crack_rescale_homogeneous_xy(b, cp, eps)
   ap.set_cutoff(pot.cutoff()+1.0)
   ap.calc_connect()
   bp.set_cutoff(pot.cutoff()+1.0)
   bp.calc_connect()
   if relax:
      pot.minim(ap, 'cg', 0.1, 100, do_pos=True, do_lat=False)
      pot.minim(bp, 'cg', 0.1, 100, do_pos=True, do_lat=False)
   pot.calc(ap, args_str="energy")
   pot.calc(bp, args_str="energy")
   print "ENERGY_DIFFERENCE %f %f %f" % (crack_strain(ap), crack_measure_g(ap), ap.energy-bp.energy)
   return ap.energy-bp.energy


def crack_find_griffith_load(a, b, pot, relax=False):
   """
   Given two configurations (a, b) which differ by one broken bond,
   find the Griffith load, that is the load at which `a` and `b` have
   the same energy accorinding to the model potential `pot`.

   Returns (strain, G, a_rescaled, b_rescaled).
   """

   eps0 = crack_strain(a)

   def func(x):
       eps, = x
       return energy_difference(a,b,p,eps,relax=relax)**2

   ds = DownhillSimplex(func, x0=[eps0], deltas=[0.001], ftol=0.01)
   eps, = ds.minimise()
   ap = crack_rescale_homogeneous_xy(a, cp, eps)
   bp = crack_rescale_homogeneous_xy(b, cp, eps)

   print 'Griffith critical strain = %f' % crack_strain(ap)
   print 'Griffith critical G = %f J/m^2' % crack_measure_g(ap)

   return (crack_strain(ap), crack_measure_g(ap), ap, bp)


def find_mapping(at_in, vectors, lattice_constant, tol=1e-4):
   """
   Find a mapping between pairs of atoms displaced by a given vector
   (or by one of a number of given vectors).
   """

   class FoundMapping:
      pass

   mapping = fzeros(at_in.n, dtype=np.int32)
   at = at_in.copy()

   # cutoff should be larger than all displacment vectors
   vlen = [farray(vector).norm()*lattice_constant for vector in vectors]
   at.set_cutoff(max(vlen)+0.5)
   print 'cutoff = ', at.cutoff
   at.calc_connect()
   
   for i in at.indices:
      try:
         for neighb in at.neighbours[i]:
            for vector in vectors:
               if ((neighb.diff/lattice_constant - vector)**2).sum() < tol:
                  print i, neighb.j, vector
                  mapping[i] = neighb.j
                  raise FoundMapping

         # corresponding atom is off the edge, so map atom onto itself
         print i, 'not found!'
         mapping[i] = i
      except FoundMapping:
         continue

   return mapping
            

crack_advance_table = {
   'Si(110)[001b]+1':  ([[-1.0/(2*np.sqrt(2)), 0.0,  0.25],
                        [-1.0/(2*np.sqrt(2)), 0.0, -0.25]],
                        [[1.0, 0.0,  0.0],
                         [0.0, 1.0,  0.0],
                         [0.0, 0.0, -1.0]]),
   'Si(110)[001b]-1':  ([[1.0/(2*np.sqrt(2)), 0.0,  0.25],
                        [1.0/(2*np.sqrt(2)), 0.0, -0.25]],
                        [[1.0, 0.0,  0.0],
                         [0.0, 1.0,  0.0],
                         [0.0, 0.0, -1.0]]),
   'Si(110)[001b]+2':  ([[-1.0/np.sqrt(2), 0.0,  0.0]],
                        [[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]]),
   'Si(110)[001b]-2':  ([[1.0/np.sqrt(2), 0.0,  0.0]],
                        [[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]])
}


def crack_unit_advance(slab, at, lattice_constant, vectors=None, transformation=None, system=None):
   """
   Given isoatomic slab and crack configurations, return copy with crack advanced by one bond
   """

   if system is not None:
      try:
         vectors, transformation = crack_advance_table[system]
      except KeyError:
         raise ValueError("unknown crack system %s" % system)

   if vectors is None:
      raise ValueError("missing vectors")

   # find mapping from old to new atom indices
   mapping = find_mapping(slab, vectors, lattice_constant)

   disp = at.pos - slab.pos
   new_disp = disp[mapping]

   # optionally apply a transformation to the new displacements
   if transformation is not None:
      new_disp = np.dot(transformation, new_disp)

   shift = [0.0, 0.0, lattice_constant/2.0]

   result = at.copy()
   result.pos[...] = slab.pos + new_disp
   return result


# displacement vectors from an atom to the atom horizontally ahead
# of it (along x), for a few different crack systems, in scaled
# coordinates. There are two possibilities in each case, one in
# front and one behind (along z axis)

crack_advance_step_frac_coords = {
    '(111)[01-1]':
        np.array([[sqrt(6.) / 4, 0., +sqrt(2.) / 4],
                  [sqrt(6.) / 4, 0., -sqrt(2.) / 4.]]),

    '(110)[001]':
        np.array([[sqrt(2.) / 4., 0., 0.25],
                  [sqrt(2.) / 4., 0., -0.25]]),

    '(110)[1-10]':
        np.array([[0.25, 0., +sqrt(2.) / 4],
                  [0.75, 0., -sqrt(2.) / 4]])
    }


def make_crack_advance_map(atoms, tol=1e-3):
    """
    Find mapping from atom indices to the index of atom one step ahead
    of them in the crack propagation direction (i.e. along +x).

    Requires 'LatticeConstant', 'CleavagePlane', and 'CrackFront' to
    be available in atoms.info dictionary.

    Returns integer array of shape (len(atoms),), and also adds a new
    array 'advance_map' into the Atoms object.
    """

    def find_atom_with_displacement(atoms, i, target_disps):
        """
        Return index of atom with relative displacement from i in target_disps

        Requires connnectivity to be calculated already
        """
        indices, offsets = atoms.neighbours.get_neighbors(i)
        diffs = (atoms.positions[indices] + np.dot(offsets, atoms.cell) -
                 atoms.positions[i])
        for j, diff in zip(indices, diffs):
            for target in target_disps:
                if all(abs(diff - target) < tol):
                    return j
        return 0

    a0 = atoms.info['LatticeConstant']
    cleavage_plane = atoms.info['CleavagePlane']
    crack_front = atoms.info['CrackFront']

    # lookup this crack system in the dictionary of crack advance steps
    cleavage_plane = MillerPlane(cleavage_plane)
    crack_front = MillerDirection(crack_front)
    key = str(cleavage_plane) + str(crack_front)
    steps = a0 * crack_advance_step_frac_coords[key]
    
    advance_map = np.zeros(len(atoms), dtype=int)

    # find biggest distance between atoms before and after a crack advance step
    max_step_length = np.sqrt((steps ** 2).sum(axis=1)).max()

    # convert from ase.Atoms to quippy.Atoms, so we can use faster
    # neighbour lists in Fortran code
    tmp_atoms = Atoms(atoms)
    tmp_atoms.set_cutoff(max_step_length + .1)
    tmp_atoms.calc_connect()

    for i in range(len(tmp_atoms)):
        advance_map[i] = find_atom_with_displacement(tmp_atoms, i, steps)
        
    # save the map inside the Atoms object, and return a copy
    atoms.new_array('advance_map', advance_map)
    return advance_map


def find_crack_tip_coordination(atoms, edge_tol=10.0,
                                strip_height=30.0, nneightol=1.3):
    """
    Return position of crack tip in `atoms`, based on atomic coordination.

    If `atoms` does not contain an `advance_map` property, then
    :func:`make_crack_advance_map` is called to generate the map.

    Parameters
    ----------
    atoms : :class:`~.Atoms' object
       The Atoms object containing the crack slab.
    edge_tol : float
       Distance from edge of system within which to exclude
       undercoodinated atoms.
    strip_height : float
       Height of strip along centre of slab in which to look
       for the track.
    nneightol : float
       Nearest neighbour tolerance, as a fraction of sum of
       covalent radii of atomic species.

    Returns
    -------
    crack_pos : array
       x, y, and z coordinates of the crack tip. Also set in ``CrackPos``
       in ``atoms.info`` dictionary.
    tip_atoms : array
       Indices of atoms near the tip Also set in ``crack_tip`` property.
    """

    old_tip_pos_y = 0
    if 'CrackPos' in atoms.info:
        old_tip_pos_y = atoms.info['CrackPos'][1]

    # Make a copy of atoms as a quippy.Atoms instance, overwriting
    # positions with time-averages values if they are available, and
    # then calculate connectivity using nneightol
    tmp_atoms = Atoms(atoms)
    if 'avgpos' in tmp_atoms.arrays:
        tmp_atoms.set_positions(tmp_atoms.arrays['avgpos'])
    tmp_atoms.set_cutoff_factor(nneightol)
    tmp_atoms.calc_connect()

    nn = tmp_atoms.n_neighb
    x = tmp_atoms.positions[:, 0]
    y = tmp_atoms.positions[:, 1]

    # find undercoordinated atoms in a central strip, and not too
    # close to the left or right edges
    left = tmp_atoms.positions[:, 0].min()
    right = tmp_atoms.positions[:, 0].max()
    uc = ((nn < 4) &
          (abs(y) < strip_height) &
          (x > left + edge_tol) &
          (x < right - edge_tol))

    # position of furthest forward undercoordinated atom ABOVE old tip position
    x_above = x[uc & (y > old_tip_pos_y)].max()

    # position of furthest forward undercoordinated atom BELOW old tip position
    x_below = x[uc & (y < old_tip_pos_y)].max()

    # rightmost undercoordinated atoms, both above and below old tip
    rightmost_uc = uc & (((y > old_tip_pos_y) & (x_above == x)) |
                         ((y <= old_tip_pos_y) & (x_below == x)))

    # we want the NEXT pair of atoms, so we use the saved mapping from
    # atom indices to the indices of atoms one unit cell to the right
    if 'advance_map' not in atoms.arrays:
        print('Generating crack advance map...')
        make_crack_advance_map(atoms)
        
    advance_map = atoms.arrays['advance_map']
    tip_atoms = advance_map[rightmost_uc]
    tip_pos = tmp_atoms.positions[tip_atoms, :].mean(axis=0)

    # Also save results in Atoms (useful for visualisation)
    atoms.info['CrackPos'] = tip_pos
    atoms.set_array('crack_tip', np.array([False]*len(atoms)))
    crack_tip = atoms.arrays['crack_tip']
    crack_tip[tip_atoms] = True
    return tip_pos


def irwin_modeI_crack_tip_stress_field(K, r, t, xy_only=True,
                                       nu=0.5, stress_state='plane strain'):
    """
    Compute Irwin singular crack tip stress field

    Parameters
    ----------
    K : float
       Mode I stress intensity factor. Units should match units of `r`.
    r : array_like
       Radial distances from crack tip. Can be a multidimensional
       array to evaluate stress field on a grid.
    t : array_like
       Angles from horzontal line y=0 ahead of crack tip,
       measured anticlockwise. Should have same shape as `r`.
    xy_only : bool
       If True (default) only xx, yy, xy and yx components will be set.
    nu : float
       Poisson ratio. Used only when ``xy_only=False``, to determine zz stresses
    stress_state : str
       One of"plane stress" or "plane strain". Used if xyz_only=False to
       determine zz stresses.
       
    Returns
    -------
    sigma : array with shape ``r.shape + (3,3)``
    """

    if r.shape != t.shape:
        raise ValueError('shapes of radial and angular arrays "r" and "t" must match')
    
    if stress_state not in ['plane strain', 'plane stress']:
        raise ValueError('stress_state should be either "plane strain" or "plane stress".')

    sigma = np.zeros(r.shape + (3, 3))
    radial = K*1./np.sqrt(2*pi*r)

    sigma[...,0,0] = radial*np.cos(t/2.0)*(1.0 - np.sin(t/2.0)*np.sin(3.0*t/2.0)) # xx
    sigma[...,1,1] = radial*np.cos(t/2.0)*(1.0 + np.sin(t/2.0)*np.sin(3.0*t/2.0)) # yy
    sigma[...,0,1] = radial*np.sin(t/2.0)*np.cos(t/2.0)*np.cos(3.0*t/2.0)         # xy
    sigma[...,1,0] = sigma[...,0,1]                                               # yx=xy

    if not xy_only and stress_state == 'plane strain':
        sigma[...,2,2] = nu*(sigma[...,0,0] + sigma[...,1,1])              # zz

    return sigma


class IrwinStressField(object):
    """
    Calculator to return Irwin near-tip stress field at atomic sites
    """
    def __init__(self, K=None, x0=None, y0=None, sxx0=0.0, syy0=0., sxy0=0., nu=0.5,
                 stress_state='plane strain'):
        self.K = K
        self.x0 = x0
        self.y0 = y0
        self.sxx0 = sxx0
        self.syy0 = syy0
        self.sxy0 = sxy0
        self.nu = nu
        self.stress_state = stress_state

    def get_stresses(self, atoms):
        K = self.K
        if K is None:
            K = get_stress_intensity_factor(atoms)

        x0, y0 = self.x0, self.y0
        if x0 is None:
            x0 = atoms.info['CrackPos'][0]
        if y0 is None:
            y0 = atoms.info['CrackPos'][1]
            
        x = atoms.positions[:, 0]
        y = atoms.positions[:, 1]
        r = np.sqrt((x - x0)**2 + (y - y0)**2)
        t = np.arctan2(y - y0, x - x0)
        
        sigma = irwin_modeI_crack_tip_stress_field(K, r, t, self.nu,
                                                   self.stress_state)
        sigma[:,0,0] += self.sxx0
        sigma[:,1,1] += self.syy0
        sigma[:,0,1] += self.sxy0
        sigma[:,1,0] += self.sxy0

        return sigma


def strain_to_G(strain, E, nu, orig_height):
    """
    Convert from strain to energy release rate G for thin strip geometry

    Parameters
    ----------
    strain : float
       Dimensionless ratio ``(current_height - orig_height)/orig_height``
    E : float
       Young's modulus relevant for a pull in y direction sigma_yy/eps_yy
    nu : float
       Poission ratio -eps_yy/eps_xx
    orig_height : float
       Unstrained height of slab

    Returns
    -------
    G : float
       Energy release rate in units consistent with input
       (i.e. in eV/A**2 if eV/A/fs units used)
    """
    return 0.5 * E / (1.0 - nu * nu) * strain * strain * orig_height


def G_to_strain(G, E, nu, orig_height):
    """
    Convert from energy release rate G to strain for thin strip geometry

    Parameters
    ----------
    G : float
       Energy release rate in units consistent with `E` and `orig_height`
    E : float
       Young's modulus relevant for a pull in y direction sigma_yy/eps_yy
    nu : float
       Poission ratio -eps_yy/eps_xx
    orig_height : float
       Unstrained height of slab

    Returns
    -------
    strain : float
       Dimensionless ratio ``(current_height - orig_height)/orig_height``
    """
    return sqrt(2.0 * G * (1.0 - nu * nu) / (E * orig_height))


def get_strain(atoms):
    """
    Return the current strain on thin strip configuration `atoms`

    Requires unstrained height of slab to be stored as ``OrigHeight``
    key in ``atoms.info`` dictionary.

    Also updates value stored in ``atoms.info``.
    """
    
    orig_height = atoms.info['OrigHeight']
    current_height = atoms.positions[:, 1].max() - atoms.positions[:, 1].min()
    strain = current_height / orig_height - 1.0
    atoms.info['strain'] = strain
    return strain


def get_energy_release_rate(atoms):
    """
    Return the current energy release rate G for `atoms

    Also updates value stored in ``atoms.info`` dictionary.
    """
    
    current_strain = get_strain(atoms)
    orig_height = atoms.info['OrigHeight']
    E = atoms.info['YoungsModulus']
    nu = atoms.info['PoissonRatio_yx']
    G = strain_to_G(current_strain, E, nu, orig_height)
    atoms.info['G'] = G
    return G


def get_stress_intensity_factor(atoms):
    """
    Return stress intensity factor K_I

    Also updates value stored in ``atoms.info`` dictionary.
    """

    G = get_energy_release_rate(atoms)
    
    E = atoms.info['YoungsModulus']
    nu = atoms.info['PoissonRatio_yx']

    Ep = E/(1-nu**2)
    K = sqrt(G/Ep)
    atoms.info['K'] = K
    return K


def fit_crack_stress_field(atoms, r_range=(0., 50.), initial_params=None, fix_params=None,
                           sigma=None, avg_sigma=None, avg_decay=0.005, calc=None, verbose=False):
    """
    Perform a least squares fit of near-tip stress field to Irwin solution

    Stresses on the atoms are fit to the Irwin K-field singular crack tip
    solution, allowingthe crack position, stress intensity factor and
    far-field stress components to vary during the fit.

    Parameters
    ----------
    atoms : :class:`~.Atoms` object
       Crack system. For the initial fit, the following keys are used
       from the :attr:`~Atoms.info` dictionary:

          - ``YoungsModulus``
          - ``PossionRatio_yx``
          - ``G`` --- current energy release rate
          - ``strain`` --- current applied strain
          - ``CrackPos`` --- initial guess for crack tip position

       The initial guesses for the stress intensity factor ``K`` are
       far-field stress ``sigma0`` are computed from
       ``YoungsModulus``, ``PoissonRatio_yx``, ``G`` and ``strain``,
       assuming plane strain in thin strip boundary conditions.

       On exit, new ``K``, ``sigma0`` and ``CrackPos`` entries are set
       in the :attr:`~Atoms.info` dictionary. These values are then
       used as starting guesses for subsequent fits.

    r_range : sequence of two floats, optional
       If present, restrict the stress fit to an annular region
       ``r_range[0] <= r < r_range[1]``, centred on the previous crack
       position (from the ``CrackPos`` entry in ``atoms.info``). If
       r_range is ``None``, fit is carried out for all atoms.

    initial_params : dict
       Names and initial values of parameters. Missing initial values
       are guessed from Atoms object.

    fix_params : dict
       Names and values of parameters to fix during the fit,
       e.g. ``{y0: 0.0}`` to constrain the fit to the line y=0

    sigma : None or array with shape (len(atoms), 3, 3)
       Explicitly provide the per-atom stresses. Avoids calling Atoms'
       calculators :meth:`~.get_stresses` method.

    avg_sigma : None or array with shape (len(atoms), 3, 3)
       If present, use this array to accumulate the time-averaged
       stress field. Useful when processing a trajectory.

    avg_decay : real
       Factor by which average stress is attenuated at each step.
       Should be set to ``dt/tau`` where ``dt`` is MD time-step
       and ``tau`` is a characteristic averaging time.

    calc : Calculator object, optional
       If present, override the calculator used to compute stresses
       on the atoms. Default is ``atoms.get_calculator``.

       To use the atom resolved stress tensor pass an instance of the 
       :class:`~quippy.elasticity.AtomResolvedStressField` class.

    verbose : bool, optional
       If set to True, print additional information about the fit.

    Returns
    -------
    params : dict with keys ``[K, x0, y0, sxx0, syy0, sxy0]``
       Fitted parameters, in a form suitable for passin
       :class:`IrwinStressField` constructor. These are the stress intensity
       factor `K`, the centre of the stress field ``(x0, y0)``, and the
       far field contribution to the stress ``(sxx0, syy0, sxy0)``.
    """

    params = {}
    if initial_params is not None:
       params.update(initial_params)

    if 'K' not in params:
       # Guess for stress intensity factor K
       if 'K' in atoms.info:
           params['K'] = atoms.info['K']
       else:
           try:
               params['K'] = get_stress_intensity_factor(atoms)
           except KeyError:
               params['K'] = 1.0*MPA_SQRT_M

    if 'sxx0' not in params or 'syy0' not in params or 'sxy0' not in params:
       # Guess for far-field stress
       if 'sigma0' in atoms.info:
          params['sxx0'], params['syy0'], params['sxy0'] = atoms.info['sigma0']
       else:
          try:
              E = atoms.info['YoungsModulus']
              nu = atoms.info['PoissonRatio_yx']
              Ep = E/(1-nu**2)
              params['syy0'] = Ep*atoms.info['strain']
              params['sxx0'] = nu*params['syy0']
              params['sxy0'] = 0.0
          except KeyError:
              params['syy0'] = 0.0
              params['sxx0'] = 0.0
              params['sxy0'] = 0.0

    if 'x0' not in params or 'y0' not in params:
       # Guess for crack position
       try:
           params['x0'], params['y0'], _ = atoms.info['CrackPos']
       except KeyError:
           params['x0'] = (atoms.positions[:, 0].min() +
                           (atoms.positions[:, 0].max() - atoms.positions[:, 0].min())/3.0)
           params['y0'] = 0.0

    # Override any fixed parameters
    if fix_params is None:
       fix_params = {}
    params.update(fix_params)

    x = atoms.positions[:, 0]
    y = atoms.positions[:, 1]
    r = np.sqrt((x - params['x0'])**2 + (y - params['y0'])**2)

    # Get local stresses
    if sigma is None:
       if calc is None:
           calc = atoms.get_calculator()
       sigma = calc.get_stresses(atoms)

    # Update avg_sigma in place
    if avg_sigma is not None:
       avg_sigma[...] = np.exp(-avg_decay)*avg_sigma + (1.0 - np.exp(-avg_decay))*sigma
       sigma = avg_sigma.copy()

    # Zero components out of the xy plane
    sigma[:,2,2] = 0.0
    sigma[:,0,2] = 0.0
    sigma[:,2,0] = 0.0
    sigma[:,1,2] = 0.0
    sigma[:,2.1] = 0.0

    mask = Ellipsis # all atoms
    if r_range is not None:
        rmin, rmax = r_range
        mask = (r > rmin) & (r < rmax)

    if verbose:
       print 'Fitting on %r atoms' % sigma[mask,1,1].shape
    
    def objective_function(params, x, y, sigma, var_params):
        params = dict(zip(var_params, params))
        if fix_params is not None:
            params.update(fix_params)
        irwin_sigma = IrwinStressField(**params).get_stresses(atoms)
        delta_sigma = sigma[mask,:,:] - irwin_sigma[mask,:,:]
        return delta_sigma.reshape(delta_sigma.size)

    # names and values of parameters which can vary in this fit
    var_params = sorted([key for key in params.keys() if key not in fix_params.keys() ])
    initial_params = [params[key] for key in var_params]

    from scipy.optimize import leastsq
    fitted_params, cov, infodict, mesg, success = leastsq(objective_function,
                                                         initial_params,
                                                         args=(x, y, sigma, var_params),
                                                         full_output=True)

    params = dict(zip(var_params, fitted_params))
    params.update(fix_params)

    # estimate variance in parameter estimates
    if cov is None:
       # singular covariance matrix
       err = dict(zip(var_params, [0.]*len(fitted_params)))
    else:
       s_sq = (objective_function(fitted_params, x, y, sigma, var_params)**2).sum()/(sigma.size-len(fitted_params))
       cov = cov * s_sq
       err = dict(zip(var_params, np.sqrt(np.diag(cov))))
    
    if verbose:
       print 'K = %.3f MPa sqrt(m)' % (params['K']/MPA_SQRT_M)
       print 'sigma^0_{xx,yy,xy} = (%.1f, %.1f, %.1f) GPa' % (params['sxx0']*GPA,
                                                              params['syy0']*GPA,
                                                              params['sxy0']*GPA)
       print 'Crack position (x0, y0) = (%.1f, %.1f) A' % (params['x0'], params['y0'])

    atoms.info['K'] = params['K']
    atoms.info['sigma0'] = (params['sxx0'], params['syy0'], params['sxy0'])
    atoms.info['CrackPos'] = np.array((params['x0'], params['y0'], atoms.cell[2,2]/2.0))

    return params, err



def find_crack_tip_stress_field(atoms, r_range=None, initial_params=None, fix_params=None,
                                sigma=None, avg_sigma=None, avg_decay=0.005, calc=None):
    """
    Find the position of the crack tip by fitting to the Irwin `K`-field solution

    Fit is carried out using :func:`fit_crack_stress_field`, and parameters
    have the same meaning as there.

    See also
    --------
    fit_crack_stress_field
    """

    params, err = fit_crack_stress_field(atoms, r_range, initial_params, fix_params, sigma,
                                         avg_sigma, avg_decay, calc)

    return np.array((params['x0'], params['y0'], atoms.cell[2,2]/2.0))
    

def plot_stress_fields(atoms, r_range=None, initial_params=None, fix_params=None,
                       sigma=None, avg_sigma=None, avg_decay=0.005, calc=None):
    """
    Fit and plot atomistic and continuum stress fields

    Firstly a fit to the Irwin `K`-field solution is carried out using
    :func:`fit_crack_stress_field`, and parameters have the same
    meaning as for that function. Then plots of the
    :math:`\sigma_{xx}`, :math:`\sigma_{yy}`, :math:`\sigma_{xy}`
    fields are produced for atomistic and continuum cases, and for the
    residual error after fitting.
    """

    from pylab import griddata, meshgrid, subplot, cla, contourf, colorbar, draw, title, clf, gca

    params, err = fit_crack_stress_field(atoms, r_range, initial_params, fix_params, sigma,
                                         avg_sigma, avg_decay, calc)

    K, x0, y0, sxx0, syy0, sxy0 = (params['K'], params['x0'], params['y0'],
                                   params['sxx0'], params['syy0'], params['sxy0'])
   
    x = atoms.positions[:, 0]
    y = atoms.positions[:, 1]

    X = np.linspace((x-x0).min(), (x-x0).max(), 500)
    Y = np.linspace((y-y0).min(), (y-y0).max(), 500)

    t = np.arctan2(y-y0, x-x0)
    r = np.sqrt((x-x0)**2 + (y-y0)**2)

    if r_range is not None:
       rmin, rmax = r_range
       mask = (r > rmin) & (r < rmax)
    else:
       mask = Ellipsis

    atom_sigma = sigma
    if atom_sigma is None:
       atom_sigma = atoms.get_stresses()

    grid_sigma = np.dstack([griddata(x[mask]-x0, y[mask]-y0, atom_sigma[mask,0,0], X, Y),
                            griddata(x[mask]-x0, y[mask]-y0, atom_sigma[mask,1,1], X, Y),
                            griddata(x[mask]-x0, y[mask]-y0, atom_sigma[mask,0,1], X, Y)])

    X, Y = meshgrid(X, Y)
    R = np.sqrt(X**2+Y**2)
    T = np.arctan2(Y, X)

    grid_sigma[((R < rmin) | (R > rmax)),:] = np.nan # mask outside fitting region

    irwin_sigma = irwin_modeI_crack_tip_stress_field(K, R, T, x0, y0)
    irwin_sigma[...,0,0] += sxx0
    irwin_sigma[...,1,1] += syy0
    irwin_sigma[...,0,1] += sxy0
    irwin_sigma[...,1,0] += sxy0
    irwin_sigma = ma.masked_array(irwin_sigma, mask=grid_sigma.mask)

    irwin_sigma[((R < rmin) | (R > rmax)),:,:] = np.nan # mask outside fitting region

    contours = [np.linspace(0, 20, 10),
                np.linspace(0, 20, 10),
                np.linspace(-10,10, 10)]

    dcontours = [np.linspace(0, 5, 10),
                np.linspace(0, 5, 10),
                np.linspace(-5, 5, 10)]

    clf()
    for i, (ii, jj), label in zip(range(3),
                                  [(0,0), (1,1), (0,1)],
                                  ['\sigma_{xx}', r'\sigma_{yy}', r'\sigma_{xy}']):
        subplot(3,3,i+1)
        gca().set_aspect('equal')
        contourf(X, Y, grid_sigma[...,i]*GPA, contours[i])
        colorbar()
        title(r'$%s^\mathrm{atom}$' % label)
        draw()

        subplot(3,3,i+4)
        gca().set_aspect('equal')
        contourf(X, Y, irwin_sigma[...,ii,jj]*GPA, contours[i])
        colorbar()
        title(r'$%s^\mathrm{Irwin}$' % label)
        draw()

        subplot(3,3,i+7)
        gca().set_aspect('equal')
        contourf(X, Y, abs(grid_sigma[...,i] -
                           irwin_sigma[...,ii,jj])*GPA, dcontours[i])
        colorbar()
        title(r'$|%s^\mathrm{atom} - %s^\mathrm{Irwin}|$' % (label, label))
        draw()



def thin_strip_displacement_y(x, y, strain, a, b):
    """
    Return vertical displacement ramp used to apply initial strain to slab

    Strain is increased from 0 to strain over distance :math:`a <= x <= b`.
    Region :math:`x < a` is rigidly shifted up/down by ``strain*height/2``.

    Here is an example of how to use this function on an artificial
    2D square atomic lattice. The positions are plotted before (left)
    and after (right) applying the displacement, and the horizontal and
    vertical lines show the `strain` (red), `a` (green) and `b` (blue)
    parameters. ::

       import matplotlib.pyplot as plt
       import numpy as np

       w = 1; h = 1; strain = 0.1; a = -0.5; b =  0.0
       x = np.linspace(-w, w, 20)
       y = np.linspace(-h, h, 20)
       X, Y = np.meshgrid(x, y)
       u_y = thin_strip_displacement_y(X, Y, strain, a, b)

       for i, disp in enumerate([0, u_y]):
           plt.subplot(1,2,i+1)
           plt.scatter(X, Y + disp, c='k', s=5)
           for y in [-h, h]:
               plt.axhline(y, color='r', linewidth=2, linestyle='dashed')
               plt.axhline(y*(1+strain), color='r', linewidth=2)
           for x, c in zip([a, b], ['g', 'b']):
               plt.axvline(x, color=c, linewidth=2)

    .. image:: thin-strip-displacement-y.png
       :width: 600
       :align: center

    Parameters
    ----------
    x : array
    y : array
       Atomic positions in unstrained slab, centered on origin x=0,y=0
    strain : float
       Far field strain to apply
    a : float
       x coordinate for beginning of strain ramp
    b : float
       x coordinate for end of strain ramp
    """

    u_y = np.zeros_like(y)
    height = y.max() - y.min()     # measure height of slab
    shift = strain * height / 2.0  # far behind crack, shift = strain*height/2

    u_y[x < a] = np.sign(y[x < a]) * shift  # region shift for x < a
    u_y[x > b] = strain * y[x > b]          # constant strain for x > b
    
    middle = (x >= a) & (x <= b)            # interpolate for a <= x <= b
    f = (x[middle] - a) / (b - a)
    u_y[middle] = (f * strain * y[middle] +
                   (1 - f) * shift * np.sign(y[middle]))
    
    return u_y


def print_crack_system(crack_direction, cleavage_plane, crack_front):
    """
    Pretty printing of crack crystallographic coordinate system 

    Specified by Miller indices for crack_direction (x),
    cleavage_plane (y) and crack_front (z), each of which should be
    a sequence of three floats
    """
    crack_direction = MillerDirection(crack_direction)
    cleavage_plane = MillerPlane(cleavage_plane)
    crack_front = MillerDirection(crack_front)

    print('Crack system              %s%s' % (cleavage_plane, crack_front))
    print('Crack direction (x-axis)  %s' % crack_direction)
    print('Cleavage plane  (y-axis)  %s' % cleavage_plane)
    print('Crack front     (z-axis)  %s\n' % crack_front)



class ConstantStrainRate(object):
    """
    Constraint which increments epsilon_yy at a constant strain rate

    Rescaling is applied only to atoms where `mask` is True (default is all atoms)
    """

    def __init__(self, orig_height, delta_strain, mask=None):
        self.orig_height = orig_height
        self.delta_strain = delta_strain
        if mask is None:
            mask = Ellipsis
        self.mask = mask

    def adjust_forces(self, positions, forces):
        pass

    def adjust_positions(self, oldpos, newpos):
        current_height = newpos[:, 1].max() - newpos[:, 1].min()
        current_strain = current_height / self.orig_height - 1.0
        new_strain = current_strain + self.delta_strain
        alpha = (1.0 + new_strain) / (1.0 + current_strain)
        newpos[self.mask, 1] = newpos[self.mask, 1]*alpha



