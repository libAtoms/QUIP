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
from quippy.farray import fzeros, farray

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
                'stress_intensity_factor'])

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


