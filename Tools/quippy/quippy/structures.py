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

import numpy as np
from quippy.atoms import *
import quippy._structures
from quippy._structures import *
from quippy.farray import fidentity, fzeros, frange, farray

import numpy as np

__all__ = quippy._structures.__all__ + ['orthorhombic_slab', 'rotation_matrix',
                                        'quartz_params', 'alpha_quartz', 'get_bulk_params',
                                        'alpha_quartz_cubic', 'get_quartz_params', 'get_bond_lengths']


def orthorhombic_slab(at, tol=1e-5, min_nrep=1, max_nrep=5, graphics=False, rot=None, periodicity=None, vacuum=None, shift=None, verbose=True):
    """Try to construct an orthorhombic cell equivalent to the
       primitive cell `at`, using supercells up to at most `max_nrep`
       repeats. Symmetry must be exact within a tolerance of `tol`. If
       `rot` is not None, we first transform `at` by the rotation
       matrix `rot`. The optional argument `periodicity` can be used to
       fix the periodicity one or more directions. It should be a three
       component vector with value zero in the unconstrained
       directions. The vector `vacuum` can be used to add vacuum in one
       or more directions. `shift` is a three component vector which
       can be used to shift the positions in the final cell. """

    def atoms_near_plane(at, n, d, tol=1e-5):
        """Return a list of atoms within a distance `tol` of the plane defined by np.dot(n, at.pos) == d"""
        pd = np.dot(n, at.pos) - d
        return (abs(pd) < tol).nonzero()[0]

    def sort_by_distance(at, ref_atom, dir, candidates):
        """Return a copy of `candidates` sorted by perpendicular distance from `ref_atom` in direction `dir`"""
        distances_candidates =  zip([at.pos[dir,i]-at.pos[dir,ref_atom] for i in candidates], candidates)
        distances_candidates.sort()
        return [p for (d, p) in distances_candidates]

    def orthorhombic_box(at):
        """Return a copy of `at` in an orthorhombic box surrounded by vacuum"""
        at = at.copy()
        at.map_into_cell()
        at.set_lattice([[2.0*(at.pos[1,:].max() - at.pos[1,:].min()), 0.0, 0.0],
                        [0.0, 2.0*(at.pos[2,:].max() - at.pos[2,:].min()), 0.0],
                        [0.0, 0.0, 2.0*(at.pos[3,:].max() - at.pos[3,:].min())]],
                       scale_positions=False)
        at.map_into_cell()
        return at

    def discard_outliers(at, indices, dirs, keep_fraction=0.5):
        """Return a copy of `indices` with the atoms with fractional coordinates along directions in `dirs`
           outside +/-`keep_fraction`/2 excluded. Lattice used is close fitting, `at.lattice`/2."""
        g = np.linalg.inv(at.lattice/2)
        t = np.dot(g, at.pos[:,indices])
        return indices[ np.logical_and(t[dirs,:] >= -keep_fraction/2.0, t[dirs,:] < keep_fraction/2.0).all(axis=1) ]

    def check_candidate_plane(at, ref_plane, cand_plane, dirs, verbose=False, label=''):
        """Check whether in-plane displacements of atoms listed in `ref_plane` match those of `cand_plane` in directions given by `dirs`"""

        # Which pair of planes has more atoms, reference or candidate?
        if len(ref_plane) < len(cand_plane):
            smaller = ref_plane
            larger  = cand_plane
        else:
            smaller = cand_plane
            larger  = ref_plane

        matches = {}
        for j in smaller:
            for k in larger:
                if at.z[k] == at.z[j] and abs(at.pos[dirs,k] - at.pos[dirs,j]).max() < tol:
                    matches[j] = k
                    break

        if verbose:
            print '   ', label, len(matches), '/', len(smaller), 'matches'

        return len(matches) == len(smaller)


    if rot is not None:
        at = transform(at, rot)

    xyz = fidentity(3)
    nrep = min_nrep-1
    max_dist = fzeros(3)

    if periodicity is not None:
        periodicity = farray(periodicity)
        periodicity = dict(zip((periodicity != 0).nonzero()[0], periodicity[periodicity != 0]))
    else:
        periodicity = {}

    if verbose:
        for (dir, p) in periodicity.iteritems():
            print 'Periodicity in direction %d fixed at %f' % (dir, p)

    if graphics:
        import atomeye
        viewer = atomeye.AtomEyeViewer()

    while sorted(periodicity.keys()) != [1,2,3]:
        nrep += 1
        if nrep > max_nrep:
            raise ValueError('Maximum size of supercell (%d) exceeded' % max_nrep)
        if verbose:
            print '\n\nSupercell %d' % nrep
        sup = supercell(at, nrep, nrep, nrep)
        box = orthorhombic_box(sup)
        box.pos[:] = box.pos - np.tile(box.pos.mean(axis=2), [box.n, 1]).T

        for dir in set([1,2,3]) - set(periodicity.keys()):

            if verbose:
                print '  Direction %d' % dir

            other_dirs = list(set([1,2,3]) - set([dir]))

            pos_index = zip(box.pos[dir,:],frange(box.n))
            pos_index.sort()

            # Find a pair of planes
            while pos_index:
                ref_pos1, ref_atom1 = pos_index.pop(0)

                # Find atom to define second plane
                while pos_index:
                    ref_pos2, ref_atom2 = pos_index.pop(0)
                    if abs(ref_pos2 - ref_pos1) > tol: break
                else:
                    continue

                ref_plane1 = atoms_near_plane(box, xyz[:,dir], box.pos[dir,ref_atom1], tol)
                ref_plane2 = atoms_near_plane(box, xyz[:,dir], box.pos[dir,ref_atom2], tol)

                # Only keep reference atoms in the centre of the cell
                ref_plane1 = discard_outliers(box, ref_plane1, other_dirs)
                ref_plane2 = discard_outliers(box, ref_plane2, other_dirs)

                if len(ref_plane1) > 2 and len(ref_plane2) > 2:
                    # Now we've got two planes, both with more than two atoms in them
                    break
            else:
                # Used up all atoms without finding two planes
                if verbose:
                    print '    No valid reference planes found.\n'
                continue

            if verbose:
                print '    Reference plane #1 through atom %d ' % ref_atom1
                print '    Reference plane #2 through atom %d at distance %r\n' % (ref_atom2, ref_pos2 - ref_pos1)

            if graphics:
                highlight = fzeros(box.n)
                highlight[ref_plane1] = 1
                highlight[ref_plane2] = 2
                box.add_property('highlight', highlight, overwrite=True)
                viewer.show(box, 'highlight')
                viewer.wait()
                raw_input('continue')

            candidates = [i for i in frange(box.n) if box.pos[dir,i] > box.pos[dir,ref_atom2] + max_dist[dir] + tol]
            candidates = sort_by_distance(box, ref_atom1, dir, candidates)

            while candidates:
                cand1 = candidates.pop(0)

                max_dist[dir] = box.pos[dir,cand1] - box.pos[dir,ref_atom1]

                if verbose:
                    print '    Trying plane through atom %d distance %r' % (cand1, max_dist[dir])

                cand_plane1 = atoms_near_plane(box, xyz[:,dir], box.pos[dir,cand1], tol)

                for cand2 in sort_by_distance(box, ref_atom1, dir, set(candidates) - set(cand_plane1)):
                    if abs((box.pos[dir,cand2] - box.pos[dir,cand1]) - (box.pos[dir,ref_atom2] - box.pos[dir,ref_atom1])) < tol:
                        if verbose:
                            print '    Found pair to plane, passing through atom %d distance %r ' % (cand2, box.pos[dir,cand2] - box.pos[dir,ref_atom1])
                        break
                else:
                    if verbose:
                        print '    Cannot find second candidate plane.\n'
                    candidates = sort_by_distance(box, ref_atom1, dir, set(candidates) - set(cand_plane1))
                    continue

                if graphics:
                    highlight[cand_plane1] = 3
                    box.highlight[:] = highlight
                    viewer.show(box, 'highlight')
                    viewer.wait()

                cand_plane2 = atoms_near_plane(box, xyz[:,dir], box.pos[dir,cand2], tol)

                if graphics:
                    highlight[cand_plane2] = 4
                    box.highlight[:] = highlight                    
                    viewer.show(box, 'highlight')
                    viewer.wait()

                    highlight[cand_plane1] = 0
                    highlight[cand_plane2] = 0

                # Remove cand_plane1 from list of candidates
                candidates = sort_by_distance(box, ref_atom1, dir, set(candidates) - set(cand_plane1))

                # Check ref_plane1 against cand_plane1 and ref_plane2 against cand_plane2 in directions
                # listed in other_dirs
                match1 = check_candidate_plane(box, ref_plane1, cand_plane1, other_dirs, verbose, 'Plane #1:')
                match2 = check_candidate_plane(box, ref_plane2, cand_plane2, other_dirs, verbose, 'Plane #2:')

                if match1 and match2:
                    periodicity[dir] = box.pos[dir,cand1] - box.pos[dir,ref_atom1]
                    if verbose:
                        print '\n  Periodicity in direction %d is %f\n' % (dir, box.pos[dir,cand1] - box.pos[dir,ref_atom1])

                    if graphics:
                        highlight[cand_plane1] = 3
                        highlight[cand_plane2] = 3
                        box.highlight[:] = highlight
                        viewer.show(box, 'highlight')
                        viewer.wait()
                        raw_input('continue...')
                    break

                if graphics:
                    raw_input('continue...')
            else:
                # Failed to find match for direction dir
                continue

    # Finally, construct new cell by selecting atoms within first unit cell
    lattice = farray(np.diag([periodicity[1], periodicity[2], periodicity[3]]))
    g = np.linalg.inv(lattice)

    nrepx, nrepy, nrepz = fit_box_in_cell(periodicity[1], periodicity[2], periodicity[3], at.lattice)

    sup = supercell(at, nrepx, nrepy, nrepz)
    sup.map_into_cell()

    # small shift to avoid coincidental cell alignments
    delta = np.tile([0.01, 0.02, 0.03], [sup.n, 1]).T
    if shift is not None and vacuum is not None:
        delta = delta + np.tile(shift, [sup.n, 1]).T
    t = np.dot(g, sup.pos) + delta

    orthorhombic = sup.select(np.logical_and(t >= -0.5, t < 0.5).all(axis=1))

    if vacuum:
        lattice = farray(np.diag(lattice.diagonal() + vacuum))

    if shift is not None and vacuum is None:
        if verbose:
            print 'Shifting positions by %s' % np.dot(lattice, shift)
        orthorhombic.pos += np.tile(np.dot(lattice, shift), [orthorhombic.n, 1]).T

    orthorhombic.set_lattice(lattice, scale_positions=False)
    orthorhombic.map_into_cell()
    return orthorhombic

def rotation_matrix(unit, y, z=None, x=None, tol=1e-5):
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


quartz_params = {'experiment': {'a': 4.9160,
                                'c': 5.4054,
                                'u': 0.4697,
                                'x': 0.4135,
                                'y': 0.2669,
                                'z': 0.1191},
                 'CASTEP_LDA': {'a': 4.87009,
                                'c': 5.36255,
                                'u': 0.46699,
                                'x': 0.41289,
                                'y': 0.27198,
                                'z': 0.11588},
                 'CASTEP_GGA': {'a': 5.02836,
                                'c': 5.51193,
                                'u': 0.48128,
                                'x': 0.41649,
                                'y': 0.24661,
                                'z': 0.13594},
                 'ASAP_JRK': {'a': 4.8403809707320216,
                              'c': 5.3285240037002248,
                              'u': 0.46417561617105912,
                              'x': 0.41174271054205958,
                              'y': 0.27872745399831672,
                              'z': 0.10973603276909905}
                 }

def alpha_quartz(a=4.9134,c=5.4052, u=0.4699, x=0.4141, y=0.2681, z=0.7854-2.0/3.0):
    """Primitive 9-atom orthorhombic alpha quartz cell"""

    from math import sqrt

    a1 = farray((0.5*a, -0.5*sqrt(3.0)*a, 0.0))
    a2 = farray((0.5*a,  0.5*sqrt(3.0)*a, 0.0))
    a3 = farray((0.0,    0.0,             c))

    lattice = fzeros((3,3))
    lattice[:,1] = a1
    lattice[:,2] = a2
    lattice[:,3] = a3

    at = Atoms(n=9,lattice=lattice)

    at.set_atoms((14,14,14,8,8,8,8,8,8))

    z += 2.0/3.0

    at.pos[:,1] =  u*a1 + 2.0/3.0*a3
    at.pos[:,2] =  u*a2 + 1.0/3.0*a3
    at.pos[:,3] = -u*a1 - u*a2
    at.pos[:,4] =  x*a1 + y*a2 + z*a3
    at.pos[:,5] = -y*a1 + (x-y)*a2  + (2.0/3.0 + z)*a3
    at.pos[:,6] = (y-x)*a1 - x*a2   + (1.0/3.0 + z)*a3
    at.pos[:,7] = y*a1 + x*a2 - z*a3
    at.pos[:,8] = -x*a1 + (y-x)*a2 + (2.0/3.0 - z)*a3
    at.pos[:,9] = (x - y)*a1 - y*a2 + (1.0/3.0 - z)*a3

    return at

def get_quartz_params(at):

    assert at.n == 9
    assert (at.z == 14).sum() == 3
    assert (at.z == 8).sum() == 6

    from quippy import get_lattice_params

    lat_params = get_lattice_params(at.lattice)
    a, c = lat_params[0], lat_params[2]
    print 'a      = ', a
    print 'c      = ', c
    print 'c/a    = ', c/a
    print 'V      = ', at.cell_volume()
    print 'V/SiO2 = ', at.cell_volume()/3.0

    frac_pos = np.dot(at.g, at.pos)
    u = frac_pos[1,1]
    x,y,z = frac_pos[:,4]
    z -= 2.0/3.0
    if z < 0.0: z += 1.0
    if z > 1.0: z -- 1.0

    print 'u      = ', u
    print 'x      = ', x
    print 'y      = ', y
    print 'z      = ', z

    return {'a':a, 'c':c, 'u':u, 'x':x, 'y':y, 'z':z}


def alpha_quartz_cubic(*args, **kwargs):
    """Non-primitive 18-atom cubic quartz cell."""

    from quippy import supercell

    a0 = alpha_quartz(*args, **kwargs)
    at = supercell(a0, 4, 4, 1)
    at.map_into_cell()

    lattice = fzeros((3,3))
    lattice[1,1] = a0.lattice[1,1]*2.0
    lattice[2,2] = a0.lattice[2,2]*2.0
    lattice[3,3] = a0.lattice[3,3]

    g = np.linalg.inv(lattice)
    t = np.dot(g, at.pos)
    cubic = at.select(np.logical_and(t >= -0.5, t < 0.5).all(axis=1))
    cubic.set_lattice(lattice)
    return cubic


def get_bond_lengths(at):
    """Return a dictionary mapping tuples (Species1, Species2) to an farray of bond-lengths"""
    at.calc_connect()
    r_ij = farray(0.0)
    res = {}
    for i in frange(at.n):
        for n in frange(at.n_neighbours(i)):
            j = at.neighbour(i, n, distance=r_ij)
            print i, j, at.z[i], at.z[j], r_ij
            minij, maxij = min((i,j)), max((i,j))
            key = (at.species[minij].stripstrings(), at.species[maxij].stripstrings())
            if not key in res: res[key] = []
            res[key].append(r_ij.astype(float))
    print res
    return dict((k,farray(v)) for (k,v) in res.iteritems())

def get_bulk_params(bulk, lat_type, verbose=True):
    """Return 6-tuple of lattice parameters a, c, u, x, y, z for
       cell `bulk` of lattice type `lat_type`"""
    a, b, c, alpha, beta, gamma = get_lattice_params(bulk.lattice)
    del b, alpha, beta, gamma
    u, x, y, z = (None, None, None, None)

    if lat_type in ('diamond', 'bcc', 'fcc'):
        if verbose:
            print '%s lattice, a=%.3f' % (lat_type, a)
    elif lat_type == 'anatase':
        u = bulk.pos[3, 5]/c
        if verbose:
            print 'anatase lattice, a=%.3f c=%.3f u=%.3f' % (a, c, u)
    elif lat_type == 'rutile':
        u = bulk.pos[1, 3]/a
        if verbose:
            print 'rutile lattice, a=%.3f c=%.3f u=%.3f' % (a, c, u)
    elif lat_type == 'alpha_quartz':
        qz = get_quartz_params(bulk)
        a, c, u, x, y, z = (qz['a'], qz['c'], qz['u'],
                            qz['x'], qz['y'], qz['z'])
        if verbose:
            print 'alpha_quartz lattice, ',
            print ('a=%.3f c=%.3f u=%.3f x=%.3f y=%.3f z=%.3f'
                   % (a, c, u, x, y ,z))
    else:
        raise ValueError('Unknown latttice type %s' % lat_type)

    return (a, c, u, x, y, z)
