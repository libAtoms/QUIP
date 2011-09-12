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
from quippy.structures import *
from quippy.farray import fidentity, fzeros, frange, farray

import numpy as np
try:
    import atomeye
except ImportError:
    pass

__all__ = ['orthorhombic_slab', 'rotation_matrix']

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
                box.show(highlight)
                atomeye.view.wait()
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
                    box.show(highlight)
                    atomeye.view.wait()

                cand_plane2 = atoms_near_plane(box, xyz[:,dir], box.pos[dir,cand2], tol)

                if graphics:
                    highlight[cand_plane2] = 4
                    box.show(highlight)
                    atomeye.view.wait()

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
                        box.show(highlight)
                        atomeye.view.wait()
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
