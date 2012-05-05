# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

"""Definition of the Atoms class

This module defines the quippy.Atoms class which stores and
manipulates a collection of atoms"""

from quippy import _atoms
from quippy._atoms import *

import os, weakref, warnings, copy, sys
from quippy.farray import frange, farray, fzeros, fvar
from quippy.dictmixin import DictMixin
from quippy.util import infer_format, parse_slice
from quippy.extendable_str import Extendable_str
from quippy import QUIPPY_TRUE, QUIPPY_FALSE
from quippy import available_modules, FortranDerivedTypes
from math import pi
import numpy as np

__all__ = _atoms.__all__ + ['AtomsReaders',
                            'AtomsWriters',
                            'atoms_reader']
                                   
if 'ase' in available_modules:
    import ase
else:
    import quippy.miniase as ase

def get_lattice_params(lattice):
    a,b,c,alpha,beta,gamma = (farray(0.), farray(0.), farray(0.),
                              farray(0.), farray(0.), farray(0.))
    _atoms.get_lattice_params(lattice,a,b,c,alpha,beta,gamma)
    return tuple(float(x) for x in (a,b,c,alpha,beta,gamma))

get_lattice_params.__doc__ = """Wrapper around _atoms.get_lattice_params()

Returns parameters of `lattice` as 6-tuple (a,b,c,alpha,beta,gamma).

""" + _atoms.get_lattice_params.__doc__

class NeighbourInfo(object):
    """Store information about a single neighbour of an atom"""

    __slots__ = ('j', 'distance', 'diff', 'cosines', 'shift')

    def __init__(self, j, distance, diff, cosines, shift):
        self.j = j
        self.distance = float(distance)
        self.diff = diff.copy()
        self.cosines = cosines.copy()
        self.shift = shift.copy()

    def __int__(self):
        return self.j

    def __repr__(self):
        return ('NeighbourInfo(j=%d,distance=%f,diff=%s,cosines=%s,shift=%s)'
                    % (self.j, self.distance, self.diff,
                       self.cosines, self.shift))

class Neighbours(object):
    """Provides access to neighbours of an atom.

       Supports both iteration over all atoms, and indexing
       e.g. at.neighbours[1] returns a list of the neighbours of the
       atom with index 1. If at.fortran_indexing is True, atom and
       neighbour indices start from 1; otherwise they are numbered
       from zero."""

    def __init__(self, at, hysteretic=False):
        self.atref = weakref.ref(at)
        self.hysteretic = hysteretic

    def is_neighbour(self, i, j):
        return (i,j) in self.pairs()

    def pairs(self):
        """Yield pairs of atoms (i,j) with i < j which are neighbours"""
        at = self.atref()
        for i, neighbour_list in zip(at.indices, self.iterneighbours()):
            for neighb in neighbour_list:
                if i < neighb.j:
                    yield (i,neighb.j)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        # Neighbours are considered to be equal if *topology* matches,
        # not distances, displacement vectors and shifts.
        return sorted(self.pairs()) == sorted(other.pairs())

    def __ne__(self, other):
        return not self.__eq__(other)

    def __iter__(self):
        return self.iterneighbours()

    def iterneighbours(self):
        """Iterate over the neighbours of all atoms"""
        at = self.atref()
        for i in at.indices:
            yield self[i]

    def __getitem__(self, i):
        at = self.atref()

        if self.hysteretic:
            connect = at.hysteretic_connect
        else:
            connect = at.connect

        if not connect.initialised:
            if self.hysteretic:
                at.calc_connect_hysteretic(connect)
            else:
                at.calc_connect(connect)

        distance = farray(0.0)
        diff = fzeros(3)
        cosines = fzeros(3)
        shift = fzeros(3, dtype=np.int32)

        res = []
        if not at.fortran_indexing:
            i = i+1 # convert to 1-based indexing

        for n in frange(at.n_neighbours(i, alt_connect=connect)):
            j = at.neighbour(i, n, distance, diff,
                             cosines, shift, alt_connect=connect)
            if not at.fortran_indexing:
                j = j-1
            res.append(NeighbourInfo(j, distance, diff, cosines, shift))

        if at.fortran_indexing:
            res = farray(res) # to give 1-based indexing
        return res

    def distances(self, Z1=None, Z2=None):
        """Distances between pairs of neighbours, optionally
           filtered by species (Z1,Z2)"""
        at = self.atref()
        for i in at.indices:
            for neighb in self[i]:
                if neighb.j > i: continue
                if Z1 is not None and Z2 is not None:
                    if sorted((at.z[i], at.z[neighb.j])) == sorted((Z1, Z2)):
                        yield neighb.distance
                else:
                    yield neighb.distance

class PropertiesWrapper(DictMixin):
    """Wrapper between quippy properties and ASE arrays"""

    def __init__(self, at):
        DictMixin.__init__(self)
        self.atref = weakref.ref(at)

    def __getitem__(self, key):
        at = self.atref()
        res = at.properties[at.name_map.get(key, key)].view(np.ndarray).T
        if res.dtype == 'int32':   # convert dtype int32 -> int
            warnings.warn('Making copy of arrays["%s"] since quippy/ASE dtypes incompatible' % key)
            res = res.astype(int)
        return res

    def __setitem__(self, key, value):
        at = self.atref()
        key = Atoms.name_map.get(key, key)
        value = np.array(value)

        if value.shape[0] != at.n:
            warnings.warn(('Assignment to arrays["%s"]: changing size '+
                           ' from %d to %d') % (key, at.n, value.shape[0]))
            for p in at.properties.keys():
                at.remove_property(p)
            at.n = value.shape[0]
            at.nbuffer = at.n
            at.ndomain = at.n

        if value.dtype.kind != 'S':
            value = value.T

        if not at.has_property(key):
            at.add_property(key, value)
        else:
            at.properties[key][...] = value

    def __delitem__(self, key):
        at = self.atref()
        warnings.warn('Deletion of arrays["%s"]' % key)
        key = at.name_map.get(key, key)
        at.remove_property(key)

    def keys(self):
        at = self.atref()
        return [at.rev_name_map.get(key, key) for key in at.properties.keys()]


class Atoms(_atoms.Atoms, ase.Atoms):

    """
    Pythonic wrapper over auto-generated Atoms class.

    Also inherits from ase.Atoms so has all ASE Atoms methods

    An atoms object contains atomic numbers, all dynamical variables
    and connectivity information for all the atoms in the simulation
    cell."""

    _cmp_skip_fields = ['own_this', 'ref_count', 'domain']

    name_map = {'positions'       : 'pos',
                'numbers'         : 'Z',
                'charges'         : 'charge'}

    rev_name_map = dict(zip(name_map.values(), name_map.keys()))

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None,
                 momenta=None, masses=None, magmoms=None, charges=None,
                 scaled_positions=None, cell=None, pbc=None, constraint=None,
                 calculator=None, info=None, n=None, lattice=None,
                 properties=None, params=None, fixed_size=None,
                 fortran_indexing=True, fpointer=None, finalise=True,
                 **readargs):

        # check for mutually exclusive options
        if cell is not None and lattice is not None:
            raise ValueError('only one of cell and lattice can be present')

        if n is None:
            n = 0
        if cell is not None:
            lattice = cell.T
        if lattice is None:
            lattice = np.eye(3)

        from quippy import Dictionary
        if properties is not None and not isinstance(properties, Dictionary):
            properties = Dictionary(properties)
        if params is not None and not isinstance(params, Dictionary):
            params = Dictionary(params)

        _atoms.Atoms.__init__(self, n=n, lattice=lattice,
                              properties=properties,
                              params=params, fixed_size=fixed_size,
                              fortran_indexing=fortran_indexing,
                              fpointer=fpointer, finalise=finalise)

        self._ase_arrays = PropertiesWrapper(self)
        self.neighbours = Neighbours(self)
        self.hysteretic_neighbours = Neighbours(self, hysteretic=True)

        # If first argument is quippy.Atoms instance, copy data from it
        if symbols is not None and isinstance(symbols, self.__class__):
            self.copy_from(symbols)
            symbols = None

        # Try to read from first argument
        if symbols is not None:
            try:
                self.read_from(symbols, **readargs)
                symbols = None
            except IOError:
                pass

        ## ASE compatibility
        remove_properties = []

        if symbols is None and numbers is None:
            if self.has_property('z'):
                numbers = self.z.view(np.ndarray)
            else:
                numbers = [0]*len(self)
                remove_properties.append('Z')

        if symbols is None and positions is None:
            if self.has_property('pos'):
                positions = self.pos.view(np.ndarray).T
            else:
                remove_properties.append('pos')

        # if symbols is None and momenta is None:
        #     if hasattr(self, 'velo'):
        #         momenta = self.velo.view(np.ndarray).T

        if symbols is None and cell is None:
            cell = self.lattice.T.view(np.ndarray)
        if symbols is None and pbc is None:
            pbc = self.get_pbc()
        if charges is None and self.has_property('charge'):
            charges = self.charge.view(np.ndarray)
            
        ase.Atoms.__init__(self, symbols, positions, numbers,
                           tags, momenta, masses, magmoms, charges,
                           scaled_positions, cell, pbc, constraint,
                           calculator)

        # remove anything that ASE added that we don't want
        for p in remove_properties:
            self.remove_property(p)

        ## end ASE compatibility

        if self.has_property('z') and not self.has_property('species'):
            self.add_property('species', ' '*10)
            if self.n != 0 and not (self.z == 0).all():
                self.set_atoms(self.z) # initialise species from z

        if info is not None:
            self.params.update(info)

        self._initialised = True

    def new_array(self, name, a, dtype=None, shape=None):
        # we overrride ase.Atoms.new_array() to allow "special" arrays
        # like "numbers", "positions" to be added more than once without
        # raising a RuntimeError
        if name in self.name_map and name in self.arrays:
            self.arrays[name] = a
            return
        ase.Atoms.new_array(self, name, a, dtype, shape)

    def set_lattice(self, lattice, scale_positions=False):
        """Change the lattice vectors, keeping the inverse lattice vectors
           up to date. Optionally map the existing atoms into the new cell
           and recalculate connectivity (by default scale_positions=False)."""
        _atoms.Atoms.set_lattice(self, lattice, scale_positions)

    def _get_cell(self):
        """Get ASE cell from QUIP lattice"""
        return self.lattice.view(np.ndarray).T

    def _set_cell(self, cell):
        """Set QUIP lattice from ASE cell"""
        self.set_lattice(cell.T, scale_positions=False)

    _cell = property(_get_cell, _set_cell)

    def _set_pbc(self, pbc):
        self.is_periodic = np.array(pbc).astype(int)

    def _get_pbc(self):
        return self.is_periodic.view(np.ndarray) == QUIPPY_TRUE

    _pbc = property(_get_pbc, _set_pbc)

    def _get_ase_arrays(self):
        """Provides access to ASE arrays, stored in QUIP properties dict"""

        return self._ase_arrays

    def _set_ase_arrays(self, value):
        """Set ASE arrays. Does not remove existing QUIP properties."""

        self._ase_arrays.update(value)

    arrays = property(_get_ase_arrays, _set_ase_arrays)

    def _get_info(self):
        """Get ASE info dictionary

        Entries are actually storred in QUIP params dictionary"""

        return self.params

    def _set_info(self, value):
        """Set ASE info dictionary.

         Updates entries in QUIP params dictionary. Note that assigning {} to
         info doesn't remove the existing params."""

        self.params.update(value)

    info = property(_get_info, _set_info)

    def iterindices(self):
        """Iterate over atoms.

        If fortran_indexing is True, returns range 1..self.n.
        Otherwise, returns range 0..(self.n-1)."""

        if self.fortran_indexing:
            return iter(frange(len(self)))
        else:
            return iter(range(len(self)))

    indices = property(iterindices)

    def iteratoms(self):
        """Iterate over atoms, calling get_atom() for each one"""
        for i in self.indices:
            yield self.get_atom(i)

    def equivalent(self, other):
        """Test for equivalence of two Atoms objects.

        Equivalence is less strong than equality.  Equality (written
        `self == other`) requires all properties and parameters to be
        equal. Equivalence requires only that the number of atoms,
        positions, atomic numbers, unit cell and periodic boundary
        conditions match.

        Note: the quippy expression a.equivalent(b) has the same
        definition as a == b in ASE. This means that a quippy.Atoms
        instance can be compared with an ase.Atoms instance using
        this method.
        """

        try:
            a = self.arrays
            b = other.arrays
            return (len(self) == len(other) and
                    (a['positions'] == b['positions']).all() and
                    (a['numbers'] == b['numbers']).all() and
                    (self._cell == other.cell).all() and
                    (self._pbc == other.pbc).all())
        except AttributeError:
            return False

    @classmethod
    def read(cls, source, format=None, **kwargs):
        """Read Atoms object from file `source` according to `format`.

        If `format` is None, filetype is inferred from filename.
        Returns a new Atoms instance; to read into an existing Atoms
        object, use the read_from() method."""

        if isinstance(source, basestring) and '@' in source:
            source, frame = source.split('@')
            frame = parse_slice(frame)
            if 'frame' in kwargs:
                raise ValueError("Conflicting frame references given: kwarg frame=%r and @-reference %r" %
                                 (kwargs['frame'], frame))
            if not isinstance(frame, int):
                raise ValueError("Frame @-reference %r does not resolve to single frame" % frame)
            kwargs['frame'] = frame

        filename, source, format = infer_format(source, format, AtomsReaders)

        opened = False
        if format in AtomsReaders:
            source = AtomsReaders[format](source, **kwargs)
            opened = True

        if isinstance(source, basestring):
            raise IOError("Don't know how to read from file '%s'" % source)
        if not hasattr(source, '__iter__'):
            raise IOError('Cannot read from %r - not an iterator' % source)

        at = iter(source).next()
        if not isinstance(at, cls):
            raise ValueError('Object %r read from  %r is not Atoms instance'
                             % (at, source))
        if opened and hasattr(source, 'close'):
            source.close()
        if filename is not None:
            at.filename = filename
        return at


    def write(self, dest=None, format=None, properties=None,
              prefix=None, **kwargs):
        if dest is None:
            # if filename is missing, save back to file from
            # which we loaded configuration
            if hasattr(self, 'filename'):
                dest = self.filename
            else:
                raise ValueError("No 'dest' and Atoms has no stored filename")

        filename, dest, format = infer_format(dest, format, AtomsWriters)
        opened = filename is not None

        if format in AtomsWriters:
            dest = AtomsWriters[format](dest, **kwargs)

        if not hasattr(dest, 'write'):
            raise ValueError('Don\'t know how to write to "%s" in format "%s"'
                             % (dest, format))

        write_kwargs = {}
        if properties is not None:
            write_kwargs['properties'] = properties
        if prefix is not None:
            write_kwargs['prefix'] = prefix
        try:
            res = dest.write(self, **write_kwargs)
        except TypeError:
            raise ValueError('destination %r doesn\'t support arguments %r'
                             % (dest, write_kwargs))

        if opened and hasattr(dest, 'close'):
            dest.close()
        return res

    def select(self, mask=None, list=None, orig_index=None):
        """Select a subset of the atoms in an Atoms object

        Use either a logical mask array or a list of atom indices to include.

        small_at = at.select([mask, list])

        """
        if mask is not None:
            mask = farray(mask)
            out = self.__class__(n=mask.sum(), lattice=self.lattice, properties={}, params={})
            _atoms.Atoms.select(out, self, mask=mask, orig_index=orig_index)
        elif list is not None:
            list = farray(list)
            out = self.__class__(n=len(list), lattice=self.lattice)
            _atoms.Atoms.select(out, self, list=list, orig_index=orig_index)
        else:
            raise ValueError('Either mask or list must be present.')
        return out

    def copy(self):
        other = self.__class__(n=self.n, lattice=self.lattice,
                               properties=self.properties, params=self.params)

        # copy any normal (not Fortran) attributes
        for k, v in self.__dict__.iteritems():
            if not k.startswith('_') and k not in other.__dict__:
                other.__dict__[k] = v

        # from _atoms.Atoms
        other.use_uniform_cutoff = self.use_uniform_cutoff
        other.cutoff = self.cutoff
        other.cutoff_break = self.cutoff_break
        other.nneightol = self.nneightol

        # from ase.Atoms
        other.constraints = copy.deepcopy(self.constraints)
        other.adsorbate_info = copy.deepcopy(self.adsorbate_info)
        return other

    def copy_from(self, other):
        """Replace contents of this Atoms object with data from `other`."""
        self.__class__.__del__(self)
        self.__init__(n=other.n, lattice=other.lattice,
                      properties=other.properties, params=other.params)

        # copy any normal (not Fortran) attributes
        for k, v in other.__dict__.iteritems():
            if not k.startswith('_') and k not in self.__dict__:
                self.__dict__[k] = v

        # from _atoms.Atoms
        self.use_uniform_cutoff = other.use_uniform_cutoff
        self.cutoff = other.cutoff
        self.cutoff_break = other.cutoff_break
        self.nneightol = other.nneightol

        # from ase.Atoms
        self.constraints = copy.deepcopy(other.constraints)
        self.adsorbate_info = copy.deepcopy(other.adsorbate_info)

    def read_from(self, source, **readargs):
        """Replace contents of this Atoms object with file `source`"""
        if isinstance(source, self.__class__):
            self.copy_from(source)
        else:
            tmp = Atoms.read(source, **readargs)
            self.shallow_copy_from(tmp)
            # tmp goes out of scope here, but reference counting
            # prevents it from being free'd.

    def __len__(self):
        return self.n

    def __getattr__(self, name):
        #print 'getattr', name
        #if name in self.properties:
        try:
            return self.properties[name]
        except KeyError:
            try:
                return self.params[name]
            except KeyError:
                raise AttributeError('Unknown Atoms attribute %s' % name)

    def __setattr__(self, name, value):
        #print 'setattr', name, value
        if not '_initialised' in self.__dict__:
            object.__setattr__(self, name, value)
        elif self.properties._fpointer is not None and name in self.properties:
            self.properties[name][...] = value
        elif self.params._fpointer is not None and name in self.params:
            if self.params.is_array(name):
                self.params[name][...] = value
            else:
                self.params[name] = value
        else:
            object.__setattr__(self, name, value)

    def md5_hash(self, ndigits):
        """Hash an atoms object with a precision of ndigits decimal
        digits.  Atomic numbers, lattice and fractional positions are
        fed to MD5 to form the hash."""

        def rounded_string_rep(a, ndigits):
            return np.array2string(a, precision=ndigits, suppress_small=True).replace('-0. ', ' 0. ')

        # Compute fractional positions, round them to ndigits, then sort them
        # for hash stability
        flat_frac_pos = np.dot(self.g,self.pos).flatten()
        flat_frac_pos.sort()

        # md5 module deprecated in Python 2.5 and later
        try:
           import hashlib
           md5 = hashlib.md5
        except ImportError:
           import md5
           md5 = md5.new

        m = md5()
        m.update(rounded_string_rep(self.lattice.flatten(), ndigits))
        m.update(str(self.z))
        m.update(rounded_string_rep(flat_frac_pos, ndigits))

        return m.hexdigest()

    def __hash__(self):
        return hash(self.md5_hash(4))

    #def __getitem__(self, i):
        # we override ase.Atoms.__getitem__ so we can raise
        # exception if we're using fortran indexing
    #    if self.fortran_indexing:
    #        raise RuntimeError('Atoms[i] inconsistent with fortran indexing')
    #    return ase.Atoms.__getitem__(self, i)

    def get_atom(self, i):
        """Return a dictionary containing the properties of the atom with
           index `i`. If fortran_indexing=True (the default), `i` should be in
           range 1..self.n, otherwise it should be in range 0..(self.n-1)."""
        if (self.fortran_indexing and (i < 1 or i > self.n)) or \
            (not self.fortran_indexing and (i < 0 or i > self.n-1)):
            raise IndexError('Atoms index out of range')
        atom = {}
        atom['_index'] = i
        atom['atoms'] = self
        for k in self.properties.keys():
            v = self.properties[k][...,i]
            if isinstance(v,np.ndarray):
                if v.dtype.kind == 'S':
                    v = ''.join(v).strip()
                elif v.shape == ():
                    v = v.item()
            atom[k.lower()] = v
        return atom

    def print_atom(self, i):
        """Pretty-print the properties of the atom with index `i`"""
        at = self.get_atom(i)
        title = 'Atom %d' % at['_index']
        title = title + '\n' + '-'*len(title)+'\n\n'
        fields = ['%-15s =  %s' % (k,at[k]) for k in sorted(at.keys())
                                            if k not in ['_index', 'atoms']]
        print title+'\n'.join(fields)

    def density(self):
        from quippy import ElementMass, N_A, MASSCONVERT

        """Density in units of :math:`g/m^3`. If `mass` property exists,
           use that, otherwise we use `z` and ElementMass table."""
        if self.has_property('mass'):
            mass = sum(self.mass)/MASSCONVERT/1.0e3
        else:
            mass = sum(ElementMass[z] for z in self.z)/MASSCONVERT/1.0e3

        return mass/(N_A*self.cell_volume()*1.0e-30)/1.0e3


    def add_property(self, name, value, n_cols=None,
                     overwrite=None, property_type=None):
        """
        Add a new property to this Atoms object.

        `name` is the name of the new property and `value` should be
        either a scalar or an array representing the value, which should
        be either integer, real, logical or string.

        If a scalar is given for `value` it is copied to every element
        in the new property.  `n_cols` can be specified to create a 2D
        property from a scalar initial value - the default is 1 which
        creates a 1D property.

        If an array is given for `value` it should either have shape
        (self.n,) for a 1D property or (n_cols,self.n) for a 2D
        property.  In this case `n_cols` is inferred from the shape of
        the `value` and shouldn't be passed as an argument.

        If `property_type` is present, then no attempt is made to
        infer the type from `value`. This is necessary to resolve
        ambiguity between integer and logical types.

        If property with the same type is already present then no error
        occurs.If `overwrite` is true, the value will be overwritten with
        that given in `value`, otherwise the old value is retained.
        """

        kwargs = {}
        if n_cols is not None: kwargs['n_cols'] = n_cols
        if overwrite is not None: kwargs['overwrite'] = overwrite

        if property_type is None:
            _atoms.Atoms.add_property(self, name, value, **kwargs)

        else:
            # override value_ref if property_type is specified

            from quippy import (T_INTEGER_A, T_REAL_A, T_LOGICAL_A, T_CHAR_A,
                                T_INTEGER_A2, T_REAL_A2, TABLE_STRING_LENGTH)

            new_property = not self.has_property(name)

            type_to_value_ref = {
               T_INTEGER_A  : 0,
               T_REAL_A : 0.0,
               T_CHAR_A  : " "*TABLE_STRING_LENGTH,
               T_LOGICAL_A : False,
               T_INTEGER_A2 : 0,
               T_REAL_A2: 0.0
               }
            try:
                value_ref = type_to_value_ref[property_type]
            except KeyError:
                raise ValueError('Unknown property_type %d' % property_type)

            if (hasattr(value, 'shape') and len(value.shape) == 2 and
                property_type != T_CHAR_A and n_cols is None):
                kwargs['n_cols'] = value.shape[0]

            _atoms.Atoms.add_property(self, name, value_ref, **kwargs)
            if new_property or overwrite:
                getattr(self, name.lower())[:] = value

    def __getstate__(self):
        es = Extendable_str()
        self.write('', estr=es, format='xyz')
        return str(es)

    def __setstate__(self, state):
        es = Extendable_str(state)
        self.read_from('', estr=es, format='xyz')

    def mem_estimate(self):
        """Estimate memory usage of this Atoms object, in bytes"""

        sizeof_table = 320
        mem = sum([p.itemsize*p.size for p in self.properties.values()])
        if self.connect.initialised:
            c = self.connect
            mem += sizeof_table*self.n*2 # neighbour1 and neighbour2 tables
            mem += 32*c.n_neighbours_total() # neighbour data
            mem += c.cell_heads.size*c.cell_heads.itemsize # cell data

        return mem



FortranDerivedTypes['type(atoms)'] = Atoms

AtomsReaders = {}
AtomsWriters = {}

def atoms_reader(source):
    """Decorator to add a new reader"""
    def decorate(func):
        global AtomsReaders
        if not source in AtomsReaders:
            AtomsReaders[source] = func
        return func
    return decorate
