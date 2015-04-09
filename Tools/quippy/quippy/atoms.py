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

import os
import weakref
import copy
import sys
import logging
from math import pi, sqrt

import numpy as np

from quippy import _atoms
from quippy._atoms import *

from quippy.oo_fortran import update_doc_string
from quippy.farray import frange, farray, fzeros, fvar
from quippy.dictmixin import DictMixin
from quippy.util import infer_format, parse_slice
from quippy.units import N_A, MASSCONVERT
from quippy.periodictable import ElementMass
from quippy.extendable_str import Extendable_str
from quippy import QUIPPY_TRUE, QUIPPY_FALSE
from quippy import available_modules, FortranDerivedTypes, get_fortran_indexing

atomslog = logging.getLogger('quippy.atoms')

__doc__ = (_atoms.__doc__ + """

This module defines the :class:`Atoms`, which stores and manipulates a
collection of atoms, as well as the :class:`Connection` class which stores
topology and neighbour lists, and the :class:`DomainDecomposition` class.

""")

__all__ = _atoms.__all__ + ['NeighbourInfo', 'get_lattice_params_', 'get_lattice_params']
                                   
if 'ase' in available_modules:
    import ase
else:
    import quippy.miniase as ase

if 'phonopy' in available_modules:
    from phonopy.structure.atoms import Atoms as PhonopyAtoms

get_lattice_params_ = get_lattice_params

def get_lattice_params(lattice):
    """
    Wrapper around Fortran :func:`get_lattice_params_`

    Returns parameters of `lattice` as 6-tuple (a,b,c,alpha,beta,gamma).
    """
    
    a,b,c,alpha,beta,gamma = (farray(0.), farray(0.), farray(0.),
                              farray(0.), farray(0.), farray(0.))
    _atoms.get_lattice_params(lattice,a,b,c,alpha,beta,gamma)
    return tuple(float(x) for x in (a,b,c,alpha,beta,gamma))

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

class Connection(_atoms.Connection):
    __doc__ = update_doc_string(_atoms.Connection.__doc__, """    
    The :class:`Connection` is a subclass of :class:`_atoms.Connection`
    which adds supports for iteration over all atoms, and indexing
    e.g. ``at.connect.neighbours[1]`` returns a list of the neighbours of
    the atom with index 1.

    When indexed with an integer from 1 to `at.n`, returns an array
    of :class:`NeighbourInfo` objects, each of which corresponds to
    a particular pair `(i,j)` and has attributes `j`, `distance`,
    `diff`, `cosines` and `shift`.

    If ``fortran_indexing`` is True, atom and neighbour indices
    start from 1; otherwise they are numbered from zero.

    If connectivity information has not already been calculated
    :meth:`calc_connect` will be called automatically. The code to
    loop over the neighbours of all atoms is quite idiomatic::

        for i in at.indices:
           for neighb in at.connect[i]:
               print (neighb.j, neighb.distance, neighb.diff,
                      neighb.cosines, neighb.shift)

    Note that this provides a more Pythonic interface to the atomic
    connectivity information than the wrapped Fortran functions
    :meth:`Atoms.n_neighbours` and :meth:`Atoms.neighbour`.
    """)

    # def __init__(self, n=None, nbuffer=None, pos=None,
    #              lattice=None, g=None, origin=None,
    #              extent=None, nn_guess=None, fill=None,
    #              fpointer=None, finalise=True):
    #     _atoms.Connection.__init__(self, n, nbuffer, pos,
    #                                lattice, g, origin,
    #                                extent, nn_guess, fill,
    #                                fpointer, finalise)
        

    def is_neighbour(self, i, j):
        return (i,j) in self.pairs()

    def pairs(self):
        """Yield pairs of atoms (i,j) with i < j which are neighbours"""
        for i, neighbour_list in zip(self.parent.indices, self.iterneighbours()):
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
        for i in self.parent.indices:
            yield self[i]

    def __getitem__(self, i):
        if not self.initialised:
            if self is self.parent.hysteretic_connect:
                self.calc_connect_hysteretic(self.parent)
            else:
                self.calc_connect(self.parent)

        distance = farray(0.0)
        diff = fzeros(3)
        cosines = fzeros(3)
        shift = fzeros(3, dtype=np.int32)

        res = []
        if not get_fortran_indexing():
            i = i+1 # convert to 1-based indexing

        for n in frange(self.n_neighbours(i)):
            j = self.neighbour(self.parent, i, n, distance, diff,
                               cosines, shift)
            if not get_fortran_indexing():
                j = j-1
            res.append(NeighbourInfo(j, distance, diff, cosines, shift))

        if get_fortran_indexing():
            res = farray(res) # to give 1-based indexing
        return res

    def distances(self, Z1=None, Z2=None):
        """Distances between pairs of neighbours, optionally
           filtered by species (Z1,Z2)"""
        for i in self.parent.indices:
            for neighb in self[i]:
                if neighb.j > i: continue
                if Z1 is not None and Z2 is not None:
                    if sorted((self.parent.z[i], self.parent.z[neighb.j])) == sorted((Z1, Z2)):
                        yield neighb.distance
                else:
                    yield neighb.distance

    def get_neighbours(self, i):
        """
        Return neighbours of atom i

        Return arrays of indices and offsets to neighbouring
        atoms. The positions of the neighbor atoms can be calculated like
        this::

            indices, offsets = atoms.connect.get_neighbors(42)
            for i, offset in zip(indices, offsets):
               print atoms.positions[i] + dot(offset, atoms.get_cell())
               
        Compatible with ase.calculators.neighborlist.NeighborList.get_neighbors(),
        providing that NeighborList is constructed with bothways=True and
        self_interaction=False.
        """
        neighbours = self[i]
        indices = np.array([n.j for n in neighbours])
        offsets = np.r_[[n.shift for n in neighbours]]
        return (indices, offsets)


    def get_neighbors(self, i):
        """
        Variant spelling of :meth:`get_neighbours`
        """
        return self.get_neighbours(i)

 
class PropertiesWrapper(DictMixin):
    """Wrapper between quippy properties and ASE arrays"""

    def __init__(self, at):
        DictMixin.__init__(self)
        self.atref = weakref.ref(at)

    def __getitem__(self, key):
        at = self.atref()
        res = at.properties[at.name_map.get(key, key)].view(np.ndarray).T
        #if res.dtype == 'int32':   # convert dtype int32 -> int
        #    res = res.astype(int)
        #    atomslog.debug('Making copy of arrays["%s"] since quippy/ASE dtypes incompatible' % key)
        return res

    def __setitem__(self, key, value):
        at = self.atref()
        atomslog.debug('Setting arrays["%s"]' % key)
        key = Atoms.name_map.get(key, key)
        value = np.array(value)

        if value.shape[0] != at.n:
            atomslog.debug(('Assignment to arrays["%s"]: changing size '+
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
        atomslog.debug('Deletion of arrays["%s"]' % key)
        key = at.name_map.get(key, key)
        at.remove_property(key)

    def keys(self):
        at = self.atref()
        return [at.rev_name_map.get(key, key) for key in at.properties.keys()]


class Atoms(_atoms.Atoms, ase.Atoms):
    __doc__ = update_doc_string(_atoms.Atoms.__doc__, """
    The :class:`Atoms` class is a Pythonic wrapper over the auto-generated
    :class:`quippy._atoms.Atoms` class. Atoms object are usually
    constructed either by reading from an input file in one of the
    :ref:`fileformats`, or by using the structure creation functions in
    the :mod:`quippy.structures` or :mod:`ase.lattice` modules.

    For example to read from an :ref:`extendedxyz` file, use::

       from quippy.atoms import Atoms
       atoms = Atoms('filename.xyz')

    Or, to create an 8-atom bulk diamond cubic cell of silicon::

       from quippy.structures import diamond
       si_bulk = diamond(5.44, 14)

    The :class:`Atoms` class is inherited from the
    :class:`ase.atoms.Atoms` so has all the ASE Atoms attributes and
    methods. This means that quippy and ASE Atoms objects are fully
    interoperable.""", signature='Atoms([symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell, pbc, constraint, calculator, info, n, lattice, properties, params, fixed_size, **read_args])')

    _cmp_skip_fields = ['own_this', 'ref_count', 'domain',
                        'connect', 'hysteretic_connect', 'source']

    name_map = {'positions'       : 'pos',
                'numbers'         : 'Z',
                'charges'         : 'charge'}

    rev_name_map = dict(zip(name_map.values(), name_map.keys()))
    
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None,
                 momenta=None, masses=None, magmoms=None, charges=None,
                 scaled_positions=None, cell=None, pbc=None, constraint=None,
                 calculator=None, info=None, n=None, lattice=None,
                 properties=None, params=None, fixed_size=None,
                 fpointer=None, finalise=True,
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
                              fpointer=fpointer, finalise=finalise)

        self._ase_arrays = PropertiesWrapper(self)

        # If first argument is quippy.Atoms instance, copy data from it
        if isinstance(symbols, self.__class__):
            self.copy_from(symbols)
            symbols = None

        # Phonopy compatibility
        if 'phonopy' in available_modules:
            if symbols is not None and isinstance(symbols, PhonopyAtoms):
                atoms = symbols
                symbols = atoms.get_chemical_symbols()
                cell = atoms.get_cell()
                positions = atoms.get_positions()
                masses = atoms.get_masses()

        # Try to read from first argument, if it's not ase.Atoms
        if symbols is not None and not isinstance(symbols, ase.Atoms):
            self.read_from(symbols, **readargs)
            symbols = None

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

        # Make sure argument to ase.Atoms constructor are consistent with
        # properties already present in this Atoms object
        if symbols is None and momenta is None and self.has_property('momenta'):
            momenta = self.get_momenta()
        if symbols is None and masses is None and self.has_property('masses'):
            masses = self.get_masses()
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

        if isinstance(symbols, ase.Atoms):
            self.copy_from(symbols)

        ## end ASE compatibility

        if self.has_property('z') and not self.has_property('species'):
            self.add_property('species', ' '*10)
            if self.n != 0 and not (self.z == 0).all():
                self.set_atoms(self.z) # initialise species from z

        if info is not None:
            self.params.update(info)

        self._initialised = True

        # synonyms for backwards compatibility
        self.neighbours = self.connect
        self.hysteretic_neighbours = self.hysteretic_connect
        

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
        """ASE info dictionary

        Entries are actually stored in QUIP params dictionary."""

        return self.params

    def _set_info(self, value):
        """Set ASE info dictionary.

        Entries are actually stored in QUIP params dictionary.
        Note that clearing Atoms.info doesn't empty params
        """

        self.params.update(value)

    info = property(_get_info, _set_info)

    def _indices(self):
        """Return array of atoms indices

        If global ``fortran_indexing`` is True, returns FortranArray containing
        numbers 1..self.n.  Otherwise, returns a standard numpuy array
        containing numbers in range 0..(self.n-1)."""

        if get_fortran_indexing():
            return farray(list(frange(len(self))))
        else:
            return np.array(list(range(len(self))))

    indices = property(_indices)

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

        .. note:: 

            The quippy expression a.equivalent(b) has the same
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
        """
        Class method to read Atoms object from file `source` according to `format`

        If `format` is None, filetype is inferred from filename.
        Returns a new Atoms instance; to read into an existing Atoms
        object, use the read_from() method.

        If `source` corresponds to a known format then it used
        to construct an appropriate iterator from the :attr:`AtomsReaders`
        dictionary. See :ref:`fileformats` for a list of supported
        file formats. 
        
        If `source` corresponds to an unknown format then it is
        expected to be an iterator returning :class:`Atoms` objects.
        """

        if isinstance(source, basestring) and '@' in os.path.basename(source):
            source, frame = source.split('@')
            if source.endswith('.db'):
                source = source+'@'+frame
                format = 'db'
            else:
                frame = parse_slice(frame)
                if 'frame' in kwargs:
                    raise ValueError("Conflicting frame references given: kwarg frame=%r and @-reference %r" %
                                     (kwargs['frame'], frame))
                if not isinstance(frame, int):
                    raise ValueError("Frame @-reference %r does not resolve to single frame" % frame)
                kwargs['frame'] = frame

        from quippy.io import AtomsReaders
        filename, source, format = infer_format(source, format, AtomsReaders)

        opened = False
        if format in AtomsReaders:
            source = AtomsReaders[format](source, format=format, **kwargs)
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
        """
        Write this :class:`Atoms` object to `dest`. If `format` is
        absent it is inferred from the file extension or type of `dest`,
        as described for the :meth:`read` method.  If `properties` is
        present, it should be a list of property names to include in the
        output file, e.g. `['species', 'pos']`.
        
        See :ref:`fileformats` for a list of supported file formats.
        """
        
        if dest is None:
            # if filename is missing, save back to file from
            # which we loaded configuration
            if hasattr(self, 'filename'):
                dest = self.filename
            else:
                raise ValueError("No 'dest' and Atoms has no stored filename")

        from quippy.io import AtomsWriters
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

    def select(self, mask=None, list=None, orig_index=True):
        """Return a new :class:`Atoms` containing a subset of the atoms in this Atoms object

        One of either `mask` or `list` should be present.  If `mask`
        is given it should be a rank one array of length `self.n`. In
        this case atoms corresponding to true values in `mask` will be
        included in the result.  If `list` is present it should be an
        arry of list containing atom indices to include in the result.

        If `orig_index` is True (default), the new object will contain
        an ``orig_index`` property mapping the indices of the new atoms
        back to the original larger Atoms object.
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
        """
        Return a copy of this :class:`Atoms` object
        """
        
        other = self.__class__(n=self.n, lattice=self.lattice,
                               properties=self.properties, params=self.params)

        # copy any normal (not Fortran) attributes
        for k, v in self.__dict__.iteritems():
            if not k.startswith('_') and k not in other.__dict__:
                other.__dict__[k] = v

        # from _atoms.Atoms
        other.cutoff = self.cutoff
        other.cutoff_skin = self.cutoff_skin
        other.nneightol = self.nneightol

        # from ase.Atoms
        other.constraints = copy.deepcopy(self.constraints)
        other.adsorbate_info = copy.deepcopy(self.adsorbate_info)
        return other

    def copy_from(self, other):
        """Replace contents of this Atoms object with data from `other`."""

        self.__class__.__del__(self)
        if isinstance(other, _atoms.Atoms):
            _atoms.Atoms.__init__(self, n=other.n, lattice=other.lattice,
                                  properties=other.properties, params=other.params)

            self.cutoff = other.cutoff
            self.cutoff_skin = other.cutoff_skin
            self.nneightol = other.nneightol
            
        elif isinstance(other, ase.Atoms):
            _atoms.Atoms.__init__(self, n=0, lattice=np.eye(3))
            ase.Atoms.__init__(self, other)

            # copy params/info dicts
            if hasattr(other, 'params'):
                self.params.update(other.params)
            if hasattr(other, 'info'):
                self.params.update(other.info)
                if 'nneightol' in other.info:
                    self.nneightol = other.info['nneightol']
                if 'cutoff' in other.info:
                    self.set_cutoff(other.info['cutoff'],
                                    other.info.get('cutoff_break'))
                if 'cutoff_factor' in other.info:
                    self.set_cutoff_factor(other.info['cutoff_factor'],
                                           other.info.get('cutoff_factor_break'))

            # create extra properties for any non-standard arrays
            standard_ase_arrays = ['positions', 'numbers', 'masses', 'charges',
                                   'momenta', 'tags', 'magmoms' ]
                
            for ase_name, value in other.arrays.iteritems():
                quippy_name = self.name_map.get(ase_name, ase_name)
                if ase_name not in standard_ase_arrays:
                    self.add_property(quippy_name, np.transpose(value))

            self.constraints = copy.deepcopy(other.constraints)
            self.adsorbate_info = copy.deepcopy(other.adsorbate_info)
            
        else:
            raise TypeError('can only copy from instances of quippy.Atoms or ase.Atoms')
        
        # copy any normal (not Fortran) attributes
        for k, v in other.__dict__.iteritems():
            if not k.startswith('_') and k not in self.__dict__:
                self.__dict__[k] = v
        

    def read_from(self, source, **readargs):
        """Replace contents of this Atoms object with Atoms read from `source`"""
        try:
            self.copy_from(source)
        except TypeError:
            tmp = Atoms.read(source, **readargs)
            self.shallow_copy_from(tmp)
            # tmp goes out of scope here, but reference counting
            # prevents it from being free'd.

    def __getattr__(self, name):
        #print 'getattr', name
        #if name in self.properties:
        if name == '_fpointer':
            raise AttributeError('Atoms object not initialised!')
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
        if (get_fortran_indexing() and (i < 1 or i > self.n)) or \
            (not get_fortran_indexing() and (i < 0 or i > self.n-1)):
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

        Here are some examples::

            a = Atoms(n=10, lattice=10.0*fidentity(3))

            a.add_property('mark', 1)                  # Scalar integer
            a.add_property('bool', False)              # Scalar logical
            a.add_property('local_energy', 0.0)        # Scalar real
            a.add_property('force', 0.0, n_cols=3)     # Vector real
            a.add_property('label', '')                # Scalar string

            a.add_property('count', [1,2,3,4,5,6,7,8,9,10])  # From list
            a.add_property('norm_pos', a.pos.norm())         # From 1D array
            a.add_property('pos', new_pos)                   # Overwrite positions with array new_pos
                                                             # which should have shape (3,10)
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
        return self.write('string')

    def __setstate__(self, state):
        self.read_from(state, format='string')

    def __reduce__(self):
        return (Atoms, (), self.__getstate__(), None, None)

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
FortranDerivedTypes['type(connection)'] = Connection
