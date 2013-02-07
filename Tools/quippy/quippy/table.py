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

import weakref
import numpy as np
import quippy
from quippy import _table
from quippy._table import *
from quippy.farray import FortranArray

__doc__ = _table.__doc__
__all__ = _table.__all__

class Table(_table.Table):
    __doc__ = _table.Table.__doc__
    _cmp_skip_fields = ['max_length', 'increment']

    def __repr__(self):
        return ('Table(n=%d,intsize=%d,realsize=%d,strsize=%d,logicalsize=%d)' %
                (self.n, self.intsize, self.realsize, self.strsize, self.logicalsize))

    def __str__(self):
        return repr(self)

    def copy(self):
        t = Table(self.intsize, self.realsize, self.strsize, self.logicalsize, self.n)
        t.append(blank_rows=self.n)
        if self.intsize != 0: t.int[...] = self.int[...]
        if self.realsize != 0: t.real[...] = self.real[...]
        if self.strsize != 0: t.str[...] = self.str[...]
        if self.logicalsize != 0: t.logical[...] = self.logical[...]
        return t

    def _get_array_shape(self, name):
        if name in ('int','real','logical'):
            return (slice(None),slice(1,self.n))
        elif name == 'str':
            return (slice(None),slice(None),slice(1,self.n))
        else:
            return None

    @classmethod
    def from_atom_indices(cls, atoms, mask=None, list=None, force_fortran_indexing=True):
        """
        Construct a new Table containing atomic indices from a list or mask

        The new table will include 4 integer columns, for indices plus
        shifts, and is suitable for passing to bfs_step,
        construct_hysteretic_region, etc.

        If force_fortran_indexing is True (the default), all atom indices
        are converted to Fortran 1-based indexing.
        """
        orig_fortran_indexing = atoms.fortran_indexing
        atoms.fortran_indexing = force_fortran_indexing

        try:
            if mask is None and list is None:
                raise ValueError('Either mask or list must be present.')

            if list is not None and force_fortran_indexing and not orig_fortran_indexing:
                # we were given 0-based atom indices, convert to 1-based indices
                list = np.array(list)+1

            if mask is not None:
                if len(mask) != len(atoms):
                    raise ValueError('size of mask must equal number of atoms')
                mask = mask.astype(bool)
                # use indexing style given by force_fortran_indexing
                list = atoms.indices[mask].nonzero()[0]

            self = cls(4,0,0,0, fortran_indexing=force_fortran_indexing)
            self.append(blank_rows=len(list))
            self.int[:,:] = 0

            if self.fortran_indexing:
                first_column = 1
            else:
                first_column = 0
            self.int[first_column,:] = list
            
        finally:
            atoms.fortran_indexing = orig_fortran_indexing

        self.atoms = weakref.ref(atoms)
        return self


    @classmethod
    def from_atom_list(cls, atoms, list, force_fortran_indexing=True):
        """
        Construct a new Table from a list of atom indices

        See Table.from_atom_indices for more details.
        """
        return cls.from_atom_indices(atoms, list=list,
                                     force_fortran_indexing=force_fortran_indexing)


    @classmethod
    def from_atom_mask(cls, atoms, mask, force_fortran_indexing=True):
        """
        Construct a new Table from an atom mask

        See Table.from_atom_indices for more details.
        """
        return cls.from_atom_indices(atoms, mask=mask,
                                     force_fortran_indexing=force_fortran_indexing)


    def to_atom_list(self, atoms=None):
        """
        Return list of atom indices that this Table represents.

        Indices returns are 0-based or 1-based depending on value of
        atoms.fortran_indexing.

        This assumes that indexing scheme of Table contents is consistent
        with self.fortran_indexing (this is true if Table was constructed
        by Table.from_atom_indices).
  
        If `atoms` is not present, the Atoms object passed to
        Table.from_atom_indices is used, or an exception is raised if this
        Table was not created in that way.
        """
        if atoms is None:
            if not hasattr(self, 'atoms'):
                raise AttributeError('Table missing "atoms" attribute, probably'+
                                     ' not created by Table.from_atom_indices()')
            atoms = self.atoms()
            if atoms is None:
                raise ValueError('weakref to Table.atoms has expired')
                
        if self.fortran_indexing:
            first_column = 1
        else:
            first_column = 0
        indices = self.int[first_column,:].copy()
        if self.fortran_indexing and not atoms.fortran_indexing:
            indices -= 1
        elif not self.fortran_indexing and atoms.fortran_indexing:
            indices += 1
        return list(indices)


    def to_atom_mask(self, atoms=None):
        """
        Return mask for atom indices that this Table represents

        Result is either an array of an FortranArray, depending on
        value of atoms.fortran_indexing.

        This assumes that indexing scheme of Table contents is consistent
        with self.fortran_indexing (this is true if Table was constructed
        by Table.from_atom_indices).

        If `atoms` is not present, the Atoms object passed to
        Table.from_atom_indices is used, or an exception is raised if this
        Table was not created in that way.        
        """

        mask = np.zeros(len(atoms), dtype=bool)
        if atoms.fortran_indexing:
            mask = mask.view(FortranArray)
        mask[self.to_atom_list(atoms)] = True
        return mask
        

from quippy import FortranDerivedTypes
FortranDerivedTypes['type(table)'] = Table
