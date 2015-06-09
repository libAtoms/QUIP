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
    def from_atom_indices(cls, atoms, mask=None, list=None):
        """
        Construct a new Table containing atomic indices from a list or mask

        The new table will include 4 integer columns, for indices plus
        shifts, and is suitable for passing to bfs_step,
        construct_hysteretic_region, etc.
        """
        if mask is None and list is None:
            raise ValueError('Either mask or list must be present.')

        if mask is not None:
            if len(mask) != len(atoms):
                raise ValueError('size of mask must equal number of atoms')
            mask = mask.astype(bool)
            list = atoms.indices[mask]

        self = cls(4,0,0,0)
        self.append(blank_rows=len(list))
        self.int[:,:] = 0

        self.int[0,:] = list

        self.atoms = weakref.ref(atoms)
        return self


    @classmethod
    def from_atom_list(cls, atoms, list):
        """
        Construct a new Table from a list of atom indices

        See Table.from_atom_indices for more details.
        """
        return cls.from_atom_indices(atoms, list=list)


    @classmethod
    def from_atom_mask(cls, atoms, mask):
        """
        Construct a new Table from an atom mask

        See Table.from_atom_indices for more details.
        """
        return cls.from_atom_indices(atoms, mask=mask)


    def to_atom_list(self, atoms=None):
        """
        Return list of atom indices that this Table represents.

        Indices returned are 0-based.

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
                
        indices = self.int[0, :].copy()
        return list(indices)


    def to_atom_mask(self, atoms=None):
        """
        Return mask for atom indices that this Table represents

        Result is either a numpy array.

        If `atoms` is not present, the Atoms object passed to
        Table.from_atom_indices is used, or an exception is raised if this
        Table was not created in that way.        
        """

        mask = np.zeros(len(atoms), dtype=bool)
        mask[self.to_atom_list(atoms)] = True
        return mask
        

from quippy import FortranDerivedTypes
FortranDerivedTypes['type(table)'] = Table
