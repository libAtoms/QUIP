# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright ST John 2017
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

from quippy import _descriptors #as _descriptors #TODO rename fortran-wrapped module to _descriptors
from quippy.oo_fortran import update_doc_string
from quippy.farray import fzeros
from ase.atoms import Atoms as ASEAtoms
from quippy.atoms import Atoms
import numpy as np

__all__ = ['Descriptor']

class DescriptorCalcResult:
    """
    Results of a descriptor calculation with gradients.
    """
    def __init__(self, data, grad, index):
        self.data = data
        self.grad = grad
        self.index = index
        #self.i_desc, self.i_atom, self.i_coord = index

from functools import wraps
def convert_atoms_types_iterable_method(method):
    """
    Decorator to transparently convert ASEAtoms objects into quippy Atoms, and
    to transparently iterate over a list of Atoms objects...
    """
    @wraps(method)
    def wrapper(self, at, *args, **kw):
        if isinstance(at, Atoms):
            return method(self, at, *args, **kw)
        elif isinstance(at, ASEAtoms):
            return method(self, Atoms(at), *args, **kw)
        else:
            return [wrapper(self, atelement, *args, **kw) for atelement in at]
    return wrapper

class Descriptor(_descriptors.Descriptor):
    __doc__ = update_doc_string(_descriptors.Descriptor.__doc__, """
    Pythonic wrapper for GAP descriptor module""",
            signature='Descriptor(args_str)')

    def __init__(self, args_str):
        """
        Initialises Descriptor object and calculate number of dimensions and
        permutations.
        """
        _descriptors.Descriptor.__init__(self, args_str)
        self._n_dim = self.dimensions()
        self._n_perm = self.n_permutations()

    n_dim = property(lambda self: self._n_dim)
    n_perm = property(lambda self: self._n_perm)

    def __len__(self):
        return self.n_dim

    def permutations(self):
        """
        Returns array containing all valid permutations of this descriptor.
        """
        perm = _descriptors.Descriptor.permutations(self, self.n_dim, self.n_perm)
        return np.array(perm).T

    @convert_atoms_types_iterable_method
    def count(self, at):
        """
        Returns how many descriptors of this type are found in the Atoms
        object.
        """
        return self.descriptor_sizes(at)[0]

    @convert_atoms_types_iterable_method
    def calc(self, at):
        """
        Calculates all descriptors of this type in the Atoms object, and
        returns the array of descriptor values. Does not compute gradients; use
        calc_grad() for that.
        """
        n_desc, n_cross = self.descriptor_sizes(at)
        data = fzeros((self.n_dim, n_desc))
        _descriptors.Descriptor.calc(self, at, data)
        return np.array(data).T

    @convert_atoms_types_iterable_method
    def calc_grad(self, at):
        """
        Calculates all descriptors of this type in the Atoms object, and
        returns a DescriptorCalcResult object that contains the descriptor
        values, the gradients, and the indices to the involved atoms.
        """
        # what do we want to get out of it?
        n_desc, n_cross = self.descriptor_sizes(at)
        data = fzeros((self.n_dim, n_desc))
        # n_cross is number of cross-terms, proportional to n_desc
        data_grad = fzeros((self.n_dim, 3*n_cross))
        data_index = fzeros((3, 3*n_cross), 'i')
        _descriptors.Descriptor.calc(self, at, data, data_grad, data_index)
        # TODO: convert to numpy arrays and expose in suitable fashion
        return DescriptorCalcResult(data, data_grad, data_index)

