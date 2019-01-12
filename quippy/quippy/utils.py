# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
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


"""
Utilities for quippy

"""

import _quippy
import quippy
import f90wrap.runtime
from copy import deepcopy as cp

__all__ = ['get_dict_arrays']


def get_dict_arrays(fdict):
    """Takes the arrays from a quippy dictionary. Copies.

    Probably fails if there are non-array elements in the dictionary"""

    if not isinstance(fdict, quippy.dictionary_module.Dictionary):
        raise TypeError('fdict argument is not a quippy.dictionary_module.Dictionary')

    arrays = {}
    for i in range(1, fdict.n + 1):
        key, error = fdict.get_key(i)
        key = key.strip().decode('ascii')
        # fixme: fails for non_array elements. Make universal: compatible with array or scalar content in dictionary
        try:    # this is an unsufficient temporary fix
            value = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                              fdict._handle, _quippy.f90wrap_dictionary__array__, key)
            arrays[key] = value.copy()
        except ValueError:
            value = fdict.get_value(key)
            try:
                # normally it is an tuple, because the error arf from fortran is converted to output
                arrays[key] = cp(value[0])
            except TypeError:
                arrays[key] = cp(value)

    return arrays

