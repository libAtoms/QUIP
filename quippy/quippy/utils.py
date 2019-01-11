"""
Utilities for quippy

"""

import _quippy
import quippy
import numpy as np
import f90wrap.runtime


def get_dict_arrays(fdict: quippy.dictionary_module.Dictionary):
    """Takes the arrays from a quippy dictionary. Copies.

    Probably fails if there are non-array elements in the dictionary"""

    if not isinstance(fdict, quippy.dictionary_module.Dictionary):
        raise TypeError('fdict argument is not a quippy.dictionary_module.Dictionary')

    arrays = {}
    for i in range(1, fdict.n + 1):
        key, error = fdict.get_key(i)
        key = key.strip().decode('ascii')
        print(key)
        # fixme: fails for non_array elements. Make universal: compatible with array or scalar content in dictionary
        try:    # this is an unsufficient temporary fix
            value = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                              fdict._handle, _quippy.f90wrap_dictionary__array__, key)
            arrays[key] = value.copy()
        except ValueError:
            pass

    return arrays

