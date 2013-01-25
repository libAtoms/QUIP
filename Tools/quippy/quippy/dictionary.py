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

from quippy import _dictionary
from quippy._dictionary import *
from quippy.oo_fortran import update_doc_string
from quippy.dictmixin import DictMixin, ParamReaderMixin
from quippy.farray import *

__all__ = _dictionary.__all__

class Dictionary(DictMixin, ParamReaderMixin, _dictionary.Dictionary):
    __doc__ = update_doc_string(_dictionary.Dictionary.__doc__, """
    
    The quippy Python :class:`Dictionary` class is designed to behave
    as much as possible like a true Python dictionary, but since it is
    implemented in Fortran it can only store a restricted range of
    data types.  Keys much be strings and values must be one of the
    types above.

    Trying to store any other type of data will raise a :exc:`ValueError`.

    For Atoms' :attr:`~quippy.atoms.Atoms.params` entries, there are
    further restrictions imposed by the implementation of the XYZ and
    NetCDF I/O routines. The only types of data that can be stored here
    are:

    - Integer
    - Real
    - String
    - Integer 3-vector
    - Real 3-vector
    - Integer 3 x 3 matrix
    - Real 3 x 3 matrix

    A :class:`Dictionary` can be created from a standard Python dictionary,
    and easily converted back::

        >>> py_dict = {'a':1, 'b':2}
        >>> fortran_dict = Dictionary(py_dict)
        >>> py_dict == dict(fortran_dict)
        True

    It also supports all the standard dictionary operations and methods::

        >>> fortran_dict['c'] = 3
	>>> fortran_dict.keys()
	['a', 'b', 'c']

    An additional feature of the quippy :class:`Dictionary` is that it
    can read and write itself to a string in the format used within XYZ
    files::

        >>> str(fortran_dict)
	'a=1 b=2 c=3'
	>>> d2 = Dictionary('a=1 b=2 c=3')
	>>> d2.keys(), d2.values()      
	(['a', 'b', 'c'], [1, 2, 3])
    """)

    _interfaces = _dictionary.Dictionary._interfaces
    _interfaces['set_value'] = [ k for k in _dictionary.Dictionary._interfaces['set_value'] if k[0] != 'set_value_s_a' ]

    _scalar_types = (T_INTEGER, T_REAL, T_COMPLEX, T_LOGICAL, T_CHAR, T_DICT)

    _array_types  = (T_INTEGER_A, T_REAL_A, T_COMPLEX_A, T_CHAR_A,
                     T_LOGICAL_A, T_INTEGER_A2, T_REAL_A2)

    def __init__(self, D=None, *args, **kwargs):
        _dictionary.Dictionary.__init__(self, *args, **kwargs)
        self._cache = {}
        self.key_cache_invalid = 1
        self._keys = []
        self._keys_lower = []
        if D is not None:
            self.read(D) # copy from D

    def keys(self):
        if self.key_cache_invalid or len(self._keys) != self.n: # HACK to solve shallow copy bug
            self._keys = [self.get_key(i).strip() for i in frange(self.n)]
            self._keys_lower = [k.lower() for k in self._keys]
            self.key_cache_invalid = 0
        return self._keys

    def has_key(self, key):
        k = self.keys() # ensure _keys_lower is up-to-date
        return key.lower() in self._keys_lower

    def get_value(self, k):
        "Return a _copy_ of a value stored in Dictionary"
        if not k in self:
            raise KeyError('Key "%s" not found ' % k)

        t, s, s2 = self.get_type_and_size(k)

        if t == T_NONE:
            raise ValueError('Key %s has no associated value' % k)
        elif t == T_INTEGER:
            v,p = self._get_value_i(k)
        elif t == T_REAL:
            v,p = self._get_value_r(k)
        elif t == T_COMPLEX:
            v,p = self._get_value_c(k)
        elif t == T_CHAR:
            v,p = self._get_value_s(k)
            v = v.strip()
        elif t == T_LOGICAL:
            v,p = self._get_value_l(k)
            v = bool(v)
        elif t == T_INTEGER_A:
            v,p = self._get_value_i_a(k,s)
        elif t == T_REAL_A:
            v,p = self._get_value_r_a(k,s)
        elif t == T_COMPLEX_A:
            v,p = self._get_value_c_a(k,s)
        elif t == T_CHAR_A:
            v,p = self._get_value_s_a2(k,s2[1], s2[2])
            v = v[...,1]  # Last index is length of string, here fixed to 1
            v.strides = (1, v.shape[0])  # Column-major storage
        elif t == T_LOGICAL_A:
            v,p = self._get_value_l_a(k,s)
            v = farray(v, dtype=bool)
        elif t == T_INTEGER_A2:
            v,p = self._get_value_i_a2(k, s2[1], s2[2])
        elif t == T_REAL_A2:
            v,p = self._get_value_r_a2(k, s2[1], s2[2])
        elif t == T_DICT:
            v,p = self._get_value_dict(k)
        else:
            raise ValueError('Unsupported dictionary entry type %d' % t)
        return v

    def get_array(self, key):
        "Return a _reference_ to an array stored in this Dictionary"""

        import _quippy, arraydata
        if key in self and self.get_type_and_size(key)[0] in Dictionary._array_types:
            a = arraydata.get_array(self._fpointer, _quippy.qp_dictionary_get_array, key)
            if self.fortran_indexing:
                a = FortranArray(a, parent=self)
            return a
        else:
            raise KeyError('Key "%s" does not correspond to an array entry' % key)

    def get_type(self, key):
        "Return an integer code for the type of the value associated with a key"

        import _quippy, arraydata
        if key in self:
            return self.get_type_and_size(key)[0]
        else:
            raise KeyError('Key "%s" not found' % key)


    def is_scalar(self, key):
        if key in self:
            return self.get_type_and_size(key)[0] in Dictionary._scalar_types
        else:
            raise KeyError('Key "%s" not found')

    def is_array(self, key):
        if key in self:
            return self.get_type_and_size(key)[0] in Dictionary._array_types
        else:
            raise KeyError('Key "%s" not found')


    def __getitem__(self, k):
        k = k.lower()

        if self.cache_invalid:
            self._cache = {}
            self.cache_invalid = 0

        try:
            v = self._cache[k]
            if v is None: raise KeyError
            return v

        except KeyError:
            if not k in self:
                raise KeyError('Key "%s" not found ' % k)

            t = self.get_type_and_size(k)[0]

            if t == T_NONE:
                raise ValueError('Key %s has no associated value' % k)

            elif t in Dictionary._scalar_types:
                self._cache[k] = None
                return self.get_value(k)

            elif t in Dictionary._array_types:
                v = self.get_array(k)
                self._cache[k] = v
                return v

            else:
                raise ValueError('Unsupported dictionary entry type %d' % t)


    def __setitem__(self, k, v):
        try:
            self.set_value(k, v)
        except TypeError:
            self.set_value(k,s2a(v,pad=None))

    def __delitem__(self, k):
        if not k in self:
            raise KeyError('Key %s not found in Dictionary' % k)
        self.remove_value(k)

    def __repr__(self):
        return ParamReaderMixin.__repr__(self)

    def __eq__(self, other):
        import logging
        if sorted(self.keys()) != sorted(other.keys()):
            logging.debug('keys mismatch: %s != %s' % (sorted(self.keys()), sorted(other.keys())))
            return False

        for key in self:
            v1, v2 = self[key], other[key]
            if isinstance(v1, np.ndarray) and isinstance(v2, np.ndarray):
                if v1.size == 0 and v2.size == 0:
                    continue
                elif v1.dtype.kind != 'f':
                    if (v1 != v2).any():
                        logging.debug('mismatch key=%s v1=%s v2=%s' % (key, v1, v2))
                        return False
                else:
                    if abs(v1 - v2).max() > self._cmp_tol:
                        logging.debug('mismatch key=%s v1=%s v2=%s' % (key, v1, v2))
                        return False
            else:
                if v1 != v2:
                    logging.debug('mismatch key=%s v1=%s v2=%s' % (key, v1, v2))
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return ParamReaderMixin.__str__(self)

    def copy(self):
        return Dictionary(self)

    def subset(self, keys, out=None, case_sensitive=None, out_no_initialise=None):
        if out is None: out = Dictionary()
        _dictionary.Dictionary.subset(self, keys, out, case_sensitive, out_no_initialise)
        return out

from quippy import FortranDerivedTypes
FortranDerivedTypes['type(dictionary)'] = Dictionary
