from __future__ import print_function, absolute_import, division

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

# Original DictMixin version from
# http://code.activestate.com/recipes/117236/
# Copyright (c) Raymond Hettinger 2002
# Licensed under the PSF License

import numpy as np

class DictMixin(object):
    '''Mixin defining all dictionary methods for classes that already have
       a minimum dictionary interface including getitem, setitem, delitem,
       and keys
       '''

    # first level definitions should be implemented by the sub-class
    def __getitem__(self, key):
        raise NotImplementedError
    def __setitem__(self, key, value):
        raise NotImplementedError
    def __delitem__(self, key):
        raise NotImplementedError
    def keys(self):
        raise NotImplementedError

    # second level definitions which assume only getitem and keys
    def has_key(self, key):
        return key in self.keys()
    def __iter__(self):
        for k in self.keys():
            yield k
    def __len__(self):
        return len(self.keys())

    # third level uses second level instead of first
    def __contains__(self, key):
        return self.has_key(key)
    def iteritems(self):
        for k in self:
            yield (k, self[k])

    # fourth level uses second and third levels instead of first
    def iterkeys(self):
        return self.__iter__()
    def itervalues(self):
        for _, v in self.iteritems():
            yield v
    def values(self):
        return list(self.itervalues())
    def items(self):
        return list(self.iteritems())
    def clear(self):
        for key in self.keys():
            del self[key]
    def setdefault(self, key, default):
        if key not in self:
            self[key] = default
            return default
        return self[key]
    def popitem(self):
        key = self.keys()[0]
        value = self[key]
        del self[key]
        return (key, value)
    def update(self, other):
        for key in other.keys():
            self[key] = other[key]
    def get(self, key, default=None):
        if key in self:
            return self[key]
        return default
    def __repr__(self):
        return repr(dict(self.items()))

def MakeFullDict(tgt):
    'Extends the dictionary interface for existing classes'
    tgt.__bases__ = tuple(list(tgt.__bases__) + [DictMixin])


from farray import *
import string, re

class ParamReaderMixin(object):
    """Mixin which adds parse(), read(), write() and asstring() to a dictionary"""

    def __copy__(self):
        return self.copy()

    def parse(self, s):
        key_quoted_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*["\{\}]([^"\{\}]+)["\{\}e+-]\s*')
        key_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*([-0-9A-Za-z_.:\[\]()e+-]+)\s*')
        key_re = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*')

        s = s.strip()

        while 1:
            # Match quoted string first, then fall through to plain key=value
            m = key_quoted_value.match(s)
            if m is None:
                m = key_value.match(s)
                if m is not None:
                    s = key_value.sub('', s, 1)
                else:
                    # Just a key with no value
                    m = key_re.match(s)
                    if m is not None:
                        s = key_re.sub('', s, 1)
            else:
                s = key_quoted_value.sub('', s, 1)

            if m is None: break # No more matches

            key = m.group(1)
            try:
                value = m.group(2)
            except IndexError:
                # default value is None
                value = None

            # Try to convert to (list of) floats, ints
            try:
                numvalue = []
                for x in string.split(value):
                    if x.find('.') == -1:
                        numvalue.append(int(float(x)))
                    else:
                        numvalue.append(float(x))
                if len(numvalue) == 1:
                    numvalue = numvalue[0] # Only one number
                elif len(numvalue) == 3:
                    numvalue = farray(numvalue) # 3-vector
                elif len(numvalue) == 9:
                    # 3x3 matrix, fortran ordering
                    numvalue = farray(numvalue).reshape((3,3), order='F')
                else:
                    raise ValueError
                value = numvalue
            except (AttributeError, ValueError):
                pass

            # Parse boolean values, e.g 'T' -> True, 'F' -> False, 'T T F' -> [True, True, False]
            if isinstance(value, str):
                str_to_bool  = {'T':True, 'F':False}

                if len(value.split()) > 1:
                    if all([x in str_to_bool.keys() for x in value.split() ]):
                        value = [str_to_bool[x] for x in value.split()]
                elif value in str_to_bool:
                    value = str_to_bool[value]

            self[key] = value

    def read(self, f):
        if isinstance(f, str):
            self.parse(f)
        elif hasattr(f, 'keys') and hasattr(f, '__getitem__'):
            self.update(f)
        else:
            try:
                for line in f:
                    self.parse(line)
            except TypeError:
                raise TypeError("Don't know how to read from object - "+str(f))


    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__, str(self))

    def __str__(self):
        return self.asstring()

    def asstring(self, sep=' '):
        if len(self) == 0: return ''

        type_val_map = {(bool,True): None,
                        (bool,False): 'F',
                        (np.bool_,True): None,
                        (np.bool_,False): 'F'}

        s = ''
        for key in self.keys():
            val = self[key]

            if hasattr(val, '__iter__'):
                val = ' '.join( str(type_val_map.get((type(x),x), x))
                        for x in farray(val).reshape(farray(val).size, order='F'))
                val.replace('[', '')
                val.replace(']', '')
            else:
                val = type_val_map.get((type(val),val), val)

            if val is None:
                s = s + '%s%s' % (key, sep)
            elif type(val) == type('') and ' ' in val:
                s = s + '%s="%s"%s' % (key, val, sep)
            else:
                s = s + '%s=%s%s' % (key, str(val), sep)

        return s.strip()

    def write(self, f):
        f.write(self.asstring(sep='\n'))

from quippy.ordereddict import OrderedDict

class PuPyDictionary(OrderedDict, ParamReaderMixin):
    """Subclass of OrderedDict for reading key/value pairs from strings or files.
       The original order of items is maintained. Values that looks like floats or ints
       or lists of floats or ints are automatically converted on reading."""

    def __init__(self, source=None):
        OrderedDict.__init__(self)
        if source is not None:
            self.read(source)

    def __repr__(self):
        return ParamReaderMixin.__repr__(self)

    def __str__(self):
        return ParamReaderMixin.__str__(self)

    def copy(self):
        return PuPyDictionary(OrderedDict.copy(self))
