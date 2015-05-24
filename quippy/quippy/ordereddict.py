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

class OrderedDict(dict):
    """Subclass of dictionary that keeps track of the order in which keys were
    added, and iterates over them in that order"""
    def __init__(self,D=None):
        self._keys = []
        dict.__init__(self)
        if D is not None:
            self.update(D)

    @staticmethod
    def frompairs(seq):
        d = OrderedDict()
        for (k,v) in seq:
            d[k] = v
        return d

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        dict.__setitem__(self, key, item)
        if key not in self._keys: self._keys.append(key)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, dict.__repr__(self))

    def clear(self):
        dict.clear(self)
        self._keys = []

    def copy(self):
        D = OrderedDict(dict.copy(self))
        D._keys = self._keys[:]
        return D

    def __copy__(self):
        return self.copy()

    def items(self):
        return zip(self._keys, self.values())


    def iteritems(self):
        for pair in zip(self._keys, self.values()):
            yield pair

    def __iter__(self):
        for key in self._keys:
            yield key

    def iterkeys(self):
        for key in self._keys:
            yield key

    def itervalues(self):
        for key in self._keys:
            yield self[key]

    def keys(self):
        return self._keys

    def popitem(self):
        try:
            key = self._keys[-1]
        except IndexError:
            raise KeyError('dictionary is empty')

        val = self[key]
        del self[key]

        return (key, val)

    def setdefault(self, key, failobj = None):
        res = dict.setdefault(self, key, failobj)
        if key not in self._keys: self._keys.append(key)
        return res

    def update(self, D):
        dict.update(self,D)
        for key in D.keys():
            if key not in self._keys: self._keys.append(key)

    def values(self):
        return map(self.get, self._keys)


    def rename(self, old_key, new_key):
        if new_key == old_key:
            return

        if new_key in self.keys():
            raise ValueError('New key already exists: %r' % new_key)

        value = self[old_key]
        old_idx = self._keys.index(old_key)
        self._keys[old_idx] = new_key

        dict.__delitem__(self, old_key)
        dict.__setitem__(self, new_key, value)
