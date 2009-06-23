class DictMixin:
    '''Mixin defining all dictionary methods for classes that already have
       a minimum dictionary interface including getitem, setitem, delitem,
       and keys 

       from http://code.activestate.com/recipes/117236/
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
    def get(self, key, default):
        if key in self:
            return self[key]
        return default
    def __repr__(self):
        return repr(dict(self.items()))

def MakeFullDict(tgt):
    'Extends the dictionary interface for existing classes'
    tgt.__bases__ = tuple(list(tgt.__bases__) + [DictMixin])
