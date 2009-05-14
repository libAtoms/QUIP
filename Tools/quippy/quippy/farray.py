"""Contains FortranArray class and utility functions for handling one-based
array indexing."""

import numpy

def frange(min,max=None,step=1):
    """Return an iterator for integers from min to max inclusive, increasing by
    step each time."""
    if max is None:
        return xrange(1,min+1,step)
    else:
        return xrange(min,max+1,step)

def fzeros(shape,dtype=float):
    return FortranArray(numpy.zeros(shape,dtype,order='F'))

def farray(seq):
    return FortranArray(numpy.array(seq,order='F'))

class FortranArray(numpy.ndarray):
    """Subclass of ndarray which uses Fortran-style one-based
    indexing. The first element is numbered one rather than zero and
    trying to access element zero will raise an IndexError
    exception. Negative indices are unchanged; -1 still refers to the
    highest index along a dimension. Slices are also Fortran-style,
    i.e. inclusive on the second element so 1:2 includes the one-based
    elements 1 and 2, equivalent to a C-style slice of 0:2. The self.flat
    iterator is still indexed from zero."""

    def additerators(self):
        self.rows = self.row_iter()
        self.cols = self.__iter__()

    def __new__(self, input_array=None, doc=None):
        """a = FortranArray(input_array=None, doc=None)

        Construct a FortanArray from input_array, optionally setting
        docstring to doc."""

	self = numpy.asarray(input_array).view(FortranArray)
        
        if doc is not None:
            self.__doc__ = doc
        self.additerators()
	return self

    @staticmethod
    def map_int(idx):
        if idx > 0:
            return idx-1
        elif idx == 0:
            raise IndexError('index 0 not permitted - FortranArrays are one-based')
        else:
            return idx

    @staticmethod
    def mapindices(indx):
        """Transform from Fortran-style one-based to C-style zero based
        indices. indx can be any object that can be used to subscript
        an ndarray - scalar integer or slice, sequence of integers and
        slices or an ndarray of integer or boolean kind."""

        one_dim = False
        if not hasattr(indx, '__iter__'):
            one_dim = True
            indx = (indx,)

        islist = isinstance(indx, list)

        if isinstance(indx, numpy.ndarray):
            indx = (indx,)

        nindx = []
        for idx in indx:
            if isinstance(idx, int) or isinstance(idx, numpy.integer):
#                print 'mapping int %d to %d' % (idx, FortranArray.map_int(idx))
                nindx.append(FortranArray.map_int(idx))
            elif isinstance(idx, slice):
                rslice_start = None
                rslice_stop = None
                rslice_step = None
                if idx.start is not None:
                    rslice_start = FortranArray.map_int(idx.start)
                if idx.stop is not None:
                    rslice_stop = idx.stop
                if idx.step is not None:
                    rslice_step = idx.step
                rslice = slice(rslice_start, rslice_stop, rslice_step)
#                print 'mapping slice %s to %s' % (idx, rslice)
                nindx.append( rslice )
            elif isinstance(idx, numpy.ndarray):
                if idx.dtype.kind == 'i':
                    if (idx == 0).any(): raise ValueError('Advanced slicing array must not contain index 0')
                    nindx.append(numpy.where(idx > 0, idx-1, idx))
                elif idx.dtype.kind == 'b':
                    nindx.append(idx)
                else:
                    raise ValueError('Advanced slicing array must be integer or boolean')
            elif idx is Ellipsis or idx is None:
                nindx.append(idx)
            else:
                raise ValueError('Unknown index object %r' % (idx,))

        if one_dim:
            return nindx[0]
        else:
            if islist:
                return nindx
            else:
                return tuple(nindx)

    def __getitem__(self, indx):
        "Overloaded __getitem__ which accepts one-based indices."
	indx = FortranArray.mapindices(indx)
	obj = numpy.ndarray.__getitem__(self, indx) 
	if (isinstance(obj, numpy.ndarray) and obj.dtype.isbuiltin):
            fa = obj.view(FortranArray)
            fa.additerators()
	    return fa
	return obj

    def __setitem__(self, indx, value):
        "Overloaded __setitem__ which accepts one-based indices."
        if not isinstance(indx, slice) and not (hasattr(indx, '__iter__') and any([isinstance(x,slice) for x in indx])):
            # if indx contains a slice then __getitem__ will be called and mapping will be done twice
            indx = FortranArray.mapindices(indx)
	numpy.ndarray.__setitem__(self, indx, value)


    def __getslice__(self, i, j):
        i = FortranArray.map_int(i)
        j = FortranArray.map_int(j)
        obj = numpy.ndarray.__getslice__(self, i, j)
	if (isinstance(obj, numpy.ndarray) and obj.dtype.isbuiltin):
            fa = obj.view(FortranArray)
            fa.additerators()
	    return fa
	return obj

    def __setslice__(self, i, j, value):
        i = FortranArray.map_int(i)
        j = FortranArray.map_int(j)
        numpy.ndarray.__setslice__(self, i, j, value)

    def nonzero(self):
        """a.nonzero()

        Return the one-based indices of the elements of a which are not zero."""
	return tuple(a + 1 for a in numpy.ndarray.nonzero(self))
	
    def count(self):
        """a.count()

        Number of nonzero elemnts of this FortranArray. Equivalent to
        len(self[self.nonzero()])."""
	return len(self[self.nonzero()])
    
    def argmin(self, axis=None, out=None):
        """a.argmin(axis=None, out=None)
    
        Return one-based indices of the minimum values along the given
        axis of `a`.
        
        Refer to `numpy.ndarray.argmax` for detailed documentation."""

	return numpy.ndarray.argmin(self,axis,out) + 1

    def argmax(self, axis=None, out=None):
        """a.argmax(axis=None, out=None)
    
        Return one-based indices of the maximum values along the given
        axis of `a`.
        
        Refer to `numpy.ndarray.argmax` for detailed documentation."""

	return numpy.ndarray.argmax(self,axis,out) + 1	

    def argsort(self, axis=-1, kind='quicksort', order=None):
        """a.argsort(axis=-1, kind='quicksort', order=None)
    
        Returns the indices that would sort this array.
        
        Refer to `numpy.argsort` for full documentation."""

	return numpy.ndarray.argsort(self,axis,kind,order) + 1

    def take(self, indices, axis=None, out=None, mode='raise'):
        """a.take(indices, axis=None, out=None, mode='raise')
    
        Return an array formed from the elements of a at the given
        one-based indices.
        
        Refer to `numpy.take` for full documentation."""

	return numpy.ndarray.take(self,FortranArray.mapindices(indices),
				  axis,out,mode)
		
    def put(self, indices, values, mode='raise'):
        """a.put(indices, values, mode='raise')
    
        Set a.flat[n] = values[n] for all n in indices.
    
        Refer to `numpy.put` for full documentation."""

	return numpy.ndarray.put(self, FortranArray.mapindices(indices), 
				 values, mode)
				 
    def __repr__(self):
        s = numpy.asarray(self).view(numpy.ndarray).__repr__()
        s = s.replace('array','FortranArray')
        s = s.replace('\n     ','\n            ')
        return s

    def __str__(self):
        return numpy.asarray(self).view(numpy.ndarray).__str__()

    def __iter__(self):
        """Iterate over this FortranArray in the way that makes sense for a
        Fortran array, treating the first dimension as the fastest-varying."""
        if self.shape == ():
            yield self.item()
        else:
            for i in frange(self.shape[0]):
                yield self[i]

    def norm2(self):
        "array of the norm**2 of each vector in a 3xN array"
        n, m = self.shape
        if n != 3:
          raise ValueError('first array dimension should be of size 3')
        out = fzeros(m)
        for i in frange(m):
          out[i] = numpy.dot(self[:,i],self[:,i])
        return out

    def norm(self):
       "Return sqrt(norm2(a))"
       return numpy.sqrt(self.norm2())

    def row_iter(self):
        """Iterator for MxN arrays to return rows [:,i] for i=1,N one by one
        as Mx1 arrays."""
        if self.shape == ():
            yield self.item()
        else:
            for i in frange(self.shape[-1]):
                yield self[...,i]       


