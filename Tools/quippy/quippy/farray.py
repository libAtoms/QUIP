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

"""Contains :class:`FortranArray` class and utility functions for handling one-based
array indexing."""

import sys, numpy

major, minor = sys.version_info[0:2]

assert (major, minor) >= (2, 4)

if (major, minor) < (2, 5):
    all = lambda seq: not False in seq
    any = lambda seq: True in seq


def frange(min,max=None,step=1):
    """Fortran equivalent of :func:`range` builtin.

    Returns an iterator for integers from ``min`` to ``max`` **inclusive**, increasing by
    ``step`` each time.

    >>> list(frange(3))
    [1, 2, 3]
    >>> list(frange(3,6,2))
    [3, 5]
    """
    if max is None:
        return xrange(1,min+1,step)
    else:
        return xrange(min,max+1,step)

def fenumerate(seq):
    """One-based equivalent of enumerate"""
    i = 1
    for s in seq:
        yield (i, s)
        i += 1

def fzeros(shape,dtype=float):
    """Create an empty :class:`FortranArray` with Fortran ordering."""
    return FortranArray(numpy.zeros(shape,dtype,order='F'))

def farray(seq, dtype=None):
    """Convert ``seq`` to a :class:`FortranArray` with Fortran ordering.

    >>> fa = farray([1,2,3])

    A copy of the data in seq will be made if necessary."""
    na = numpy.array(seq,order='F', dtype=dtype)
    return FortranArray(na)

def fidentity(n):
    """Return the ``n`` dimensional identity matrix as a :class:`FortranArray`."""

    return farray(numpy.identity(n))

def fvar(seq):
    """
    Create one or more rank-0 instances of :class:`FortranArray` and
    inject them into the global namespace.

    This is a convenience function useful for making arrays to use
    as intent(in,out) arguments to a Fortran function. A single
    string argument causes variables with one-letter names to be
    created. The new arrays are also returned.
    
    The following examples are equivalent::

       >>> fvar("abc")
       (FortranArray(0.0), FortranArray(0.0), FortranArray(0.0))
       >>> a, b, c = fvar(['a','b','c'])
    """
    
    import inspect
    frame = inspect.currentframe().f_back
    try:
        res = tuple([farray(0.0) for s in seq])
        for s, t in zip(seq, res):
            frame.f_globals[s] = t
        return res
    finally:
        del frame

class FortranArray(numpy.ndarray):
    """Subclass of :class:`numpy.ndarray` which uses Fortran-style one-based indexing.

    The first element is numbered one rather than zero and
    trying to access element zero will raise an :exc:`IndexError`
    exception. Negative indices are unchanged; -1 still refers to the
    highest index along a dimension. Slices are also Fortran-style,
    i.e. inclusive on the second element so 1:2 includes the one-based
    elements 1 and 2, equivalent to a C-style slice of 0:2. The ``self.flat``
    iterator is still indexed from zero."""

    __array_priority__ = 100.0

    def __array_finalize__(self, obj):
        self.cols = self.col_iter()
        self.rows = self.row_iter()

    def __array_wrap__(self, out, context=None):
        return out.view(FortranArray)

    def __new__(cls, input_array=None, doc=None):
        """Construct a FortanArray from input_array

        a = FortranArray(input_array=None, doc=None)

        If doc is not None, docstring of new array is set to doc."""

	self = numpy.asarray(input_array)
        self = self.view(FortranArray)
        
        if doc is not None:
            self.__doc__ = doc

        self.cols = self.col_iter()
        self.rows = self.row_iter()
	return self

    def __eq__(self, other):
        obj = numpy.ndarray.__eq__(self, other)
        if isinstance(obj, numpy.ndarray):
            return obj.view(FortranArray)
        else:
            return obj

    def __ne__(self, other):
        obj = numpy.ndarray.__ne__(self, other)
        if isinstance(obj, numpy.ndarray):
            return obj.view(FortranArray)
        else:
            return obj

    @staticmethod
    def map_int(idx):
        if idx > 0:
            return idx-1
        elif idx == 0:
            raise IndexError('index 0 not permitted - FortranArrays are one-based')
        else:
            return idx

    def mapindices(self, indx):
        """Transform from Fortran-style one-based to C-style zero based indices.

        indx can be any object that can be used to subscript
        an ndarray - scalar integer or slice, sequence of integers and
        slices or an ndarray of integer or boolean kind."""


        def nested(L):
            return [x for x in L if isinstance(x,list)] != []

        one_dim = False
        if not hasattr(indx, '__iter__'):
            one_dim = True
            indx = (indx,)

        islist = isinstance(indx, list)
        doext = False

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
            elif isinstance(idx, list):
                islist = True
                if 0 in idx: raise ValueError('Advanced slicing list must not contain index 0')
                nindx.append([FortranArray.map_int(i) for i in idx])
            else:
                raise ValueError('Unknown index object %r' % (idx,))

        if one_dim and not islist:
            if len(self.shape) > 1:
                res = (Ellipsis,nindx[0])
            else:
                res = nindx[0]
        else:
            if islist:
                #if not Ellipsis in nindx:
                #    if doext:
                #        pass
                #    elif nested(nindx):
                #        res = [Ellipsis] + nindx
                #    else:
                #        res = [Ellipsis] + [nindx]
                #else:
                res = nindx
            else:
                if len(self.shape) > len(nindx):
                    res = tuple([Ellipsis] + nindx)
                else:
                    res = tuple(nindx)

        return res

    def __getitem__(self, indx):
        "Overloaded __getitem__ which accepts one-based indices."
	indx = self.mapindices(indx)
	obj = numpy.ndarray.__getitem__(self, indx)
        if isinstance(obj, numpy.ndarray):
            fa = obj.view(FortranArray)
            return fa
        return obj

    def __setitem__(self, indx, value):
        "Overloaded __setitem__ which accepts one-based indices."

        domap = True
        doext = False
        if isinstance(indx, slice): domap = False
        if indx is Ellipsis: domap = False
        if isinstance(indx, list): doext = True
        if isinstance(indx, numpy.ndarray): doext = True

        if hasattr(indx, '__iter__'):
            if any([isinstance(x,slice) or x is Ellipsis for x in indx]): domap = False
            if any([isinstance(x,numpy.ndarray) for x in indx]): doext = True
            if any([isinstance(x,list) for x in indx]): doext = True
            if len(indx) != len(self.shape): domap = False
        elif isinstance(indx, int):
            if len(self.shape) != 1:
                domap = False
                indx = (Ellipsis, indx)

        if doext:
            domap = True
        
        if domap:
#            print 'mapping', indx
            indx = self.mapindices(indx)
            
	numpy.ndarray.__setitem__(self, indx, value)

    def __getslice__(self, i, j):
        "Overloaded __getslice__ which accepts one-based indices."
        if i != 0:
            i = FortranArray.map_int(i)
        obj = numpy.ndarray.__getslice__(self, i, j)
        if isinstance(obj, numpy.ndarray):
            fa = obj.view(FortranArray)
            return fa

    def __setslice__(self, i, j, value):
        "Overloaded __setslice__ which accepts one-based indices."
        if i != 0:
            i = FortranArray.map_int(i)
        numpy.ndarray.__setslice__(self, i, j, value)

    def nonzero(self):
        """Return the one-based indices of the elements of a which are not zero."""
	return tuple(a + 1 for a in numpy.ndarray.nonzero(self))
	
    def count(self):
        """Number of nonzero elemnts of this FortranArray.

        Equivalent to len(self[self.nonzero()])."""
	return len(self[self.nonzero()])
    
    def argmin(self, axis=None, out=None):
        """Return one-based indices of the minimum values along the given  axis of `a`.
        
        Refer to `numpy.ndarray.argmax` for detailed documentation."""
	if axis is not None and axis > 0:
	    axis -= 1
	return numpy.ndarray.argmin(self,axis,out) + 1

    def argmax(self, axis=None, out=None):
        """Return one-based indices of the maximum values along the given axis of `a`.
        
        Refer to `numpy.ndarray.argmax` for detailed documentation."""
	if axis is not None and axis > 0:
	    axis -= 1
	return numpy.ndarray.argmax(self,axis,out) + 1	

    def argsort(self, axis=None, kind='quicksort', order=None):
        """Returns the indices that would sort this array.
        
        Refer to `numpy.argsort` for full documentation."""

        if axis is not None and axis > 0:
            axis -= 1
	return numpy.ndarray.argsort(self,axis,kind,order) + 1

    def take(self, indices, axis=None, out=None, mode='raise'):
        """Return an array formed from the elements of a at the given
        one-based indices.
        
        Refer to `numpy.take` for full documentation."""

        if axis is not None and axis > 0:
            axis -= 1
	return numpy.ndarray.take(self,self.mapindices(indices),
				  axis,out,mode)
		
    def put(self, indices, values, mode='raise'):
        """Set a.flat[n] = values[n] for all n in indices.
    
        Refer to `numpy.put` for full documentation."""

	return numpy.ndarray.put(self, self.mapindices(indices), 
				 values, mode)
				 
    def __repr__(self):
        s = repr(numpy.asarray(self).view(numpy.ndarray))
        s = s.replace('array','FortranArray')
        s = s.replace('\n     ','\n            ')
        return s
        

    def __str__(self):
        return str(numpy.asarray(self).view(numpy.ndarray))

    def __iter__(self):
        """Iterate over this :class:`FortranArray` treating first dimension as fastest varying.

        Calls fast :meth:`ndarray.__iter__` for a 1D array."""

        if len(self.shape) > 1:
            return self.col_iter()
        else:
            return numpy.ndarray.__iter__(numpy.asarray(self).view(numpy.ndarray))


    def row_iter(self):
        """Iterate over this :class:`FortranArray` treating first dimension as fastest varying"""
        if self.shape == ():
            yield self.item()
        else:
            for i in frange(self.shape[0]):
                obj = numpy.ndarray.__getitem__(self, i-1)
                if (isinstance(obj, numpy.ndarray) and obj.dtype.isbuiltin):
                    fa = obj.view(FortranArray)
                    yield fa
                else:
                    yield obj
                
    def norm2(self):
        """Squared norm of a 1D or 2D array.

        For a 1D array, returns dot(self,self)
        For a 2D array, must have shape (3,n). Returns array a where a[i] = dot(self[:,i],self[:,i])"""
        if len(self.shape) == 2:
            n, m = self.shape
            if n != 3:
                raise ValueError('first array dimension should be of size 3')
            out = fzeros(m)
            for i in frange(m):
                out[i] = numpy.dot(self[:,i],self[:,i])
            return out
        elif len(self.shape) == 1:
            return numpy.dot(self,self)
        elif len(self.shape) == 0:
            return self.item()
        else:
            raise ValueError("Don't know how to take norm2 of array with shape %s" % str(self.shape))
            

    def norm(self):
       "Return ``sqrt(norm2(self))``"
       return numpy.sqrt(self.norm2())

    def col_iter(self):
        """Iterator for ``MxN`` arrays to return ``cols [...,i]`` for ``i=1,N`` one by one as ``Mx1`` arrays."""
        if self.shape == ():
            yield self.item()
        else:
            for i in frange(self.shape[-1]):
                obj = numpy.ndarray.__getitem__(self, (Ellipsis, i-1)).view(FortranArray)
                yield obj

    def all(self, axis=None, out=None):
        """One-based analogue of :meth:`numpy.ndarray.all`"""
        
	if axis is not None and axis > 0:
	    axis -= 1
	obj = numpy.ndarray.all(self, axis, out)
        if isinstance(obj, numpy.ndarray):
            obj = obj.view(FortranArray)
        return obj

    def any(self, axis=None, out=None):
        """One-based analogue of :meth:`numpy.ndarray.any`"""
        
	if axis is not None and axis > 0:
	    axis -= 1
	obj = numpy.ndarray.any(self, axis, out).view(FortranArray)
        if isinstance(obj, numpy.ndarray):
            obj = obj.view(FortranArray)
        return obj

    def sum(self, axis=None, dtype=None, out=None):
        """One-based analogue of :meth:`numpy.ndarray.sum`"""
	if axis is not None and axis > 0:
	    axis -= 1
	obj = numpy.ndarray.sum(self, axis, out).view(FortranArray)
        if isinstance(obj, numpy.ndarray):
            obj = obj.view(FortranArray)
        return obj

    def mean(self, axis=None, dtype=None, out=None):
        """One-based analogue of :meth:`numpy.ndarray.mean`"""
	if axis is not None and axis > 0:
	    axis -= 1
	obj = numpy.ndarray.mean(self, axis, out).view(FortranArray)
        if isinstance(obj, numpy.ndarray):
            obj = obj.view(FortranArray)
        return obj


    def stripstrings(self):
        """Return contents as string (0- and 1-dimensional arrays) or array of strings with
        trailing spaces removed (2-dimensional arrays).

        Raises :exc:`ValueError` if this :class:`FortranArray` does not have a string datatype.
        """

        if self.dtype.kind != 'S': raise ValueError('dtype.kind must be "S"')
        if len(self.shape) == 0:
            return self.item()
        elif len(self.shape) == 1:
            return ''.join(self).strip()
        else:
            return farray([''.join(x).strip() for x in self])


def padded_str_array(d, length):
    """Return :class:`FortranArray` with shape ``(length, len(d))``,
       filled with rows from ``d`` padded with spaces"""
    res = fzeros((length, len(d)), 'S')
    res[...] = ' '
    for i, line in fenumerate(d):
        res[1:len(line),i] = list(line)
    return res


def convert_farray_to_ndarray(func):
   """Decorator to convert all occurences of farray in arguments to ndarray.
      If result contains ndarray, they will be converted back to farray."""
   from functools import wraps
   @wraps(func)

   def nfunc(*args, **kwargs):
       from numpy  import ndarray
       from farray import FortranArray

       nargs = []
       for a in args:
           if isinstance(a, FortranArray):
               nargs.append(a.view(ndarray))
           else:
               nargs.append(a)

       nkwargs = {}
       for k, v in kwargs.iteritems():
           if isinstance(v, FortranArray):
               nkwargs[k] = v.view(ndarray)
           else:
               nkwargs[k] = v

       nargs = tuple(nargs)
       res = func(*nargs, **nkwargs)
       if res is None: return
       
       one_result = not isinstance(res, tuple)
       if one_result:
          res = (res,)

       nres = []
       for r in res:
           if isinstance(r, ndarray):
               nres.append(r.view(FortranArray))
           else:
               nres.append(r)

       if one_result:
          return nres[0]
       else:
          return tuple(nres)
          

   return nfunc


def convert_ndarray_to_farray(func):
   """Decorator to convert all occurences of ndarray in arguments to farray.
   If results contains farray, they will be converted back to ndarray. """
   from functools import wraps
   @wraps(func)

   def nfunc(*args, **kwargs):
       from numpy  import ndarray
       from farray import FortranArray

       nargs = []
       for a in args:
           if isinstance(a, ndarray):
               nargs.append(a.view(FortranArray))
           else:
               nargs.append(a)

       nkwargs = {}
       for k, v in kwargs.iteritems():
           if isinstance(v, ndarray):
               nkwargs[k] = v.view(FortranArray)
           else:
               nkwargs[k] = v

       nargs = tuple(nargs)
       res = func(*nargs, **nkwargs)
       if res is None: return

       one_result = not isinstance(res, tuple)
       if one_result:
          res = (res,)

       nres = []
       for r in res:
           if isinstance(r, FortranArray):
               nres.append(r.view(ndarray))
           else:
               nres.append(r)

       if one_result:
          return nres[0]
       else:
          return tuple(nres)

   return nfunc

def f2n(x):
    """Return ``x.view(numpy.ndarray)``"""
    return x.view(numpy.ndarray)

def n2f(x):
    """Return ``x.view(FortranArray)``"""
    return x.view(FortranArray)
    
def unravel_index(x, dims):
    """1-based version of numpy.unravel_index"""

    return tuple([n+1 for n in numpy.unravel_index(x-1, dims)])


def s2a(d, pad=None):
    """Convert from list of strings to array of shape (pad, len(d)).
    If pad is not specified, we use the length of longest string in d"""
    if pad is None:
        pad = max([len(s) for s in d])
    res = fzeros((pad, len(d)), 'S')
    res[...] = ' '
    for i, line in fenumerate(d):
        res[1:len(line),i] = list(line)
    return res

def a2s(a):
    """Convert from array of shape (pad, len(d)) to list of strings d"""
    return [ ''.join(a[:,i]).strip() for i in frange(a.shape[1]) ]

def loadtxt(filename):
    return farray(numpy.loadtxt(filename))

def loadcsv(filename):
    """Read CSV formatted data file and return dictionary of farrays, using first
       row in file as column labels for dictionary keys."""
    data = numpy.loadtxt(filename,delimiter=',',skiprows=1)
    cols = open(filename,'rU').readline().strip().split(',')
    data = dict(zip(cols, (farray(data[:,i]) for i in range(data.shape[1]))))
    return data
