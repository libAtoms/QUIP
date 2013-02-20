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

"""Object oriented interface on top of f2py generated wrappers.

This module adds support for derived types to f2py generated modules,
providing the Fortran source code was built using
:mod:`f2py_wrapper_gen`."""

import logging
import inspect
import sys
import imp
import weakref
import re
from types import MethodType

import numpy as np

from quippy.farray import *
from quippy.util import args_str, is_interactive_shell
from quippy.dictmixin import PuPyDictionary
from quippy import fortran_indexing # default value
import quippy.arraydata as arraydata

__all__ = ['FortranDerivedType', 'FortranDerivedTypes', 'wrap_all']

major, minor = sys.version_info[0:2]

if (major, minor) < (2, 5):
    all = lambda seq: not False in seq
    any = lambda seq: True in seq

wraplog = logging.getLogger('quippy.oo_fortran')

fortran_class_prefix = ''

py_keywords = ['and',       'del',       'from',      'not',       'while',
               'as',        'elif',      'global',    'or',        'with',
               'assert',    'else',      'if',        'pass',      'yield',
               'break',     'except',    'import',    'print',
               'class',     'exec',      'in',        'raise',
               'continue',  'finally',   'is',        'return',
               'def',       'for',       'lambda',    'try']

py_keywords_map = dict([(k,k+'_') for k in py_keywords])

numpy_to_fortran = {
    'd': 'real(dp)',
    'i': 'integer',
    'b': 'logical',
    'c': 'complex(dp)',
    'S': 'character',
    'f': 'real(dp)'
    }

FortranDerivedTypes = {}
FortranRoutines = {}

from quippy import QUIPPY_TRUE, QUIPPY_FALSE

from numpy.f2py.crackfortran import badnames

special_names = {
    'assignment(=)' : 'assignment',
    }

rev_special_names = dict((v,k) for (k,v) in special_names.iteritems())


def is_scalar_type(t):
    if t in numpy_to_fortran.values(): return True
    if t.startswith('character'): return True
    if t.startswith('integer'): return True
    return False

def fortran_equivalent_type(t):
    tmap = {'int': 'integer',
            'float': 'real(dp)',
            'bool': 'logical',
            'str': 'character',
            'complex' : 'complex(dp)'
            }

    if type(t).__name__ in tmap:
        return (tmap[type(t).__name__], 'scalar', 1)
    elif isinstance(t, FortranDerivedType):
        for typename,cls in FortranDerivedTypes.iteritems():
            if isinstance(t, cls):
                return (typename, 'scalar', 1)
        else:
            raise TypeError('Unknown derived type %s' % t)
    elif hasattr(t,'__iter__') or hasattr(t,'dtype'):
        a = np.array(t)
        if a.dtype.kind in numpy_to_fortran:
            return (numpy_to_fortran[a.dtype.kind], a.shape, a.dtype.itemsize)
        else:
            raise TypeError('Unknown array type %s' % a.dtype.kind)
    else:
        raise TypeError('Unknown type %s' % type(t))

def python_equivalent_type(t):
    tmap = { 'integer': 'int',
             'real(dp)': 'float',
             'logical': 'bool',
             'character': 'str',
             'complex(dp)' : 'complex',
             '' : 'None'}
    t = tmap.get(t,t)
    if t.startswith('character('):
        t = 'str'
    return t
    

def type_is_compatible(spec, arg):
    if 'optional' in spec['attributes'] and arg is None:
        return True

    try:
        arg_type, dims, itemsize = fortran_equivalent_type(arg)
    except TypeError:
        return False

    if dims == (): dims = 'scalar'

    spec_type = spec['type'].lower()
    if spec_type.startswith('character'):
        spec_type = 'character'

    # First check type matches
    if arg_type != spec_type: return False

    # Now check dimensions
    dimattrs = [s for s in spec['attributes'] if s.startswith('dimension')]

    if dims == 'scalar':
        return len(dimattrs) == 0
    else:
        try:
            (dimattr,) = dimattrs # asserts len(dimattr) == 1
        except ValueError:
            return False

        adims = dimattr[dimattr.index('(')+1:dimattr.index(')')].split(',')

        # See if dimensions are all numeric
        try:
            adims = [int(a) for a in adims]
        except ValueError:
            # Just check number of dimensions
            if spec_type != 'character':
                return len(adims) == len(dims)
            else:
                return ((itemsize != 1 and len(dims) == len(adims)) or   # arrays of strings, dtype='S'
                        (itemsize == 1 and (len(dims) == len(adims)+1 or len(dims) == len(adims))))   # arrays of characters, dtype='S1'
        else:
            # Check each dimension matches
            return all([x == y for (x,y) in zip(adims,dims)])

def process_in_args(args, kwargs, inargs, prefix):
    # Process positional arguments
    newargs = []
    modargs = []
    for arg, spec in zip(args,inargs):
        if arg is not None and spec['type'].startswith('type'):
            cls = FortranDerivedTypes[spec['type'].lower()]
            if isinstance(arg, cls):
                newargs.append(arg._fpointer)
                modargs.append(None)
            else:
                try:
                    raise NotImplementedError
                    newarg = cls(arg)
                    newargs.append(newarg._fpointer)
                    modargs.append(newarg)
                except:
                    TypeError('Argument %s should be of type %s but got incompatible type %s' % (arg, spec['type'], type(arg)))

        elif arg is not None and spec['type'].startswith('character'):
            # if arg is a list of strings of unequal length, pad with spaces
            if isinstance(arg, list) and any([len(x) != len(arg[0]) for x in arg]):
                a = s2a(arg).T
                newargs.append(a)
                modargs.append(None)
            else:
                newargs.append(arg)
                modargs.append(None)
        else:
            newargs.append(arg)
            modargs.append(None)

    # Process keyword arguments and args_str arguments.

    # If we're expecting an args_str argument, then any kwarg is
    # permitted. Otherwise kwarg name must match one of the expected
    # names. If kwarg name matches one of expected values, but type of
    # argument does not, then we assume it's an arg_str argument.

    kwarg_lookup = dict([(badnames.get(x['name'].lower(),x['name'].lower()),x) for x in inargs])
    got_args_str = prefix+'args_str' in kwarg_lookup

    newkwargs = {}
    modkwargs = {}
    args_str_kwargs = {}
    for k,a in kwargs.iteritems():
        k = prefix+k
        if got_args_str:
            if k != prefix+'args_str' and (k not in kwarg_lookup or (k in kwarg_lookup and not type_is_compatible(kwarg_lookup[k], a))):
                args_str_kwargs[k[len(prefix):]] = a
                continue
        if k not in kwarg_lookup:
            raise ValueError('Unknown keyword argument %s' % k)
        if kwarg_lookup[k]['type'].startswith('type'):
            cls = FortranDerivedTypes[kwarg_lookup[k]['type'].lower()]
            if isinstance(a, cls):
                newkwargs[k] = a._fpointer
                modkwargs[k] = None
            elif a is None:
                continue
            else:
                try:
                    raise NotImplementedError
                    na = cls(a)
                    newkwargs[k] = na._fpointer
                    modkwargs[k] = na
                except:
                    raise TypeError('Argument %s should be of type %s, but got incompatible type %s' % (k,a,kwarg_lookup[k],type(a)))
        else:
            if a is None: continue
            newkwargs[k] = a
            modkwargs[k] = None

    # Construct final args_str by merging args_str argument with args_str_kwargs
    if got_args_str:
        args_str_final = ''
        if prefix+'args_str' in newkwargs:
            if isinstance(newkwargs[prefix+'args_str'], basestring):
                args_str_final = args_str_final + " " + newkwargs[prefix+'args_str']
            else:
                args_str_final = args_str_final + " " + args_str(newkwargs[prefix+'args_str'])
        if args_str_kwargs != {}:
            args_str_final = args_str_final + " " + args_str(args_str_kwargs)
        if args_str_final != '':
            newkwargs[prefix+'args_str'] = args_str_final

    return tuple(newargs), tuple(modargs), newkwargs, modkwargs

def process_results(res, args, kwargs, inargs, outargs, prefix, fortran_indexing):
    newres = []

    madeseq = False
    if res is not None and len(outargs) <= 1:
        madeseq = True
        res = (res,)

    # intent(out) arguments form result tuple
    if res is not None:
        for r, spec in zip(res,outargs):
            if spec['type'].startswith('type'):
                newres.append(FortranDerivedTypes[spec['type'].lower()](fpointer=r,finalise=True))
            elif isinstance(r, np.ndarray):
                if fortran_indexing:
                    # Convert to one-based FortranArray
                    r = r.view(FortranArray)
                newres.append(r)
            elif spec['type'] == 'logical':
                newres.append(r == QUIPPY_TRUE)
            else:
                newres.append(r)

    # update any objects in args or kwargs affected by this call
    for arg, spec in zip(args, inargs):
        if (isinstance(arg, FortranDerivedType) and not 'fintent(in)' in spec['attributes']):
            arg._update()

    type_lookup = dict([(badnames.get(x['name'].lower(),x['name'].lower()),x) for x in inargs])
    for k,a in kwargs.iteritems():
        k = prefix + k
        if (isinstance(a, FortranDerivedType) and not 'fintent(in)' in type_lookup[k]['attributes']):
            a._update()

    if res is None:
        return None
    elif madeseq:
        return newres[0]
    else:
        return tuple(newres)


class FortranDerivedType(object):
    """Abstract base class for all fortran derived types.

    This class is an abstract base class for all Fortran
    derived-types. It contains all the magic necessary
    to communicate with a module wrapped using the
    :mod:`f2py_wrapper_gen` module.

    The constructor invokes a Fortran interface or routine to initialise
    the object. When the reference count for instances of this class drops to zero,
    the destructor frees the associated Fortran derived type.
    """

    _modobj = None
    _moddoc = None
    _classdoc = None
    _routines = {}
    _subobjs = {}
    _arrays = {}
    _interfaces = {}
    _elements = {}
    _cmp_skip_fields = []
    _cmp_tol = 1e-8
    _fortran_indexing = fortran_indexing

    _prefix = ''

    def __init__(self, *args, **kwargs):

        call_init = True
        self._fpointer = None
        self._finalise = True
        self._update_hooks = []
        self._subobjs_cache = {}

        if 'fpointer' in kwargs:
            if kwargs['fpointer'] is not None:
                call_init = False
                self._fpointer = kwargs['fpointer']
            del kwargs['fpointer']

        if 'finalise' in kwargs:
            self._finalise = kwargs['finalise']
            del kwargs['finalise']

        if 'fortran_indexing' in kwargs:
            self.fortran_indexing = kwargs['fortran_indexing']
            #print 'setting %s fortran_indexing id %d to %r' % (self.__class__.__name__, id(self), self.fortran_indexing)
            del kwargs['fortran_indexing']

        orig_args = args
        orig_kwargs = kwargs.copy()

        if call_init:
            if '__init__' in self._interfaces:
                self._fpointer = self._runinterface('__init__', *args, **kwargs)
            elif '__init__' in self._routines:
                self._fpointer = self._runroutine('__init__', *args, **kwargs)

        wraplog.debug('Constructing %s(fpointer=%r, finalise=%r)' % (self.__class__.__name__, self._fpointer, self._finalise))
        self._update()

    def __len__(self):
        return self.n        

    def _set_fortran_indexing(self, fortran_indexing):
        self._fortran_indexing = fortran_indexing
        self._subobjs_cache = {}

    def _get_fortran_indexing(self):
        return self._fortran_indexing

    fortran_indexing = property(_get_fortran_indexing, _set_fortran_indexing)

    def is_same_fortran_object(self, other):
        """Test if `self` and `other` point to the same Fortan object."""
        if not hasattr(other, '_fpointer'):
            return False
        return (self._fpointer == other._fpointer).all()

    def shallow_copy(self):
        """Return a shallow copy of `self`."""
        other = self.__class__()
        other.shallow_copy_from(self)
        return other

    def shallow_copy_from(self, other):
        """Transform `self` into a shallow copy of `other`."""

        # free any memory currently associated with `self`
        self.__class__.__del__(self)

        # copy any normal (not Fortran) attributes
        for k, v in other.__dict__.iteritems():
            if not k.startswith('_') and k not in self.__dict__:
                self.__dict__[k] = v

        self._fpointer = other._fpointer
        if hasattr(self, 'ref_count'):
            self.ref_count = self.ref_count + 1

    def __del__(self):
        try:
            wraplog.debug('Freeing %s(fpointer=%r, finalise=%r)' % (self.__class__.__name__, self._fpointer, self._finalise))
        except:
            pass
            #print 'Freeing %s(fpointer=%r, finalise=%r)' % (self.__class__.__name__, self._fpointer, self._finalise)

        if hasattr(self, '_initialised'):
            del self._initialised
        self._subobjs_cache = {}
        if self._fpointer is not None and not (self._fpointer == 0).all() and self._finalise and '__del__' in self._routines:
            fobj, doc = self._routines['__del__']
            fobj(self._fpointer)
            self._fpointer = None


    def __repr__(self):
        if self._fpointer is None:
            return '<%s object at 0x%x fpointer=None>' % (self.__class__.__name__, id(self))
        else:
            return '<%s object at 0x%x fpointer=%r>' % (self.__class__.__name__, id(self), tuple(self._fpointer[self._fpointer != 0]))

    def __str__(self):
        items = []
        for k in dir(self):
            if k.startswith('_'): continue
            v = getattr(self,k)
            if isinstance(v, MethodType): continue
            if isinstance(v, FortranDerivedType):
                items.append('%s=%s()' % (k,v.__class__.__name__))
            else:
                if isinstance(v, str):
                    v = v.strip()
                items.append('%s=%s' % (k,v))
        return '%s(%s)' % (self.__class__.__name__, ',\n'.join(items))

    def __eq__(self, other):
        # Class must agree
        if not issubclass(self.__class__, other.__class__) and not issubclass(other.__class__, self.__class__):
            wraplog.debug('class mismatch %s not subclass of %s (or vice versa)' % (self.__class__.__name__, other.__class__.__name__))
            return False

        # Don't compare uninitialised objects
        if hasattr(self, 'initialised'):
            if not self.initialised and not other.initialised:
                return True
            if ((self.initialised and not other.initialised) or
                (other.initialised and not self.initialised)):
                return False

        # Compare elements and sub-objects
        for el in self._elements.keys() + self._subobjs.keys():
            if el in self._cmp_skip_fields:
                continue

            if getattr(self, el) != getattr(other, el):
                wraplog.debug('element mismatch %s' % el)
                return False

        # Compare arrays
        for array in self._arrays:
            a = getattr(self, array)
            b = getattr(other, array)
            if a is None and b is None: continue
            if a is None or b is None: return False
            if a.dtype.kind != 'f':
                if not (a == b).all():
                    wraplog.debug('array mismatch %s' % array)
                    return False
            else:
                if (abs(a-b) > self._cmp_tol).any():
                    wraplog.debug('real array mismatch %s' % array)
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def _update(self):
        """
        Automatically invoked whenever this object needs to be
        updated. This happens when it is first created, when it is
        direcly passed to a Fortran routine with ``intent(in,out)`` or
        when it is a component of another object which is passed with
        ``intent(in,out)``.

        This method also invokes all the hook in the
        :attr:`_update_hooks` list. This can be used in subclasses to
        allow a customised response. For example this mechanism is used
        in :class:`quippy.extras.Atoms` to update Atoms properties.
        """
        
        wraplog.debug('updating %s at 0x%x' % (self.__class__, id(self)))
        if self._fpointer is None: return
        for hook in self._update_hooks:
            hook(self)

    def add_hook(self, hook, *args):
        self._update_hooks.append(hook)

    def remove_hook(self, hook, *args):
        self._update_hooks.remove(hook)

    def _get_array_shape(self, name):
        """
        This method can be used to override Fortran's idea of the shape of 
        arrays within derived types, for example to present only a partial
        view of an array. This is used in :class:`quippy.table.Table`
        to allow the sizes of arrays within the Table class to correspond to the
        current extent of the Table, rather than the size of the allocated storage
        which will usually be larger.

        If this method returns ``None`` then the full array is presented, otherwise
        the return value should be a tuple `(N_1, N_2, ..., N_d)` where `d`
        is the number of dimensions of the array and `N_1, N_2` etc. are the lengths
        of each dimension.
        """
        return None

    def _runroutine(self, name, *args, **kwargs):
        """
        Internal method used to invoke the Fortran routine `name`

        `name` must be a valid key in :attr:`_routines`. Wrapper
        methods which simply call :meth:`_runroutine()` are
        automatically generated in subclasses of
        :class:`FortranDerivedType` by :func:`wrap_all()`.

        Input arguments which are instances of a subclass of
        :class:`FortranDerivedType` are replaced by their
        :attr:`_fpointer` integer attribute.

        If there is an keyword argument with the name `args_str`
        then unexpected keyword arguments are permitted. All the
        undefined keyword arguments are collected together to form a
        dictionary which is converted to string form and used as the the
        `arg_string` argument, providing rudimentary support for
        variable numbers and types of arguments. For example::

            p = Potential('IP SW', xml_string)
	    p.calc(at, virial=True, energy=True)

        is equivalent to::

            p = Potential('IP SW', xml_string)
	    p.calc(at, args_str="virial energy")
	
        The return value us made up of a tuple of the arguments to the
        Fortran routine which are ``intent(out)``. Pointers to Fortran
        derived-type instances are replaced with new instances of the
        appropriate subclass of :class:`FortranDerivedType`. Arrays
        are converted to use one-based indexing using
        :class:`~quippy.farray.FortranArray`.
        """
        
        if not name.startswith('__init__') and self._fpointer is None:
            raise ValueError('%s object not initialised.' % self.__class__.__name__)

        if not name in self._routines:
            raise NameError('Unknown fortran routine: %s' % name)

        fobj, doc = self._routines[name]

        #inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'])
        #outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])

        inargs  = [ x for x in doc['args'] if not 'intent(out)' in x['attributes'] ]
        outargs = [ x for x in doc['args'] if 'intent(out)' in x['attributes'] ]        
        
        if not name.startswith('__init__'):
            # Put self at beginning of args list
            args = tuple([self] + list(args))
        newargs, modargs, newkwargs, modkwargs = process_in_args(args, kwargs, inargs, self._prefix)

        try:
            res = fobj(*newargs, **newkwargs)
        except RuntimeError:
            raise
            try:
                exctype, value, tb = sys.exc_info()

                if is_interactive_shell():
                    error_str = 'Fortran routine %s(%s):\n %s' % (name, ', '.join(list('%r' % a for a in args) + ['%s=%r' % (k,v) for (k,v) in kwargs.iteritems()]), str(value).strip())
                else:
                    # Remove Fortran traceback
                    error_str = str(value).strip().split('\n')[-1].strip()

                raise exctype, RuntimeError(error_str), tb

            finally:
                del tb
        except:
            raise

        if not name.startswith('__init__'):
            newres = process_results(res, args, kwargs, inargs, outargs, self._prefix, self.fortran_indexing)
        else:
            newres = res

        return newres


    def _runinterface(self, name, *args, **kwargs):
        """
        Internal method used to invoke the appropriate routine
        within the Fortran interface `name`. If no routine is found 
        matching the names and types of the arguments provided then
        an :exc:`TypeError` exception is raised.

        Arguments and results are handled in the same way as 
        :func:`_runroutine`.
        """

        if not name in self._interfaces:
            raise ValueError('Unknown interface %s.%s' % (self.__class__.__name__, name))

        wraplog.debug('Interface %s %s' % (self.__class__.__name__, name))

        for rname, spec, routine in self._interfaces[name]:

            wraplog.debug('Trying candidate routine %s' % rname)

            #inargs = filter(lambda x: not 'intent(out)' in x['attributes'], spec['args'])
            inargs  = [ x for x in spec['args'] if not 'intent(out)' in x['attributes'] ]

            if not name.startswith('__init__'):
                inargs = inargs[1:] # remove self argument

            oblig_inargs = [x for x in inargs if  not 'optional' in x['attributes'] ]
            opt_inargs   = [x for x in inargs if 'optional' in x['attributes'] ]

            # Check number of arguments is compatible
            # Should be in range len(oblig_inargs) <= L <= len(oblig_inargs) + len(opt_inargs)
            if (len(args)+len(kwargs) < len(oblig_inargs) or
                len(args)+len(kwargs) > len(oblig_inargs) + len(opt_inargs)):
                wraplog.debug('Number of arguments incompatible: %d must be in range %d <= n <= %d' %
                              (len(args), len(oblig_inargs), len(oblig_inargs)+len(opt_inargs)))
                continue

            newinargs = (oblig_inargs + opt_inargs)[:len(args)]

            # Check types and dimensions are compatible
            if not all([type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]):
                wraplog.debug('Types and dimensions of args incompatible %s %s %s' %
                              (newinargs, args, [type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]))
                continue

            # Check keyword arguments, if incompatible continue to next
            # candidate interface
            if kwargs:
                innames = [badnames.get(x['name'],x['name']).lower() for x in oblig_inargs + opt_inargs]

            if not all([self._prefix+key.lower() in innames[len(args):] for key in kwargs.keys()]):
                wraplog.debug('Unexpected keyword argument valid=%s got=%s' % (kwargs.keys(), innames[len(args):]))
                continue

            try:
                for key, arg in kwargs.iteritems():
                    (inarg,) = [x for x in inargs if badnames.get(x['name'],x['name']).lower() == self._prefix+key.lower()]
                    if not type_is_compatible(inarg, arg):
                        wraplog.debug('Types and dimensions of kwarg %s incompatible' % key)
                        raise ValueError

            except ValueError:
                continue

            wraplog.debug('calling '+rname)
            return routine(self, *args, **kwargs)


        raise TypeError('No matching routine found in interface %s' % name)



    def _runinterface_simple(self, name, *args, **kwargs):
        if not name in self._interfaces:
            raise ValueError('Unknown interface %s.%s' % (self.__class__.__name__, name))

        wraplog.debug('Interface %s %s' % (self.__class__.__name__, name))

        for rname, spec, routine in self._interfaces[name]:
            wraplog.debug('Trying candidate routine %s' % rname)
            try:
                return routine(self, *args, **kwargs)
            except Exception, e:
                if isinstance(e, RuntimeError):
                    raise
                else:
                    continue

        raise TypeError('No matching routine found in interface %s' % name)


def nested_dict_update(mod, w):
    if type(mod) == list:
        if type(w) == list:
            mod += w
        else:
            raise RuntimeError('First argument is list, but second is not.')
    elif type(mod) == str:
        mod += '\n' + w
    else:

        for key, value in w.iteritems():
            if key in mod:
                nested_dict_update(mod[key], value)
            else:
                mod[key] = value


def flatten_list_of_dicts(dictlist):
    if type(dictlist) == dict:
        return dictlist

    modout = { }
    for d in dictlist:
        nested_dict_update(modout, d)
    return modout


def wrap_all(fobj, spec, mods, merge_mods, short_names, prefix, package, modules_name_map, fortran_indexing):
    """
    Returns tuple `(classes, routines, params)` suitable for
    importing into top-level package namespace. `topmod` should be an
    f2py-generated module containing `fortran` objects, and `spec`
    should be the specification dictionary generated by
    :func:`f2py_wrapper_gen.wrap_mod`.  `mods` is a list of the names
    of Fortran modules to wrap, and `short_names` is a dictionary
    mapping shortened Fortran derived-type names to their canonical
    form.

    `classes` and `routines` are lists of `(name, value)` tuples
    where `value` is a newly defined subclass of
    :class:`FortranDerivedType` or newly wrapped routine respectively.
    `params` is a Python dictionary of Fortran parameters (constants).

    Here's how this function is used in quippy's :file:`__init.py__` to
    import the new classes, routines and params into the top-level quippy
    namespace::
   
       classes, routines, params = wrap_all(_quippy, spec, spec['wrap_modules'], spec['short_names'])

       for name, cls in classes:
 	  setattr(sys.modules[__name__], name, cls)

       for name, routine in routines:
   	  setattr(sys.modules[__name__], name, routine)

       sys.modules[__name__].__dict__.update(params)
    """

    all_classes = []
    leftover_routines = []
    interfaces = {}
    top_level_routines = []
    all_params = {}
    pymods = {}

    for mod, modfile in mods:
        wraplog.debug('Module '+str(mod))
        if mod in merge_mods:
            curspec = flatten_list_of_dicts([spec[x] for x in merge_mods[mod]])
        else:
            curspec = spec[mod]

        pymodname = package+'.'+modules_name_map.get(mod,mod)

        classes, routines, mod_params = wrapmod(fobj, curspec, modname=mod, modfile=modfile,
                                                short_names=short_names,
                                                params=all_params, prefix=prefix,
                                                fortran_indexing=fortran_indexing,
                                                pymodname=pymodname)
        all_classes.extend(classes)

        pymod = imp.new_module(pymodname)
        pymod.__doc__ = process_docstring(curspec['doc'])
        pymod.__package__ = package
        
        for name, cls in classes:
            cls.__module__ = pymod.__name__
            setattr(pymod, name, cls)

        pymod.__dict__.update(mod_params)
        pymods[modules_name_map.get(mod,mod)] = pymod
        
        # make classes and params symbols public
        # (available with 'from ... import *')
        setattr(pymod, '__all__', [c[0] for c in classes]+mod_params.keys())
                
        leftover_routines.append((mod, curspec, routines, modfile))

    # add orphaned routines to classes and generate top-level interfaces
    for mod, modspec, routines, modfile in leftover_routines:
        for routine in routines:
            rspec = modspec['routines'][routine]

            if (len(rspec['args']) > 0 and
                rspec['args'][0]['type'].lower() in FortranDerivedTypes and
                rspec['args'][0]['name'].lower() == prefix+'this'):

                cls = FortranDerivedTypes[rspec['args'][0]['type'].lower()]

                lclsname = cls.__name__.lower()[len(fortran_class_prefix):]

                if routine.startswith(lclsname+'_'):
                    method_name = routine[len(lclsname)+1:]
                elif lclsname in short_names and routine.startswith(short_names[lclsname]+'_'):
                    method_name = routine[len(short_names[lclsname])+1:]
                else:
                    method_name = routine

                method_name = py_keywords_map.get(method_name, method_name)
                method_name = special_names.get(method_name, method_name)

                wrapped_routine = wraproutine(fobj, modspec, routine, cls.__name__+'.'+method_name, prefix,
                                              fortran_indexing=fortran_indexing, skip_this=True, modfile=modfile)
                FortranRoutines[cls.__name__+'.'+method_name] = wrapped_routine
                setattr(cls, method_name, wrapped_routine)
                wraplog.debug('  added method %s to class %s' % (method_name, cls.__name__))

            else:
                wrapped_routine = wraproutine(fobj, modspec, routine, routine, prefix,
                                              fortran_indexing=fortran_indexing, skip_this=True, modfile=modfile)
                for intf_name,intf_spec in modspec['interfaces'].iteritems():
                    if routine in intf_spec['routines']:
                        if not intf_name in interfaces: interfaces[intf_name] = (mod, intf_spec, [], modfile)
                        interfaces[intf_name][2].append((routine, rspec, wrapped_routine, modfile))
                        wraplog.debug('  added routine %s to top-level interface %s' % (routine, intf_name))
                        break
                else:
                    FortranRoutines[routine] = wrapped_routine
                    top_level_routines.append((routine,wrapped_routine))
                    wrapped_routine.__module__ = pymods[modules_name_map.get(mod,mod)].__name__
                    setattr(pymods[modules_name_map.get(mod,mod)], routine, wrapped_routine)
                    pymods[modules_name_map.get(mod,mod)].__all__.append(routine)
                    wraplog.debug('  added top-level routine %s' % routine)
                    
    # remap interface names which clash with keywords and skip overloaded operators
    interfaces = dict([(py_keywords_map.get(k,k),v)
                         for (k,v) in interfaces.iteritems() ])
    interfaces = dict([(special_names.get(k,k),v)
                         for (k,v) in interfaces.iteritems() ])

    for name, (mod, intf_spec, routines, modfile) in interfaces.iteritems():
        wrapped_interface = wrapinterface(name, intf_spec, routines, prefix, modfile=modfile)
        FortranRoutines[name] = wrapped_interface
        top_level_routines.append((name, wrapped_interface))
        wrapped_interface.__module__ = pymods[modules_name_map.get(mod,mod)].__name__
        setattr(pymods[modules_name_map.get(mod,mod)], name, wrapped_interface)
        pymods[modules_name_map.get(mod,mod)].__all__.append(name)
        wraplog.debug('  added interface routine %s' % routine)

    return pymods



def wrapmod(modobj, moddoc, modname, modfile, short_names, params, prefix, fortran_indexing, pymodname):

    wrapmethod = lambda name: lambda self, *args, **kwargs: self._runroutine(name, *args, **kwargs)
    wrapinterface = lambda name: lambda self, *args, **kwargs: self._runinterface(name, *args, **kwargs)
    def wrapinit(name, doc=None):
        init = lambda self, *args, **kwargs: FortranDerivedType.__init__(self, *args, **kwargs)
        init.__doc__ = doc
        return init

    routines = [x for x in moddoc['routines'].keys()]
    types =  [x for x in moddoc['types'].keys()]

    classes = []
    for cls in types:
        lcls = cls.lower()

        wraplog.debug('  Class %s' % cls)

        method_names = [ x for x in moddoc['routines'].keys()
                    if (len(moddoc['routines'][x]['args']) > 0 and
                        moddoc['routines'][x]['args'][0]['type'].lower() == 'type(%s)' % lcls) ]

        # Preferentially use initialise_ptr, finalise_ptr if available since this
        # does proper reference counting.
        constructors = \
            [ x for x in method_names if x.startswith('%s_initialise_ptr' % lcls) ] +\
            [ x for x in method_names if (x.startswith('%s_initialise' % lcls) or
                                     x.startswith('%s_init' % lcls) or
                                     x.startswith('%s_allocate' % lcls))]

        destructors = \
            [ x for x in method_names if x.startswith('%s_finalise_ptr' % lcls)]+\
            [ x for x in method_names if x.startswith('%s_finalise' % lcls)]


        if lcls in short_names:
            scls = short_names[lcls]
            constructors += \
                [ x for x in method_names if x.startswith('%s_initialise_ptr' % scls)
                  and 'intent(out)' in moddoc['routines'][x]['args'][0]['attributes'] ] + \
                [ x for x in method_names if(x.startswith('%s_initialise' % scls) or
                                        x.startswith('%s_allocate' % scls))
                  and 'intent(out)' in moddoc['routines'][x]['args'][0]['attributes'] ]

            destructors += \
                [ x for x in method_names if x.startswith('%s_finalise_ptr' % scls) ] + \
                [ x for x in method_names if x.startswith('%s_finalise' % scls) ]

        if (len(constructors) == 0):
            wraplog.debug("Can't find constructor for type %s. Skipping class" % cls)
            continue

        if (len(destructors) == 0):
            wraplog.debug("Can't find destructor for type %s. Skipping class" % cls)
            continue

        destructor = destructors[0]
        constructor = constructors[0]

        # Name of class, in TitleCase
        tcls = fortran_class_prefix+cls[0].upper()+cls[1:]
        tcls = tcls[len(fortran_class_prefix):]

        if constructor:
            constructor_name = constructor
            constructor = add_doc(wrapinit(constructor), getattr(modobj, prefix+constructor),
                                  moddoc['routines'][constructor], constructor, tcls, prefix,
                                  format='numpydoc', skip_this=True, modfile=modfile)

        _routines = {}
        _interfaces = {}
        methods = {}

        if constructor:
            _routines['__init__'] = (getattr(modobj, prefix+constructor_name), moddoc['routines'][constructor_name])
            
        if destructor:
            _routines['__del__'] = (getattr(modobj, prefix+destructor), moddoc['routines'][destructor])

        interfaces = {}

        for name in method_names:
            fullname = name
            if name[:len(lcls)+1] == lcls+'_':
                name = name[len(lcls)+1:]

            if lcls in short_names:
                if name[:len(short_names[lcls])+1] == short_names[lcls]+'_':
                    name = name[len(short_names[lcls])+1:]

            name = py_keywords_map.get(name, name)
            name = special_names.get(name, name)

            fobj = getattr(modobj, prefix+fullname)
            doc = moddoc['routines'][fullname]
            if fullname in constructors:
                name = '__init__'+name
                wraplog.debug('    adding constructor %s.%s' % (tcls, name))

            if fullname in destructors:
                del routines[routines.index(fullname)]
                continue

            _routines[name] = (fobj, doc)
            func = add_doc(wrapmethod(name), fobj, doc, fullname, tcls+'.'+name,
                           prefix, format='numpydoc', skip_this=True, modfile=modfile)
            func.__module__ = pymodname
            methods[name] = func

            for intf_name,value in moddoc['interfaces'].iteritems():
                intf_routines = value['routines']
                if fullname in intf_routines:
                    if not intf_name in interfaces:
                        interfaces[intf_name] = []

                    # regenerate interface documentation in sphinx format and with correct names
                    intf_func = add_doc(wrapmethod(name), fobj, doc, fullname, intf_name,
                                        prefix, format='sphinx', skip_this=True, modfile=modfile)
                    intf_func.__module__ = pymodname
                    interfaces[intf_name].append((name, moddoc['routines'][fullname], intf_func))

            wraplog.debug('    adding method %s' % name)
            del routines[routines.index(fullname)]

        # only keep interfaces with more than one routine in them
        for name,value in interfaces.iteritems():
            orig_name = name
            if name.lower() == 'finalise': continue
            if name.lower() == 'initialise' or \
                    name.lower() == 'initialise_ptr': name = '__init__'
            name = py_keywords_map.get(name, name)
            name = special_names.get(name, name)
            
            if len(value) > 1:
                _interfaces[name] = value
            else:
                # if there's just one routine we want to update the routine, prepending the interface docstring
                rname, rspec, rfunc = value[0]
                fobj, rspec = _routines[rname]
                rspec['doc'] = '\n'.join(moddoc['interfaces'][orig_name]['doc']) + '\n' + rspec['doc']

                del _routines[rname]
                _routines[name] = (fobj, rspec)
                if name == '__init__':
                    constructor = add_doc(wrapinit(name), fobj, rspec, rname, name, prefix,
                                          format='numpydoc', skip_this=True, modfile=modfile)
                    __init__ = constructor
                else:
                    func = add_doc(wrapmethod(name), fobj, rspec, rname, name, prefix,
                                   format='numpydoc', skip_this=True, modfile=modfile)
                    func.__module__ = pymodname
                    methods[name] = func                    

        for intf, value in _interfaces.iteritems():

            # if there's a routine with same name as interface, move it first
            if intf in methods:
                methods['__'+intf] = methods[intf]

            docname = intf
            if docname.endswith('_'): docname=docname[:-1]

            func = wrapinterface(intf)
            if intf == '__init__': docname = 'initialise'
            docname = rev_special_names.get(docname, docname)
            
            doc = '\n'.join(moddoc['interfaces'][docname]['doc']) + '\n\n'

            for rname,spec,routine in value:
                if spec['args_str']:
                    doc += args_str_table(spec)

            doc += 'Wrapper around Fortran interface ``%s`` containing multiple routines:\n\n' % intf

            for rname,spec,routine in value:
                if rname in methods:
                    methods['_'+rname] = methods[rname]
                    del methods[rname]

                routine_lines = routine.__doc__.split('\n')
                signature_line = routine_lines[0]
                if intf == '__init__':
                   signature_line = signature_line.replace('.'+rname,'')
                #signature_line = signature_line.replace(rname, intf)
                
                doc +=        ('   .. function :: %s\n' % signature_line +
                    '\n'.join(['      %s'   % line for line in routine_lines[1:]])) + '\n\n'

            if intf != '__init__':
                func.__doc__ = doc
                func.__module__ = pymodname
                methods[intf] = func
            else:
                constructor = wrapinit(constructor_name, doc)

        constructor_doc_lines = constructor.__doc__.split('\n')
        constructor_doc_lines.append('Class is wrapper around Fortran type ``%s`` defined in file :svn:`%s`.\n' % (cls, modfile))

        classdoc = '\n'.join([constructor_doc_lines[0]] + ['\n'] +
                             process_docstring(moddoc['types'][cls]['doc']).split('\n') +
                             ['\n'] + constructor_doc_lines[1:])

        new_cls = type(object)(tcls, (FortranDerivedType,),
                              {'__doc__': classdoc,
                               '_moddoc': None,
                               '_modobj': None,
                               '_classdoc': None,
                               '_routines': _routines,
                               '_subobjs': {},
                               '_arrays': {},
                               '_interfaces': _interfaces,
                               '_elements': {},
                               })
        FortranDerivedTypes['type(%s)' % cls.lower()] = new_cls
        classes.append((tcls, new_cls))

        new_cls._classdoc = moddoc['types'][cls]
        new_cls._moddoc = moddoc
        new_cls._modobj = modobj
        new_cls._prefix = prefix

        for name, func in methods.iteritems():
            setattr(new_cls, name, func)
        new_cls.__init__ = constructor

        # Add properties for get and set routines of scalar type, sub
        # objects, and arrays.

        for name, el in moddoc['types'][cls]['elements'].iteritems():
            if name == 'thetype': name = 'type' # special case: f2py misparses name 'type'
            name = name.lower() # lower case all attributes

            if 'get' in el and 'set' in el:
                if is_scalar_type(el['type']):
                    wraplog.debug('    adding scalar property %s get=%s set=%s' % (name, el['get'], el['set']))

                    new_cls._elements[name] = (getattr(modobj, el['get']), getattr(modobj, el['set']), el['type'])

                    # If element clashes with name of routine, move routine first
                    if hasattr(new_cls, name):
                        setattr(new_cls, name+'_', getattr(new_cls, name))

                    setattr(new_cls, name, property(fget=wrap_get(name),
                                                    fset=wrap_set(name),
                                                    doc=process_docstring(el['doc'])))
                #elif not 'pointer' in el['attributes']:
                else:
                    wraplog.debug('    adding property %s get=%s set=%s' % (name, el['get'], el['set']))
                    new_cls._subobjs[name] = (el['type'],
                                              getattr(modobj, el['get']),
                                              getattr(modobj, el['set']))
                    setattr(new_cls, name, property(fget=wrap_obj_get(name),
                                                    fset=wrap_obj_set(name),
                                                    doc=process_docstring(el['doc'])))


            # array members
            if 'array' in el:
                wraplog.debug('    adding array %s' % name)
                new_cls._arrays[name] = (getattr(modobj, el['array']), el['doc'], el['type'])
                setattr(new_cls, '_'+name, property(fget=wrap_array_get(name,reshape=False),
                                                    fset=wrap_array_set(name,reshape=False),
                                                    doc=process_docstring(el['doc'])))
                setattr(new_cls, name, property(fget=wrap_array_get(name),
                                                fset=wrap_array_set(name),
                                                doc=process_docstring(el['doc'])))

            # arrays of derived types (1D only for now..)
            if 'array_getitem' in el and 'array_setitem' and 'array_len' in el:
                wraplog.debug('    adding derived type array %s' % name)
                new_cls._arrays[name] = (getattr(modobj, el['array_getitem']),
                                         getattr(modobj, el['array_setitem']),
                                         getattr(modobj, el['array_len']),
                                         el['doc'], el['type'])

                setattr(new_cls, name, property(fget=wrap_derived_type_array_get(name),
                                                fset=wrap_derived_type_array_set(name),
                                                doc=process_docstring(el['doc'])))
                                                 




    # Try to evaluate params
    evaldict = np.__dict__
    evaldict['huge'] = lambda x: np.finfo('d').max
    evaldict['cmplx'] = lambda x,y,kind: complex(x,y)
    evaldict['farray'] = farray
    mod_params = {}
    for (name, ftype, attributes, value) in moddoc['parameters']:
        code = value.replace('(/','farray([')
        code = code.replace('/)','])')
        code = code.replace('_dp', '')
        try:
            mod_params[name] = eval(code, evaldict)
            evaldict[name] = mod_params[name]
            wraplog.debug('  adding parameter %s' % name)
        except NameError:
            wraplog.debug('  ignorning NameError in parameter %s = %s' % (name, code))

    return (classes, routines, mod_params)


def wrap_get(name):

    def func(self):
        if self._fpointer is None:
            raise ValueError('%s object not initialised.' % self.__class__.__name__)
        if not name in self._elements:
            raise ValueError('Unknown element %s in class %s' % (name, self.__class.__name))
        res = self._elements[name][0](self._fpointer)

        if self._elements[name][2] == 'logical':
            return res != QUIPPY_FALSE
        else:
            return res

    return func

def wrap_set(name):

    def func(self, value):
        if self._fpointer is None:
            raise ValueError('%s object not initialised.' % self.__class__.__name__)
        if not name in self._elements:
            raise ValueError('Unknown element %s in class %s' % (name, self.__class.__name))

        if self._elements[name][2] == 'logical':
            value = value and QUIPPY_TRUE or QUIPPY_FALSE

        res = self._elements[name][1](self._fpointer, value)
        self._update()
        return res

    return func

def wrap_obj_get(name):
    def func(self):
        if self._fpointer is None:
            raise ValueError('%s object not initialised.' % self.__class__.__name__)
        cls, getfunc, setfunc = self._subobjs[name]
        p = getfunc(self._fpointer)
        if not (p == 0).all():
            try:
                obj = self._subobjs_cache[tuple(p)]
            except KeyError:
                obj = self._subobjs_cache[tuple(p)] = FortranDerivedTypes[cls.lower()](fpointer=p,finalise=False,
                                                                                       fortran_indexing=self.fortran_indexing)
                obj.parent = weakref.proxy(self)
            return obj

    return func

def wrap_obj_set(name):
    def func(self, value):
        if self._fpointer is None:
            raise ValueError('%s object not initialised.' % self.__class__.__name__)
        cls, getfunc, setfunc = self._subobjs[name]
        setfunc(self._fpointer, value._fpointer)
        
    return func

def wrap_array_get(name, reshape=True):
    def func(self):
        try:
            arrayfunc, doc, arraytype = self._arrays[name]
        except KeyError:
            arrayfunc, doc, arraytype = self._arrays['_'+name]
        try:
            a = arraydata.get_array(self._fpointer, arrayfunc)
            if self.fortran_indexing:
                a = FortranArray(a, doc, parent=self)
            nshape = self._get_array_shape(name)
            if reshape and nshape is not None:
                a = a[nshape]
            return a
        except ValueError:
            return None

    return func

def wrap_array_set(name, reshape=True):
    def func(self, value):
        try:
            arrayfunc, doc, arraytype = self._arrays[name]
        except KeyError:
            arrayfunc, doc, arraytype = self._arrays['_'+name]
        try:
            a = arraydata.get_array(self._fpointer, arrayfunc)
            if self.fortran_indexing: a = FortranArray(a, doc)
            nshape = self._get_array_shape(name)
            if reshape and nshape is not None:
                a = a[nshape]
            if arraytype == 'logical':
                value = np.where(value, QUIPPY_TRUE, QUIPPY_FALSE)
            a[...] = value

        except ValueError:
            return None

    return func

def wrap_derived_type_array_get(name):
    def func(self):
        getfunc, setfunc, lenfunc, doc, arraytype = self._arrays[name]
        return FortranDerivedTypeArray(self, getfunc, setfunc,
                                       lenfunc, doc, arraytype,
                                       self.fortran_indexing)
    return func

def wrap_derived_type_array_set(name):
    def func(self, value):
        getfunc, setfunc, lenfunc, doc, arraytype = self._arrays[name]
        a = FortranDerivedTypeArray(self, getfunc, setfunc,
                                    lenfunc, doc, arraytype,
                                    self.fortran_indexing)
        for i,v in zip(a.indices, value):
            a[i] = v
        
    return func        

class FortranDerivedTypeArray(object):
    def __init__(self, parent, getfunc, setfunc, lenfunc, doc, arraytype, fortran_indexing):
        self.parent = weakref.proxy(parent)
        self.getfunc = getfunc
        self.setfunc = setfunc
        self.lenfunc = lenfunc
        self.doc = doc
        self.arraytype = arraytype
        self.fortran_indexing = fortran_indexing

    def iterindices(self):
        if self.fortran_indexing:
            return iter(frange(len(self)))
        else:
            return iter(range(len(self)))

    indices = property(iterindices)

    def iteritems(self):
        for idx in self.indices:
            yield self[idx]

    def __iter__(self):
        return self.iteritems()

    def __len__(self):
        return self.lenfunc(self.parent._fpointer)

    def __getitem__(self, i):
        if not self.fortran_indexing:
            i += 1
        pp = self.getfunc(self.parent._fpointer, i)
        try:
            obj = p._subobjs_cache[tuple(pp)]
        except KeyError:
            obj = p._subobjs_cache[tuple(pp)] = FortranDerivedTypes[self.arraytype.lower()](fpointer=pp,finalise=False,
                                                                                            fortran_indexing=self.fortran_indexing)
        return obj

    def __setitem__(self, i, value):
        if not self.fortran_indexing:
            i += 1
        self.setfunc(self.parent._fpointer, i, value._fpointer)


# regexp magic to convert old f90doc documentation syntax to something more compatible with ResT
doc_subs = [(re.compile(r'([A-Za-z0-9_]+)%([A-Za-z0-9_]+)'), r'\1.\2'), # Fortran attribute access % -> Python .
            (re.compile(r"'(.*?)'"), r"``\1``"), # Single quoted strings to fixed width font
            (re.compile(r"(\w+)_(\W)"), r"\1\_\2"), # Escape underscores at the end of words...
            (re.compile(r"\$(.*?)(\\_)(.*?)\$"), r'$\1_\3$'), # ... except in math mode
            (re.compile(r"\$(.*?)\$"), r":math:`\1`") # LaTeX $...$ -> :math:`...`
            ]

def normalise_type_case(type):
    classnames = [c.__name__ for c in FortranDerivedTypes.values()]
    lower_to_normed = dict((c.lower(), c) for c in classnames)
    return lower_to_normed.get(type.lower(), type)

def process_docstring(doc):
    # turn verbatim and displayed maths sections into Sphinx blockquotes
    indent = 0
    lines = doc.split('\n')
    in_verbatim = False
    in_math = False
    in_literal = False
    first_line_in_verbatim = False
    for i in range(len(lines)):
        if lines[i].startswith('>'):
            if not in_verbatim:
                in_verbatim = True
                first_line_in_verbatim = True
                indent += 3
                if lines[i-1] != '':
                    if lines[i-1].endswith(':'):
                        lines[i-1] += ':\n'
                    else:
                        lines[i-1] += ' ::\n'
                elif not lines[i-2].endswith('::'):
                    if lines[i-2].endswith(':'):
                        lines[i-2] += ':'
                    else:
                        lines[i-2] += ' ::'
            lines[i] = lines[i].replace('>', '')
        elif in_verbatim:
            if first_line_in_verbatim:
                first_line_in_verbatim = False
                if lines[i-1] != '':
                    lines[i-1] += '\n'
            in_verbatim = False
            indent -= 3

        if i > 1 and lines[i-1].endswith('::') and not in_literal:
            literal_start_col = len(lines[i-1]) - len(lines[i-1].lstrip())
            in_literal = True

        if in_literal and lines[i]:
            if len(lines[i]) > literal_start_col and lines[i][literal_start_col] != ' ':
                in_literal = False

        if '\\begin{displaymath}' in lines[i]:
            lines[i] = lines[i].replace('\\begin{displaymath}', '\n.. math::\n')
            indent += 3
            in_math = True

        if '\\end{displaymath}' in lines[i]:
            lines[i] = lines[i].replace('\\end{displaymath}', '')
            indent -= 3
            in_math = False

        if not in_verbatim and not in_literal and not in_math:
            for regexp, repl in doc_subs:
                lines[i] = regexp.sub(repl, lines[i])

        lines[i] = ' '*indent + lines[i]
        
    try:
       doc = inspect.cleandoc('\n'.join(lines))
    except AttributeError:
       doc = '\n'.join(lines)

    return doc


def add_doc(func, fobj, doc, fullname, name, prefix, format='numpydoc', skip_this=False, modfile=None):
    func.__name__ = fullname
    func._fobj = fobj

    if doc is None:
        func.__doc__ = None
        return func

    if format not in ['numpydoc', 'sphinx']:
        raise ValueError('Unsuported format %s' % format)

    if format == 'numpydoc':
        arg_header = 'Parameters\n----------\n'
        ret_header = 'Returns\n-------\n'
        ref_header = 'References\n----------\n'
    else:
        arg_header = ''
        ret_header = ''
        ref_header = ''

    d = fobj.__doc__
    L = d.split('\n')
    arg_lines_in = L[3:]

    arg_lines = []
    ret_lines = []

    signature_line = L[1]

    if '.' in name:
        signature_line = signature_line.replace(name[:name.index('.')].lower()+'_', '')
    signature_line = signature_line.replace(prefix, '')

    if '=' in signature_line:
        signature_line = signature_line[:signature_line.index('=')] + signature_line[signature_line.index('='):].replace(fullname, name)
        retvar = signature_line[:signature_line.index('=')].strip()
        signature_line = signature_line[signature_line.index('=')+1:].strip()
    else:
        signature_line = signature_line.replace(fullname.lower()+'_', '')
        retvar = None
        signature_line = signature_line.replace(fullname, name)

    if skip_this:
        signature_line = signature_line.replace('this,', '')
        signature_line = signature_line.replace('this[,', '[')
        signature_line = signature_line.replace('(this)', '')


    final_doc = signature_line + '\n\n'
    
    if doc['doc']:
        final_doc += process_docstring(doc['doc']) + '\n\n'

    got_args_str = False
    if arg_lines_in:
        for arg in doc['args']:
            argname = arg['name'].lower()
            if argname == 'args_str':
                got_args_str = True
            if argname in badnames: argname = badnames[argname]

            for i, line in enumerate(arg_lines_in):
                if line.startswith('  %s :' % argname): break
            else:
                raise ValueError('%s not found in lines %r' % (argname, arg_lines_in))

            line = line.replace(argname,argname[len(prefix):])
            argname, type = line.split(':', 2)
            argname = argname.strip()
            if argname == 'this' and skip_this:
                continue

            if arg['type'].startswith('type('):
                type = arg['type'][5:-1]
                type = normalise_type_case(type)
                type = ':class:`~.%s` object' % type
            if 'optional' in arg['attributes']:
                type += ', optional'
                

            if argname == retvar:
                if format == 'numpydoc':
                    ret_lines.append('%s : %s' % (argname, type))
                    if arg['doc']:
                        ret_lines.extend(['    ' + line for line in process_docstring(arg['doc']).split('\n') ])
                else:
                    ret_lines.append(':returns: **%s** -- %s %s' % (argname, type, arg['doc']))
            else:
                if format == 'numpydoc':
                    arg_lines.append('%s : %s' % (argname, type))
                    if arg['doc']:
                        arg_lines.extend(['    ' + line for line in process_docstring(arg['doc']).split('\n') ])
                else:
                    arg_lines.append(':param %s: %s' % (argname, '\n'.join(process_docstring(arg['doc']).split('\n'))))
                    if type:
                        arg_lines.append(':type %s: %s' % (argname, type))

    if format != 'sphinx':
        final_doc += args_str_table(doc)

    if arg_lines:
        final_doc += arg_header + '\n'.join(arg_lines) + '\n\n'
    if ret_lines:
        final_doc += ret_header + '\n'.join(ret_lines) + '\n\n'

    final_doc += ref_header + '\nRoutine is wrapper around Fortran routine ``%s`` defined in file :svn:`%s`.' % (fullname, modfile)        
            
    func.__doc__ = final_doc
    return func


def args_str_table(spec):
    if not spec['args_str']:
        return ''

    value_map = {'T': 'True',
                 'F': 'False'}

    names = ['Name']
    types = ['Type']
    defaults = ['Default']
    docs = ['Comments']

    for name in sorted(spec['args_str'].keys()):
        arg = spec['args_str'][name]
        default = 'None'
        if arg['value']:
            default = value_map.get(arg['value'], arg['value'])
        type = python_equivalent_type(arg['type'])
        doc = arg['doc']

        names.append(name.strip())
        types.append(type.strip())
        defaults.append(default.strip())
        doc = doc.strip()
        doc = doc[0].capitalize() + doc[1:] # capitalise first letter only
        docs.append(doc)

    max_name_len = max(len(name) for name in names)
    max_type_len = max(len(type) for type in types)
    max_default_len = max(len(default) for default in defaults)
    
    cols = (max_name_len, max_type_len, max_default_len, 40)
    
    args_str_lines = ['.. rubric:: args_str options','']
    fmt = '%%-%ds %%-%ds %%-%ds %%-%ds' % cols

    for i, (name, type, default, doc) in enumerate(zip(names, types, defaults, docs)):
        if i == 0:
            args_str_lines.append(fmt % ('='*cols[0], '='*cols[1], '='*cols[2], '='*cols[3]))            
        doc_words = doc.split()
        while doc_words:
            doc_line = ''
            while doc_words:
                word = doc_words.pop(0)
                if len(doc_line) + 1 + len(word) > cols[3]:
                    doc_words.insert(0, word) # put it back
                    break
                else:
                    doc_line = doc_line + ' ' + word

            args_str_lines.append(fmt % (name, type, default, doc_line.strip()))
            name = type = default = ''
        if i == 0 or i == len(names)-1:
            args_str_lines.append(fmt % ('='*cols[0], '='*cols[1], '='*cols[2], '='*cols[3]))            

    args_str_lines.extend(['', ''])
    
    return '\n'.join(args_str_lines)



def wraproutine(modobj, moddoc, name, shortname, prefix, fortran_indexing=True, skip_this=False, modfile=None):
    doc = moddoc['routines'][name]
    fobj = getattr(modobj, prefix+name)

    inargs  = [ x for x in doc['args'] if not 'intent(out)' in x['attributes'] ]
    outargs = [ x for x in doc['args'] if 'intent(out)' in x['attributes'] ]        

    #inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'])
    #outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])

    def func(*args, **kwargs):
        newargs, modargs, newkwargs, modkwargs = process_in_args(args, kwargs, inargs, prefix)

        try:
            res = fobj(*newargs, **newkwargs)
        except RuntimeError:
            try:
                exctype, value, tb = sys.exc_info()

                if is_interactive_shell():
                    error_str = 'Fortran routine %s(%s):\n %s' % (name, ', '.join([repr(a) for a in args] + ['%s=%r' % (k,v) for (k,v) in kwargs.iteritems()]), str(value).strip())
                else:
                    # Remove Fortran traceback
                    error_str = str(value).strip().split('\n')[-1].strip()

                raise exctype, RuntimeError(error_str), tb

            finally:
                del tb
        except:
            raise

        newres = process_results(res, args, kwargs, inargs, outargs, prefix, fortran_indexing)

        return newres

    return add_doc(func, fobj, doc, name, shortname, prefix,
                   format='numpydoc', skip_this=skip_this, modfile=modfile)


def wrapinterface(name, intf_spec, routines, prefix, modfile=None):

    def func(*args, **kwargs):
        wraplog.debug('Interface %s ' % name)

        for rname, spec, routine, modfile in routines:
            wraplog.debug('Trying candidate routine %s' % rname)

            inargs  = [ x for x in spec['args'] if not 'intent(out)' in x['attributes'] ]
            oblig_inargs = [x for x in inargs if  not 'optional' in x['attributes'] ]
            opt_inargs   = [x for x in inargs if 'optional' in x['attributes'] ]

            #inargs = filter(lambda x: not 'intent(out)' in x['attributes'], spec['args'])
            #oblig_inargs = filter(lambda x: not 'optional' in x['attributes'], inargs)
            #opt_inargs = filter(lambda x: 'optional' in x['attributes'], inargs)

            # Check number of arguments is compatible
            # Should be in range len(oblig_inargs) <= L <= len(oblig_inargs) + len(opt_inargs)
            if (len(args)+len(kwargs) < len(oblig_inargs) or
                len(args)+len(kwargs) > len(oblig_inargs) + len(opt_inargs)):
                if wraplog.getEffectiveLevel() <= logging.DEBUG:
                    wraplog.debug('Number of arguments incompatible: %d must be in range %d <= n <= %d' %
                                  (len(args), len(oblig_inargs), len(oblig_inargs)+len(opt_inargs)))
                continue

            newinargs = (oblig_inargs + opt_inargs)[:len(args)]

            # Check types and dimensions are compatible
            if not all([type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]):
                if wraplog.getEffectiveLevel() <= logging.DEBUG:
                    wraplog.debug('Types and dimensions of args incompatible %s %s %s' %
                                  (newinargs, args, [type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]))
                continue

            # Check keyword arguments, if incompatible continue to next
            # candidate interface
            if kwargs:
                innames = [badnames.get(x['name'],x['name']).lower() for x in oblig_inargs + opt_inargs]

            if not all([prefix+key.lower() in innames[len(args):] for key in kwargs.keys()]):
                if wraplog.getEffectiveLevel() <= logging.DEBUG:
                    wraplog.debug('Unexpected keyword argument valid=%s got=%s' % (kwargs.keys(), innames[len(args):]))
                continue

            try:
                for key, arg in kwargs.iteritems():
                    (inarg,) = [x for x in inargs if badnames.get(x['name'],x['name']).lower() == prefix+key.lower()]
                    if not type_is_compatible(inarg, arg):
                        if wraplog.getEffectiveLevel() <= logging.DEBUG:
                            wraplog.debug('Types and dimensions of kwarg %s incompatible' % key)
                        raise ValueError
            except ValueError:
                continue

            wraplog.debug('calling '+rname)
            return routine(*args, **kwargs)

        raise TypeError('No matching routine found in interface %s' % name)


    def func_simple(*args, **kwargs):
        wraplog.debug('Interface %s' % name)

        for rname, spec, routine, modfile in routines:
            wraplog.debug('Trying candidate routine %s' % rname)
            try:
                return routine(self, *args, **kwargs)
            except Exception, e:
                if isinstance(e, RuntimeError):
                    raise
                else:
                    continue

        raise TypeError('No matching routine found in interface %s' % name)


    doc = '\n'.join(intf_spec['doc']) + '\n\n'

    doc += 'Routine is wrapper around Fortran interface ``%s`` containing multiple routines:\n\n' % name
    
    for rname, spec, routine, modfile in routines:

        # regenerate routine documentation, using sphinx format rather than numpydoc format
        tmp_routine = add_doc(routine, routine._fobj, spec, rname, name, prefix, format='sphinx', modfile=modfile)
        
        routine_lines = tmp_routine.__doc__.split('\n')
        doc +=        ('  .. function :: %s\n' % routine_lines[0] +
            '\n'.join(['     %s'   % line for line in routine_lines[1:]])) + '\n'

        #doc = doc.replace("``%s``"%rname,"``%s``"%rname.upper())
        #doc = doc.replace(rname,name)
        #doc = doc.replace("``%s``"%rname.upper(),"``%s``"%rname)

    func.__doc__ = doc
    return func


def update_doc_string(doc, extra, sections=None, signature=None):
    """
    Insert `extra` in the docstring `doc`, before the first matching section

    Searches for each section heading in the list `sections` in turn.
    If sections is not given, the default is `['Parameters', 'See also']`.
    If not sections are found, extra text is appended to end of docstring.

    Returns the new doc string, suitable for assigning to cls.__doc__
    """

    if sections is None:
       sections = ['Parameters', 'See also']

    try:
       doc = inspect.cleandoc(doc)
       extra = inspect.cleandoc(extra)
    except AttributeError:
       pass

    extra = '\n' + extra + '\n'
    
    lines = doc.split('\n')

    if signature is not None:
        lines[0] = signature
    
    for section in sections:
       indices = [i for i, line in enumerate(lines) if line == section ]
       if len(indices) == 1:
          break
    else:
        indices = [len(lines)-1] # insert at end

    index, = indices
    doc = '\n'.join([line.rstrip() for line in lines[:index] + extra.split('\n') + lines[index:]])
    doc = doc.replace('\n\n\n', '\n\n')

    return doc
