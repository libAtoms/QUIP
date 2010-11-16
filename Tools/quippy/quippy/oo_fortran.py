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

"""Object oriented interface on top of f2py generated wrappers."""

import logging, numpy, sys
import arraydata
from farray import *
import numpy
from types import MethodType
from util import args_str, is_interactive_shell
from dictmixin import PuPyDictionary

major, minor = sys.version_info[0:2]

if (major, minor) < (2, 5):
    all = lambda seq: not False in seq
    any = lambda seq: True in seq


fortran_class_prefix = 'Fortran'

py_keywords = ['and',       'del',       'from',      'not',       'while',    
               'as',        'elif',      'global',    'or',        'with',     
               'assert',    'else',      'if',        'pass',      'yield',    
               'break',     'except',    'import',    'print',              
               'class',     'exec',      'in',        'raise',              
               'continue',  'finally',   'is',        'return',             
               'def',       'for',       'lambda',    'try']

numpy_to_fortran = {
    'd': 'real(dp)',
    'i': 'integer',
    'b': 'logical',
    'c': 'complex(dp)',
    'S': 'character',
    'f': 'real(dp)'
    }

FortranDerivedTypes = {}

from quippy import QUIPPY_TRUE, QUIPPY_FALSE

from numpy.f2py.crackfortran import badnames

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
      return ('type(%s)' % t.__class__.__name__.lower(), 'scalar', 1)
   elif hasattr(t,'__iter__') or hasattr(t,'dtype'):
      a = numpy.array(t)
      if a.dtype.kind in numpy_to_fortran:
         return (numpy_to_fortran[a.dtype.kind], a.shape, a.dtype.itemsize)
      else:
         raise TypeError('Unknown array type %s' % a.dtype.kind)      
   else:
      raise TypeError('Unknown type %s' % type(t))
          
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
   for arg, spec in zip(args,inargs):
      if arg is not None and spec['type'].startswith('type'):
         if isinstance(arg, FortranDerivedTypes[spec['type'].lower()]):
            newargs.append(arg._fpointer)
         else:
             raise TypeError('Argument %s should be of type %s' % (arg, spec['type']))
      elif arg is not None and spec['type'].startswith('character'):
          # if arg is a list of strings of unequal length, pad with spaces
          if isinstance(arg, list) and any([len(x) != len(arg[0]) for x in arg]):
              newargs.append(s2a(arg).T)
          else:
              newargs.append(arg)
      else:
         newargs.append(arg)

   # Process keyword arguments and args_str arguments.

   # If we're expecting an args_str argument, then any kwarg is
   # permitted. Otherwise kwarg name must match one of the expected
   # names. If kwarg name matches one of expected values, but type of
   # argument does not, then we assume it's an arg_str argument.

   kwarg_lookup = dict([(badnames.get(x['name'].lower(),x['name'].lower()),x) for x in inargs])
   got_args_str = prefix+'args_str' in kwarg_lookup

   newkwargs = {}
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
         if isinstance(a, FortranDerivedTypes[kwarg_lookup[k]['type'].lower()]):
            newkwargs[k] = a._fpointer
         elif a is None:
             continue
         else:
             raise TypeError('Argument %s=%s should be of type %s' % (k,a,kwarg_lookup[k]))
      else:
          if a is None: continue
          newkwargs[k] = a

   # Construct final args_str by merging args_str argument with args_str_kwargs
   if got_args_str:
       args_str_final = {}
       if prefix+'args_str' in newkwargs:
           args_str_final.update(PuPyDictionary(newkwargs[prefix+'args_str']))
       if args_str_kwargs != {}:
           args_str_final.update(args_str_kwargs)
       if args_str_final != {}:
           newkwargs[prefix+'args_str'] = args_str(args_str_final)

   return tuple(newargs), newkwargs

def process_results(res, args, kwargs, inargs, outargs, prefix):
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
          elif isinstance(r, numpy.ndarray):
              # Convert to one-based FortranArray
              newres.append(r.view(FortranArray))
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
   """Abstract base class for all fortran derived types."""

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

   _prefix = ''

   def __init__(self, *args, **kwargs):

       call_init = True
       self._fpointer = None
       self._finalise = True
       self._subobjs_cache = {}
       
       if 'fpointer' in kwargs:
           if kwargs['fpointer'] is not None:
               call_init = False
               self._fpointer = kwargs['fpointer']
           del kwargs['fpointer']

       if 'finalise' in kwargs:
           self._finalise = kwargs['finalise']
           del kwargs['finalise']

       orig_args = args
       orig_kwargs = kwargs.copy()

       if call_init:
           if '__init__' in self._interfaces:
               self._fpointer = self._runinterface('__init__', *args, **kwargs)
           elif '__init__' in self._routines:
               self._fpointer = self._runroutine('__init__', *args, **kwargs)

       self._update()

   def __del__(self):
       #print 'del %s(fpointer=%s, finalise=%d)' % (self.__class__.__name__, self._fpointer, self._finalise)

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
      import logging

      # Class name must match
      if self.__class__ != other.__class__:
          logging.debug('class mismatch %s != %s')
          return False

      # Don't compare uninitialised objects
      if hasattr(self, 'initialised') and not self.initialised:
          return True

      # Compare elements and sub-objects
      for el in self._elements.keys() + self._subobjs.keys():
          if el in self._cmp_skip_fields:
              continue
          
          if getattr(self, el) != getattr(other, el):
              logging.debug('element mismatch %s' % el)
              return False

      # Compare arrays
      for array in self._arrays:
          a = getattr(self, array)
          b = getattr(other, array)
          if a is None and b is None: continue
          if a is None or b is None: return False
          if a.dtype.kind != 'f':
              if not (a == b).all():
                  logging.debug('array mismatch %s' % array)
                  return False
          else:
              if (abs(a-b) > self._cmp_tol).any():
                  logging.debug('real array mismatch %s' % array)
                  return False
      return True

   def __ne__(self, other):
      return not self.__eq__(other)

   def _update(self):
       logging.debug('updating %s at 0x%x' % (self.__class__, id(self)))
       if self._fpointer is None: return
       self._update_hook()

   def _get_array_shape(self, name):
       return None

   def _update_hook(self):
       pass

   def _runroutine(self, name, *args, **kwargs):
       if not name.startswith('__init__') and self._fpointer is None:
           raise ValueError('%s object not initialised.' % self.__class__.__name__)

       if not name in self._routines:
           raise NameError('Unknown fortran routine: %s' % name)

       fobj, doc = self._routines[name]

       inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'])
       outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])

       if not name.startswith('__init__'):
           # Put self at beginning of args list
           args = tuple([self] + list(args))
       newargs, newkwargs = process_in_args(args, kwargs, inargs, self._prefix)
       
       try:
           res = fobj(*newargs, **newkwargs)
       except RuntimeError:
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
           newres = process_results(res, args, kwargs, inargs, outargs, self._prefix)
       else:
           newres = res
      
       return newres


   def _runinterface(self, name, *args, **kwargs):

       if not name in self._interfaces:
           raise ValueError('Unknown interface %s.%s' % (self.__class__.__name__, name))
       
       logging.debug('Interface %s %s' % (self.__class__.__name__, name))

       for rname, spec, routine in self._interfaces[name]:

           logging.debug('Trying candidate routine %s' % rname)

           inargs = filter(lambda x: not 'intent(out)' in x['attributes'], spec['args'])

           if not name.startswith('__init__'):
               inargs = inargs[1:] # remove self argument
           
           oblig_inargs = filter(lambda x: not 'optional' in x['attributes'], inargs)
           opt_inargs = filter(lambda x: 'optional' in x['attributes'], inargs)

           # Check number of arguments is compatible
           # Should be in range len(oblig_inargs) <= L <= len(oblig_inargs) + len(opt_inargs)
           if (len(args)+len(kwargs) < len(oblig_inargs) or
               len(args)+len(kwargs) > len(oblig_inargs) + len(opt_inargs)):
               logging.debug('Number of arguments incompatible: %d must be in range %d <= n <= %d' %
                             (len(args), len(oblig_inargs), len(oblig_inargs)+len(opt_inargs)))
               continue

           newinargs = (oblig_inargs + opt_inargs)[:len(args)]

           # Check types and dimensions are compatible
           if not all([type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]):
               logging.debug('Types and dimensions of args incompatible %s %s %s' %
                             (newinargs, args, [type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]))
               continue

           # Check keyword arguments, if incompatible continue to next
           # candidate interface
           if kwargs:
               innames = [badnames.get(x['name'],x['name']).lower() for x in oblig_inargs + opt_inargs]

           if not all([self._prefix+key.lower() in innames[len(args):] for key in kwargs.keys()]):
               logging.debug('Unexpected keyword argument valid=%s got=%s' % (kwargs.keys(), innames[len(args):]))
               continue

           try:
               for key, arg in kwargs.iteritems():
                   (inarg,) = [x for x in inargs if badnames.get(x['name'],x['name']).lower() == self._prefix+key.lower()]
                   if not type_is_compatible(inarg, arg):
                       logging.debug('Types and dimensions of kwarg %s incompatible' % key)
                       raise ValueError
                   
           except ValueError:
               continue

           logging.debug('calling '+rname)
           return routine(self, *args, **kwargs)


       raise TypeError('No matching routine found in interface %s' % name)


def wrap_all(fobj, spec, mods, short_names, prefix):
   all_classes = []
   leftover_routines = []
   standalone_routines = []
   params = {}

   for mod in mods:
       logging.debug('Module %s' % mod)
       classes, routines, params = wrapmod(fobj, spec[mod],
                                           short_names=short_names,
                                           params=params, prefix=prefix)
       all_classes.extend(classes)
       leftover_routines.append((spec[mod], routines))

   # Now try to add orphaned routines to classes where possible
   for modspec, routines in leftover_routines:
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

            if method_name in py_keywords: name = name+'_'
         
            setattr(cls, method_name, wraproutine(fobj, modspec, routine, cls.__name__+'.'+method_name, prefix))
            logging.debug('  added method %s to class %s' % (method_name, cls.__name__))
                  
         else:
            standalone_routines.append((routine,wraproutine(fobj, modspec, routine, routine, prefix)))
            logging.debug('  added stand-alone routine %s' % routine)

   return all_classes, standalone_routines, params



def wrapmod(modobj, moddoc, short_names, params, prefix):

   wrapmethod = lambda name: lambda self, *args, **kwargs: self._runroutine(name, *args, **kwargs)
   wrapinterface = lambda name: lambda self, *args, **kwargs: self._runinterface(name, *args, **kwargs)
   wrapinit = lambda name: lambda self, *args, **kwargs: FortranDerivedType.__init__(self, *args, **kwargs)

   routines = [x for x in moddoc['routines'].keys()]
   types =  [x for x in moddoc['types'].keys()]

   classes = []
   for cls in types:
      lcls = cls.lower()

      logging.debug('  Class %s' % cls)

      methods = [ x for x in moddoc['routines'].keys()
                  if (len(moddoc['routines'][x]['args']) > 0 and 
                      moddoc['routines'][x]['args'][0]['type'].lower() == 'type(%s)' % lcls) ]

      # Preferentially use initialise_ptr, finalise_ptr if available since this
      # does proper reference counting.
      constructors = \
          [ x for x in methods if x.startswith('%s_initialise_ptr' % lcls) ] +\
          [ x for x in methods if (x.startswith('%s_initialise' % lcls) or
                                   x.startswith('%s_allocate' % lcls))]

      destructors = \
          [ x for x in methods if x.startswith('%s_finalise_ptr' % lcls)]+\
          [ x for x in methods if x.startswith('%s_finalise' % lcls)]


      if lcls in short_names:
          scls = short_names[lcls]
          constructors += \
              [ x for x in methods if x.startswith('%s_initialise_ptr' % scls)
                and 'intent(out)' in moddoc['routines'][x]['args'][0]['attributes'] ] + \
              [ x for x in methods if(x.startswith('%s_initialise' % scls) or
                                      x.startswith('%s_allocate' % scls))
                and 'intent(out)' in moddoc['routines'][x]['args'][0]['attributes'] ]

          destructors += \
              [ x for x in methods if x.startswith('%s_finalise_ptr' % scls) ] + \
              [ x for x in methods if x.startswith('%s_finalise' % scls) ]

      if (len(constructors) == 0):
          logging.debug("Can't find constructor for type %s. Skipping class" % cls)
          continue
          
      if (len(destructors) == 0):
          logging.debug("Can't find destructor for type %s. Skipping class" % cls)
          continue

      destructor = destructors[0]
      constructor = constructors[0]

      tcls = fortran_class_prefix+cls[0].upper()+cls[1:]
      new_cls = type(object)(tcls, (FortranDerivedType,),
                            {'__doc__': moddoc['doc'],
                             '_moddoc': None,
                             '_modobj': None,
                             '_classdoc': None,
                             '_routines': {},
                             '_subobjs': {},
                             '_arrays': {},
                             '_interfaces': {},
                             '_elements': {},
                             })
      FortranDerivedTypes['type(%s)' % cls.lower()] = new_cls
      classes.append((tcls, new_cls))

      tcls = tcls[len(fortran_class_prefix):]

      new_cls._classdoc = moddoc['types'][cls]
      new_cls._moddoc = moddoc
      new_cls._modobj = modobj
      new_cls._prefix = prefix

      if constructor:
          new_cls._routines['__init__'] = (getattr(modobj, prefix+constructor), moddoc['routines'][constructor])
          new_cls.__init__ = add_doc(wrapinit(constructor), getattr(modobj, prefix+constructor),
                                     moddoc['routines'][constructor], constructor, tcls, prefix)

      if destructor:
          new_cls._routines['__del__'] = (getattr(modobj, prefix+destructor), moddoc['routines'][destructor])

      interfaces = {}

      for name in methods:
         fullname = name
         if name[:len(lcls)+1] == lcls+'_':
            name = name[len(lcls)+1:]

         if lcls in short_names:
             if name[:len(short_names[lcls])+1] == short_names[lcls]+'_':
                name = name[len(short_names[lcls])+1:]

         if name in py_keywords: name = name+'_'

         fobj = getattr(modobj, prefix+fullname)
         doc = moddoc['routines'][fullname]
         if fullname in constructors:
             name = '__init__'+name
             logging.debug('    adding constructor %s.%s' % (tcls, name))

         if fullname in destructors:
             del routines[routines.index(fullname)]
             continue
         
         new_cls._routines[name] = (fobj, doc)
         func = add_doc(wrapmethod(name), fobj, doc, name, tcls+'.'+name, prefix)
         setattr(new_cls, name, func)
         
         for intf_name,value in moddoc['interfaces'].iteritems():
            intf_routines = value['routines']
            if fullname in intf_routines:
               if not intf_name in interfaces: interfaces[intf_name] = []
               interfaces[intf_name].append((name, moddoc['routines'][fullname], getattr(new_cls,name)))

         logging.debug('    adding method %s' % name)
         del routines[routines.index(fullname)]

      # only keep interfaces with more than one routine in them
      # if there's just one routine we want to copy the docstring
      has_initialise_ptr = 'initialise_ptr' in interfaces.keys()

      for name,value in interfaces.iteritems():
         if name.lower() == 'finalise': continue
         if name.lower() == 'initialise' or \
                 name.lower() == 'initialise_ptr': name = '__init__'
         if name in py_keywords: name = name + '_'

         if len(value) > 1: 
             new_cls._interfaces[name] = value

      for intf, value in new_cls._interfaces.iteritems():

         # if there's a routine with same name as interface, move it first
         if hasattr(new_cls, intf):
            setattr(new_cls,intf+'_', getattr(new_cls, intf))

         docname = intf
         if docname.endswith('_'): docname=docname[:-1]

         func = wrapinterface(intf)
         if intf == '__init__': docname = 'initialise'
         doc = '\n'.join(moddoc['interfaces'][docname]['doc']) + '\n\n'
         
         for rname,spec,routine in value:
             if hasattr(new_cls, rname):
                 setattr(new_cls, '_'+rname, getattr(new_cls, rname))
                 delattr(new_cls, rname)
             
             if intf == '__init__':
                 doc += routine.__doc__.replace('.'+rname,'')
             else:
                 doc += routine.__doc__.replace(rname,intf)

         if intf != '__init__':
             func.__doc__ = doc
             setattr(new_cls, intf, func)
         else:
             func  = wrapinit(constructor)
             func.__doc__ = doc
             new_cls.__init__ = func

      # Add properties for get and set routines of scalar type, sub
      # objects, and arrays.

      for name, el in moddoc['types'][cls]['elements'].iteritems():
         if name == 'thetype': name = 'type' # special case: f2py misparses name 'type'
         name = name.lower() # lower case all attributes

         if 'get' in el and 'set' in el:
            if is_scalar_type(el['type']):
               logging.debug('    adding scalar property %s get=%s set=%s' % (name, el['get'], el['set']))

               new_cls._elements[name] = (getattr(modobj, el['get']), getattr(modobj, el['set']), el['type'])
                  
               # If element clashes with name of routine, move routine first
               if hasattr(new_cls, name):
                   setattr(new_cls, name+'_', getattr(new_cls, name))
                   
               setattr(new_cls, name, property(fget=wrap_get(name),
                                               fset=wrap_set(name),
                                               doc=el['doc']))
            #elif not 'pointer' in el['attributes']:
            else:
                logging.debug('    adding property %s get=%s set=%s' % (name, el['get'], el['set']))
                new_cls._subobjs[name] = (el['type'], 
                                          getattr(modobj, el['get']), 
                                          getattr(modobj, el['set']))
                setattr(new_cls, name, property(fget=wrap_obj_get(name),
                                                fset=wrap_obj_set(name),
                                                doc=el['doc']))
                        

         # array members
         if 'array' in el:
            logging.debug('    adding array %s' % name)
            new_cls._arrays[name] = (getattr(modobj, el['array']), el['doc'])
            setattr(new_cls, '_'+name, property(fget=wrap_array_get(name,reshape=False),doc=el['doc']))
            setattr(new_cls, name, property(fget=wrap_array_get(name),doc=el['doc']))

            
   

   # Try to evaluate params
   evaldict = numpy.__dict__
   evaldict['huge'] = lambda x: numpy.finfo('d').max
   evaldict['cmplx'] = lambda x,y,kind: complex(x,y)
   evaldict['farray'] = farray
   for (name, ftype, attributes, value) in moddoc['parameters']:
       code = value.replace('(/','farray([')
       code = code.replace('/)','])')
       code = code.replace('_dp', '')
       try:
           params[name] = eval(code, evaldict)
           evaldict[name] = params[name]
           logging.debug('  adding parameter %s' % name)
       except NameError:
           logging.debug('  ignorning NameError in parameter %s = %s' % (name, code))

   
   return (classes, routines, params)


def wrap_get(name):
                               
   def func(self):
       if self._fpointer is None:
           raise ValueError('%s object not initialised.' % self.__class__.__name__)
       if not name in self._elements:
           raise ValueError('Unknown element %s in class %s' % (name, self.__class.__name)) 
       res = self._elements[name][0](self._fpointer)

       if self._elements[name][2] == 'logical':
           return res == QUIPPY_TRUE
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
             obj = self._subobjs_cache[tuple(p)] = FortranDerivedTypes[cls.lower()](fpointer=p,finalise=False)
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
            arrayfunc, doc = self._arrays[name]
        except KeyError:
            arrayfunc, doc = self._arrays['_'+name]
        try:
            a = arraydata.get_array(self._fpointer, arrayfunc)
            nshape = self._get_array_shape(name)
            if reshape and nshape is not None:
                return FortranArray(a,doc)[nshape]
            else:
                return FortranArray(a,doc)
        except ValueError:
            return None
            
    return func


def add_doc(func, fobj, doc, fullname, name, prefix):

   if doc is None:
      func.__doc__ = fobj.__doc__
   else:
       
      d = fobj.__doc__
      L = d.split('\n')
      arg_lines = L[3:]

      L[0] = L[0].replace(fullname, name)
      L[0] = L[0].replace(prefix,'')
      if '.' in name:
          L[0] = L[0].replace(name[:name.index('.')].lower()+'_', '')
          L[1] = L[1].replace(name[:name.index('.')].lower()+'_', '')
      L[1] = L[1].replace(prefix, '')
      if '=' in L[1]:
          L[1] = L[1][:L[1].index('=')] + L[1][L[1].index('='):].replace(fullname, name)
      else:
          L[1] = L[1].replace(fullname.lower()+'_', '')
          L[1] = L[1].replace(fullname, name)

      if arg_lines:
          indent = max([len(s) for s in arg_lines])
          fmt = '  %%-%ds %%s' % indent

          for arg in doc['args']:
             name = arg['name'].lower()
             if name in badnames: name = badnames[name]

             for i, line in enumerate(arg_lines):
                 if line.startswith('  %s :' % name): break
             else:
                 raise ValueError('%s not found in lines %r' % (name, arg_lines))

             arg_lines[i] = arg_lines[i].replace(name,name[len(prefix):])

             if arg['type'].startswith('type('):
                arg_lines[i] = '  %s : %s object' % (name[len(prefix):], arg['type'][5:-1])

             if arg['doc'] != '':
                arg_lines[i] = (fmt % (arg_lines[i].strip(),arg['doc'])).replace('\n','\n   '+indent*' ')

      func.__doc__ = '%s\n\n%s' % ('\n'.join(L[:3]+arg_lines), doc['doc'])
      func.__name__ = fullname
      func._fobj = fobj

   return func

def wraproutine(modobj, moddoc, name, shortname, prefix):
   doc = moddoc['routines'][name]
   fobj = getattr(modobj, prefix+name)
   inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'])
   outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])
                               
   def func(*args, **kwargs):
       newargs, newkwargs = process_in_args(args, kwargs, inargs, prefix)

       try:
           res = fobj(*newargs, **newkwargs)
       except RuntimeError:
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

       newres = process_results(res, args, kwargs, inargs, outargs, prefix)
      
       return newres

   return add_doc(func, fobj, doc, name, shortname, prefix)



