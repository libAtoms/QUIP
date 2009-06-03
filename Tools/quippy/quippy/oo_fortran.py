"""Object oriented interface on top of f2py generated wrappers."""

import logging, numpy, sys
from arraydata import arraydata
from farray import *
import numpy
from types import MethodType

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
    'complex': 'complex',
    'S': 'character',
    'f': 'real(dp)'
    }

FortranDerivedTypes = {}

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
           }

   if type(t).__name__ in tmap:
      return (tmap[type(t).__name__], 'scalar')
   elif isinstance(t, FortranDerivedType):
      return ('type(%s)' % t.__class__.__name__.lower(), 'scalar')
   elif hasattr(t,'__iter__') or hasattr(t,'dtype'):
      a = numpy.array(t)
      if a.dtype.kind in numpy_to_fortran:
         return (numpy_to_fortran[a.dtype.kind], a.shape)
      else:
         raise TypeError('Unknown array type %s' % a.dtype.kind)      
   else:
      raise TypeError('Unknown type %s' % type(t))
          
def type_is_compatible(spec, arg):
    arg_type, dims = fortran_equivalent_type(arg)
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
            return len(adims) == len(dims)
        else:
            # Check each dimension matches
            return all([x == y for (x,y) in zip(adims,dims)])

def process_in_args(args, kwargs, inargs):
   newargs = []
   for arg, spec in zip(args,inargs):
      if spec['type'].startswith('type'):
         if isinstance(arg, FortranDerivedTypes[spec['type'].lower()]):
            newargs.append(arg._fpointer)
         else:
            raise TypeError('Argument %s should be of type %s' % (arg, spec['type']))
      else:
         newargs.append(arg)

   type_lookup = dict([(badnames.get(x['name'].lower(),x['name'].lower()),x['type']) for x in inargs])
   newkwargs = {}
   for k,a in kwargs.iteritems():
      if not k in type_lookup: raise ValueError('Unknown keyword argument %s' % k)
      if type_lookup[k].startswith('type'):
         if isinstance(a, FortranDerivedTypes[type_lookup[k].lower()]):
            newkwargs[k] = a._fpointer
         elif a is None:
             continue
         else:
            raise TypeError('Argument %s=%s should be of type %s' % (k,a,type_lookup[k]))
      else:
          if a is None: continue
          newkwargs[k] = a

   return tuple(newargs), newkwargs

def process_results(res, args, kwargs, inargs, outargs):
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
          else:
             newres.append(r)

   # update any objects in args or kwargs affected by this call
   for arg, spec in zip(args, inargs):
      if (isinstance(arg, FortranDerivedType) and 'pointer' in spec['attributes'] and
          not 'fintent(in)' in spec['attributes']):
          arg._update()

   type_lookup = dict([(badnames.get(x['name'].lower(),x['name'].lower()),x) for x in inargs])
   for k,a in kwargs.iteritems():
      if (isinstance(a, FortranDerivedType) and 'pointer' in type_lookup[k]['attributes']
          and not 'fintent(in)' in type_lookup[k]['attributes']):
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
   _elements = {}

   def __init__(self, *args, **kwargs):

       call_init = True
       self._fpointer = None
       self._finalise = True
       
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
       #print 'del %s(fpointer=0x%x, finalise=%d)' % (self.__class__.__name__, self._fpointer, self._finalise)

       if self._fpointer and self._finalise and '__del__' in self._routines:
           fobj, doc = self._routines['__del__']
           fobj(self._fpointer)
           self._fpointer = None
 

   def __repr__(self):
       return '%s(fpointer=0x%x, finalise=%d)' % (self.__class__.__name__, self._fpointer, self._finalise)

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


   def _update(self):
       logging.debug('updating %s at 0x%x' % (self.__class__, id(self)))
       if self._fpointer is None: return
       
       for name in self._arrays:
          try:
             delattr(self, name)
          except AttributeError:
             pass

       for name, (arrayfunc, doc) in self._arrays.iteritems():
          dtype, shape, loc = arrayfunc(self._fpointer)
          dtype = dtype.strip()
          if shape.any() and loc != 0:
             if dtype in numpy_to_fortran.keys():
                 nshape = self._get_array_shape(name)
                 if nshape is not None:
                     setattr(self, '_'+name, FortranArray(arraydata(shape,dtype,loc), doc))
                     #setattr(self, '_'+name, arraydata(shape,dtype,loc))
                     setattr(self, name, getattr(self, '_'+name)[nshape])
                 else:
                     setattr(self, name, FortranArray(arraydata(shape,dtype,loc), doc))
                     #setattr(self, name, arraydata(shape,dtype,loc))
             else:
                # access to arrays of derived type not yet implemented
                continue
          else:
             logging.debug('Ignoring uninitialised array %s (shape=%s, dtype=%s, loc=%s)' % (name, shape, dtype, loc))
             continue

       for name in self._subobjs:
          try:
             delattr(self, name)
          except AttributeError:
             pass

       for name, (cls, getfunc, setfunc) in self._subobjs.iteritems():
          if name == 'thetype': name = 'type' # special case: f2py misparses name 'type'
          if not cls.lower() in FortranDerivedTypes:
             logging.debug('Unknown class %s' % cls)
             continue
          p = getfunc(self._fpointer)
          if p != 0:
             savedoc = getattr(self, name).__doc__
             setattr(self, name, FortranDerivedTypes[cls.lower()](fpointer=p,finalise=False))
             getattr(self,name).__doc__ = savedoc

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
       newargs, newkwargs = process_in_args(args, kwargs, inargs)
       
       res = fobj(*newargs, **newkwargs)

       if not name.startswith('__init__'):
           newres = process_results(res, args, kwargs, inargs, outargs)
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
               innames = [badnames.get(x['name'],x['name']) for x in oblig_inargs + opt_inargs]

           if not all([key in innames[len(args):] for key in kwargs.keys()]):
               logging.debug('Unexpected keyword argument valid=%s got=%s' % (kwargs.keys, innames[len(args):]))
               continue

           try:
               for key, arg in kwargs.iteritems():
                   (inarg,) = [x for x in inargs if badnames.get(x['name'],x['name']) == key]
                   if not type_is_compatible(inarg, arg):
                       logging.debug('Types and dimensions of kwarg %s incompatible' % key)
                       raise ValueError
                   
           except ValueError:
               continue

           logging.debug('calling '+rname)
           return routine(self, *args, **kwargs)


       raise TypeError('No matching routine found in interface %s' % name)


def wrap_all(topmod, spec, mods, short_names):
   all_classes = []
   leftover_routines = []
   standalone_routines = []
   params = {}

   for mod in mods:
       try:
           fmod = getattr(topmod, mod)
       except AttributeError:
           fmod = None
       logging.debug('Module %s' % mod)
       fmod.__doc__ = spec[mod]['doc']
       classes, routines, params = wrapmod(fmod, spec[mod],
                                           short_names=short_names,
                                           params=params)
       all_classes.extend(classes)
       leftover_routines.append((fmod, spec[mod], routines))

   # Now try to add orphaned routines to classes where possible
   for modobj, modspec, routines in leftover_routines:
      for routine in routines:
         rspec = modspec['routines'][routine]

         if (len(rspec['args']) > 0 and 
             rspec['args'][0]['type'].lower() in FortranDerivedTypes and
             rspec['args'][0]['name'].lower() == 'this'):
            
            cls = FortranDerivedTypes[rspec['args'][0]['type'].lower()]

            lclsname = cls.__name__.lower()[len(fortran_class_prefix):]

            if routine[:len(lclsname)+1] == lclsname+'_':
               method_name = routine[len(lclsname)+1:]
            else:
               method_name = routine

            if method_name in py_keywords: name = name+'_'
         
            setattr(cls, method_name, wraproutine(modobj, modspec, routine, cls.__name__+'.'+method_name))
            logging.debug('  added method %s to class %s' % (method_name, cls.__name__))
                  
         else:
            standalone_routines.append((routine,wraproutine(modobj, modspec, routine, routine)))
            logging.debug('  added stand-alone routine %s' % routine)

   return all_classes, standalone_routines, params



def wrapmod(modobj, moddoc, short_names, params):

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


      constructors = [ x for x in methods if (x.startswith('%s_initialise' % lcls) or
                                              x.startswith('%s_allocate' % lcls))]

      destructors = [ x for x in methods if x.startswith('%s_finalise' % lcls)]


      if lcls in short_names:
          constructors += [ x for x in methods if(x.startswith('%s_initialise' % short_names[lcls]) or
                                                  x.startswith('%s_allocate' % short_names[lcls]))
                            and 'intent(out)' in moddoc['routines'][x]['args'][0]['attributes'] ]

          destructors += [ x for x in methods if x.startswith('%s_finalise' % short_names[lcls]) ]

      if (len(constructors) == 0):
          logging.warning("Can't find constructor for type %s. Skipping class" % cls)
          continue
          
      if (len(destructors) == 0):
          logging.warning("Can't find destructor for type %s. Skipping class" % cls)
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


      new_cls._classdoc = moddoc['types'][cls]
      new_cls._moddoc = moddoc
      new_cls._modobj = modobj

      if constructor:
          new_cls._routines['__init__'] = (getattr(modobj, constructor), moddoc['routines'][constructor])
          new_cls.__init__ = add_doc(wrapinit(constructor), getattr(modobj, constructor),
                                     moddoc['routines'][constructor], constructor, tcls)

      if destructor:
          new_cls._routines['__del__'] = (getattr(modobj, destructor), moddoc['routines'][destructor])

      interfaces = {}

      for name in methods:
         fullname = name
         if name[:len(lcls)+1] == lcls+'_':
            name = name[len(lcls)+1:]

         if lcls in short_names:
             if name[:len(short_names[lcls])+1] == short_names[lcls]+'_':
                name = name[len(short_names[lcls])+1:]

         if name in py_keywords: name = name+'_'

         fobj = getattr(modobj, fullname)
         doc = moddoc['routines'][fullname]
         if fullname in constructors:
             name = '__init__'+name
             logging.debug('    adding constructor %s.%s' % (tcls, name))
                          

         new_cls._routines[name] = (fobj, doc)
         func = add_doc(wrapmethod(name), fobj, doc, name, tcls+'.'+name)
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
      for name,value in interfaces.iteritems():
         if name.lower() == 'initialise': name = '__init__'
         if len(value) > 1: 
             new_cls._interfaces[name] = value
         else:
             pass
#             value[0][2].__doc__  ='\n'.join(moddoc['interfaces'][name]['doc']) + '\n\n' + value[0][2].__doc__

      for intf, value in new_cls._interfaces.iteritems():
         # if there's a routine with same name as interface, move it first
         if hasattr(new_cls, intf):
            setattr(new_cls,intf+'_', getattr(new_cls, intf))

         func = wrapinterface(intf)
         docname = intf
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


      # Constructor documentation could be in initialise interface
      #if 'initialise' in moddoc['interfaces']:
      #    if moddoc['interfaces']['initialise']['doc'] != []:
      #        d['__init__'].__doc__ = '\n'.join(moddoc['interfaces']['initialise']['doc'])+ '\n\n' + d['__init__'].__doc__

      # Add property() calls for get and set routines of scalar type,
      # and placeholders for dynamically created sub objects
      
      for name, el in moddoc['types'][cls]['elements'].iteritems():
         if name == 'thetype': name = 'type' # special case: f2py misparses name 'type'
         name = name.lower() # lower case all attributes

         if 'get' in el and 'set' in el:
            if is_scalar_type(el['type']):
               logging.debug('    adding property %s get=%s set=%s' % (name, el['get'], el['set']))

               new_cls._elements[name] = (getattr(modobj, el['get']), getattr(modobj, el['set']))
                  
               setattr(new_cls, '_get_%s' % name, wrap_get(name))
               setattr(new_cls, '_set_%s' % name, wrap_set(name))

               setattr(new_cls, name, property(fget=getattr(new_cls,'_get_%s'%name),
                                               fset=getattr(new_cls,'_set_%s'%name),
                                               doc=el['doc']))
            else:
                new_cls._subobjs[name] = (el['type'], 
                                          getattr(modobj, el['get']), 
                                          getattr(modobj, el['set']))
                setattr(new_cls, name,FortranDerivedType())
                getattr(new_cls, name).__doc__ = el['doc']

         # Placeholders for array members
         if 'array' in el:
            logging.debug('    adding array %s' % name)
            new_cls._arrays[name] = (getattr(modobj, el['array']), el['doc'])
            setattr(new_cls, name, FortranArray(doc=el['doc']))
   

   # Try to evaluate params
   evaldict = numpy.__dict__
   evaldict['huge'] = lambda x: numpy.finfo('d').max
   evaldict['cmplx'] = lambda x,y,kind: complex(x,y)
   evaldict['farray'] = farray
   for (name, ftype, attributes, value) in moddoc['parameters']:
       code = value.replace('(/','farray([')
       code = code.replace('/)','])')
       code = code.replace('_dp', '')
       params[name] = eval(code, evaldict)
       evaldict[name] = params[name]
       logging.debug('  adding parameter %s' % name)

   
   return (classes, routines, params)


def wrap_get(name):
                               
   def func(self):
       if self._fpointer is None:
           raise ValueError('%s object not initialised.' % self.__class__.__name__)
       if not name in self._elements:
           raise ValueError('Unknown element %s in class %s' % (name, self.__class.__name)) 
       res = self._elements[name][0](self._fpointer)
       return res

   return func

def wrap_set(name):
                              
   def func(self, value):
       if self._fpointer is None:
           raise ValueError('%s object not initialised.' % self.__class__.__name__)
       if not name in self._elements:
           raise ValueError('Unknown element %s in class %s' % (name, self.__class.__name)) 
       res = self._elements[name][1](self._fpointer, value)
       self._update()
       return res

   return func


def add_doc(func, fobj, doc, fullname, name):

   if doc is None:
      func.__doc__ = fobj.__doc__
   else:
      d = fobj.__doc__
      L = d.split('\n')
      arg_lines = L[3:]

      L[0] = L[0].replace(fullname, name)
      if '=' in L[1]:
          L[1] = L[1][:L[1].index('=')] + L[1][L[1].index('='):].replace(fullname, name)
      else:
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
             raise ValueError('%s not found in lines' % name)

         if arg['type'].startswith('type('):
            arg_lines[i] = '  %s : %s object' % (name, arg['type'][5:-1])

         if arg['doc'] != '':
            arg_lines[i] = (fmt % (arg_lines[i].strip(),arg['doc'])).replace('\n','\n   '+indent*' ')

      func.__doc__ = '%s\n\n%s' % ('\n'.join(L[:3]+arg_lines), doc['doc'])
      func.__name__ = fullname
      func._fobj = fobj

   return func

def wraproutine(modobj, moddoc, name, shortname):
   doc = moddoc['routines'][name]
   fobj = getattr(modobj, name)
   inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'])
   outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])
                               
   def func(*args, **kwargs):
      newargs, newkwargs = process_in_args(args, kwargs, inargs)

      res = fobj(*newargs, **newkwargs)
      newres = process_results(res, args, kwargs, inargs, outargs)
      
      return newres

   return add_doc(func, fobj, doc, name, shortname)



