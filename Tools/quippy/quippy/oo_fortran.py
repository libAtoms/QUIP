import logging, numpy, sys
from arraydata import arraydata
from farray import FortranArray
import numpy

major, minor = sys.version_info[0:2]

if (major, minor) < (2, 5):
    all = lambda seq: not False in seq
    any = lambda seq: True in seq


#logging.root.setLevel(logging.DEBUG)

py_keywords = ['and',       'del',       'from',      'not',       'while',    
               'as',        'elif',      'global',    'or',        'with',     
               'assert',    'else',      'if',        'pass',      'yield',    
               'break',     'except',    'import',    'print',              
               'class',     'exec',      'in',        'raise',              
               'continue',  'finally',   'is',        'return',             
               'def',       'for',       'lambda',    'try']

scalar_types = ['real(dp)', 'integer', 'logical', 'complex', 'character']
numpy_scalar_types = ['d', 'i', 'i', 'complex', 'S']

numpy_to_fortran = dict(zip(numpy_scalar_types, scalar_types))
numpy_to_fortran['f'] = 'real(dp)' # will be converted up by f2py

class FortranDerivedType(object):
   """Abstract base class for all fortran derived types."""
   pass

FortranDerivedTypes = {}

from numpy.f2py.crackfortran import badnames

def fix_badnames(x):
    return badnames.get(x,x)

def is_scalar_type(t):
   if t in scalar_types: return True
   if t.startswith('character'): return True
   if t.startswith('integer'): return True
   return False


def wrap_all(topmod, spec, mods, short_names, default_init_args, base_classes={}):
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
                                           default_init_args=default_init_args,
                                           params=params, base_classes=base_classes)
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
               
            if routine[:len(cls.__name__)+1] == cls.__name__.lower()+'_':
               method_name = routine[len(cls.__name__)+1:]
            else:
               method_name = routine

            if method_name in py_keywords: name = name+'_'
         
            setattr(cls, method_name, wraproutine(modobj, modspec, routine, cls.__name__+'.'+method_name))
            logging.debug('  added method %s to class %s' % (method_name, cls.__name__))
                  
         else:
            standalone_routines.append((routine,wraproutine(modobj, modspec, routine, routine)))
            logging.debug('  added stand-alone routine %s' % routine)

   return all_classes, standalone_routines, params


def wrapmod(modobj, moddoc, short_names, default_init_args, params, base_classes):

   routines = [x for x in moddoc['routines'].keys()]
   types =  [x for x in moddoc['types'].keys()]

   classes = []
   for cls in types:
      lcls = cls.lower()
      if '%s_initialise' % lcls in routines: 
         constructor = '%s_initialise' % lcls
      elif '%s_allocate' % lcls in routines:
         constructor = '%s_allocate' % lcls
      elif lcls in short_names:
         if '%s_initialise' % short_names[lcls] in routines:
            constructor = '%s_initialise'  % short_names[lcls]
         elif '%s_allocate' % short_names[lcls] in routines:
            constructor = '%s_allocate'  % short_names[lcls]
      else:
         cands = filter(lambda x: x.startswith('%s_initialise' % lcls), routines)
         if cands:
            constructor = cands[0] # pick the first one, for now
         else:
            logging.debug("Can't find constructor for type %s" % cls)
            continue

      if '%s_finalise' % lcls in routines:
         destructor = '%s_finalise' % lcls
      elif lcls in short_names and '%s_finalise' % short_names[lcls] in routines:
         destructor = '%s_finalise' % short_names[lcls]
      else:
         logging.debug("Can't find destructor for type %s" % cls)
         continue

      logging.debug('  Class %s' % cls)

      methods = [ x for x in moddoc['routines'].keys()
                  if (len(moddoc['routines'][x]['args']) > 0 and 
                      moddoc['routines'][x]['args'][0]['type'].lower() == 'type(%s)' % lcls) ]
      
      d = {}
      d['_classdoc'] = moddoc['types'][cls]
      d['__init__'] =  wrapinit(cls, modobj, moddoc, constructor, default_init_args)
      d['__del__'] = wrapdel(cls, modobj, moddoc, destructor)
      d['__str__'] = wraprepr(cls)
      d['_update'] = _update
      d['_arrays'] = {}
      d['_subobjs'] = {}
      interfaces = {}

      for name in methods:
         fullname = name
         if name[:len(lcls)+1] == lcls+'_':
            name = name[len(lcls)+1:]

         if lcls in short_names:
             if name[:len(short_names[lcls])+1] == short_names[lcls]+'_':
                name = name[len(short_names[lcls])+1:]

         if name in py_keywords: name = name+'_'
         
         d[name] = wraproutine(modobj, moddoc, fullname, cls+'.'+name)

         for intf_name,value in moddoc['interfaces'].iteritems():
            intf_routines = value['routines']
            if fullname in intf_routines:
               if not intf_name in interfaces: interfaces[intf_name] = []
               interfaces[intf_name].append((cls+'.'+name, moddoc['routines'][fullname], d[name]))

         logging.debug('    adding method %s' % name)
         del routines[routines.index(fullname)]

      d['__doc__'] = moddoc['doc']

      # only keep interfaces with more than one routine in them
      # if there's just one routine we want to copy the docstring
      d['_interfaces'] = {}
      for name,value in interfaces.iteritems():
         if len(value) > 1: 
             d['_interfaces'][name] = value
         else:
             value[0][2].__doc__  ='\n'.join(moddoc['interfaces'][name]['doc']) + '\n\n' + value[0][2].__doc__

      for intf, value in d['_interfaces'].iteritems():
         # if there's a routine with same name as interface, move it first
         if intf in d:
            d[intf+'_'] = d[intf]
         d[intf] = wrap_interface(cls, intf, value, moddoc['interfaces'][intf]['doc'])

      # Constructor documentation could be in initialise interface
      if 'initialise' in moddoc['interfaces']:
          if moddoc['interfaces']['initialise']['doc'] != []:
              d['__init__'].__doc__ = '\n'.join(moddoc['interfaces']['initialise']['doc'])+ '\n\n' + d['__init__'].__doc__

      # Add property() calls for get and set routines of scalar type,
      # and placeholders for dynamically created sub objects
      
      for name, el in moddoc['types'][cls]['elements'].iteritems():
         if name == 'thetype': name = 'type' # special case: f2py misparses name 'type'
         name = name.lower() # lower case all attributes

         if 'get' in el and 'set' in el:
            if is_scalar_type(el['type']):
               logging.debug('    adding property %s get=%s set=%s' % (name, el['get'], el['set']))

               d['_get_%s' % name] = wrap_get(modobj, el['get'])
               d['_set_%s' % name] = wrap_set(modobj, el['set'])

               d[name] = property(fget=d['_get_%s'%name],fset=d['_set_%s'%name],doc=el['doc'])
            else:
               d['_subobjs'][name] = (el['type'], 
                                     getattr(modobj, el['get']), 
                                     getattr(modobj, el['set']))
               d[name] = FortranDerivedType()
               d[name].__doc__ = el['doc']

         # Placeholders for array members
         if 'array' in el:
            logging.debug('    adding array %s' % name)
            d['_arrays'][name] = (getattr(modobj, el['array']), el['doc'])
            d[name] = FortranArray(doc=el['doc'])
   
         if cls in base_classes:
            new_cls = type(object)(cls[0].upper()+cls[1:], (FortranDerivedType,base_classes[cls]), d)
         else:
             new_cls = type(object)(cls[0].upper()+cls[1:], (FortranDerivedType,), d)
         FortranDerivedTypes['type(%s)' % cls.lower()] = new_cls
         classes.append((cls[0].upper()+cls[1:], new_cls))

   # Try to evaluate params
   evaldict = numpy.__dict__
   evaldict['huge'] = lambda x: numpy.finfo('d').max
   evaldict['cmplx'] = lambda x,y,kind: complex(x,y)
   evaldict.update(params)
   for (name, ftype, attributes, value) in moddoc['parameters']:
       code = value.replace('(/','array([')
       code = code.replace('/)','])')
       code = code.replace('_dp', '')
       params[name] = eval(code, evaldict)
       evaldict[name] = params[name]
       logging.debug('  adding parameter %s' % name)

   
   return (classes, routines, params)


def wrap_get(modobj, name):
   fobj = getattr(modobj, name)
                               
   def func(self):
      res = call_fortran(fobj, (self._p,))
      return res

   return func

def wrap_set(modobj, name):
   fobj = getattr(modobj, name)
                               
   def func(self, value):
      res = call_fortran(fobj, (self._p,value))
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

   return func

def process_in_args(args, kwargs, inargs):
   newargs = []
   for arg, spec in zip(args,inargs):
      if spec['type'].startswith('type'):
         if isinstance(arg, FortranDerivedTypes[spec['type'].lower()]):
            newargs.append(arg._p)
         else:
            raise TypeError('Argument %s should be of type %s' % (arg, spec['type']))
      else:
         newargs.append(arg)

   type_lookup = dict([(fix_badnames(x['name'].lower()) ,x['type']) for x in inargs])
   newkwargs = {}
   for k,a in kwargs.iteritems():
      if not k in type_lookup: raise ValueError('Unknown keyword argument %s' % k)
      if type_lookup[k].startswith('type'):
         if isinstance(a, FortranDerivedTypes[type_lookup[k].lower()]):
            newkwargs[k] = a._p
         else:
            raise TypeError('Argument %s=%s should be of type %s' % (k,a,type_lookup[k]))
      else:
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
             newres.append(FortranDerivedTypes[spec['type'].lower()](p=r,finalise=True))
          else:
             newres.append(r)

   # update any objects in args or kwargs affected by this call
   for arg, spec in zip(args, inargs):
      if (isinstance(arg, FortranDerivedType) and 'pointer' in spec['attributes'] and
          not 'fintent(in)' in spec['attributes']):
          arg._update()

   type_lookup = dict([(fix_badnames(x['name'].lower()),x) for x in inargs])
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

def wraproutine(modobj, moddoc, name, shortname):
   doc = moddoc['routines'][name]
   fobj = getattr(modobj, name)
   inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'])
   outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])
                               
   def func(*args, **kwargs):
      newargs, newkwargs = process_in_args(args, kwargs, inargs)

      res = call_fortran(fobj, newargs, newkwargs)
      newres = process_results(res, args, kwargs, inargs, outargs)
      
      return newres

   return add_doc(func, fobj, doc, name, shortname)

def wraprepr(cls):
   def func(self):
      items = []
      for k in dir(self):
         if k.startswith('_'): continue
         v = getattr(self,k)
         if type(v).__name__ == 'instancemethod': continue
         if isinstance(v, FortranDerivedType): 
            items.append('%s=%s()' % (k,v.__class__.__name__))
         else:
            items.append('%s=%s' % (k,v))
      return '%s(%s)' % (cls, ',\n'.join(items))

   return func

def wrapinit(cls, modobj, moddoc, name, default_init_args):
   doc = moddoc['routines'][name]
   fobj = getattr(modobj, name)
   inargs  = filter(lambda x: not 'intent(out)' in x['attributes'], doc['args'][1:])
   inargs.append({'name':'p','type':'integer','attributes':[]})
   inargs.append({'name':'finalise','type':'logical','attributes':[]})
   outargs = filter(lambda x: 'intent(out)' in x['attributes'], doc['args'])
                               
   def func(self, *args, **kwargs):
      logging.debug('init %s %s %s ' % (cls, args, kwargs))
      newargs, newkwargs = process_in_args(args, kwargs, inargs)

      if newargs == () and cls in default_init_args:
         newargs = default_init_args[cls]

      if not 'p' in newkwargs:
         newkwargs2 = newkwargs.copy()
         if 'finalise' in newkwargs2: del newkwargs2['finalise']
         self._p = call_fortran(fobj, newargs, newkwargs2)
      else:
         self._p = newkwargs['p']

      self._finalise = True
      if 'finalise' in newkwargs:
         self._finalise = newkwargs['finalise']

      self._update()

   return add_doc(func, fobj, doc, name, cls)


def wrapdel(cls, modobj, moddoc, name):
   doc = moddoc['routines'][name]
   fobj = getattr(modobj, name)
                               
   def func(self):
      logging.debug('del %s' % cls)

      if self._finalise:
         call_fortran(fobj, (self._p,))
      self._p = None
 
   return add_doc(func, fobj, doc, name, 'del')


def fortran_equivalent_type(t):
   tmap = {'int': 'integer', 
           'float': 'real(dp)',
           'bool': 'logical',
           'str': 'character'}

   if type(t).__name__ in tmap:
      return (tmap[type(t).__name__], 'scalar')
   elif hasattr(t,'__iter__'):
      a = numpy.array(t)
      if a.dtype.kind in numpy_to_fortran:
         return (numpy_to_fortran[a.dtype.kind], a.shape)
      else:
         raise TypeError('Unknown array type %s' % a.dtype.kind)      
   elif isinstance(t, FortranDerivedType):
      return ('type(%s)' % t.__class__.__name__.lower(), 'scalar')
   else:
      raise TypeError('Unknown type %s' % type(t))
          
def type_is_compatible(spec, arg):
    arg_type, dims = fortran_equivalent_type(arg)

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

def wrap_interface(cls, name, routines, doc):

   def func(*args, **kwargs):
      logging.debug('Interface %s %s' % (cls, name))

      for rname, spec, routine in routines:
         
         inargs = filter(lambda x: not 'intent(out)' in x['attributes'], spec['args'])
         oblig_inargs = filter(lambda x: not 'optional' in x['attributes'], inargs)
         opt_inargs = filter(lambda x: 'optional' in x['attributes'], inargs)

         # Check number of arguments is compatible
         # Should be in range len(oblig_inargs) <= L <= len(oblig_inargs) + len(opt_inargs)
         if (len(args)+len(kwargs) < len(oblig_inargs) or
             len(args)+len(kwargs) > len(oblig_inargs) + len(opt_inargs)):
            continue

         newinargs = (oblig_inargs + opt_inargs)[:len(args)]

         # Check types and dimensions are compatible
         if not all([type_is_compatible(spec, a) for (spec,a) in zip(newinargs, args)]):
            continue

         # Check keyword arguments, if incompatible continue to next
         # candidate interface
         if kwargs:
            innames = [fix_badnames(x['name']) for x in oblig_inargs + opt_inargs]

            if not all([key in innames[len(args):] for key in kwargs.keys()]):
               continue

            try:
               for key, arg in kwargs.iteritems():
                  (inarg,) = [x for x in inargs if fix_badnames(x['name']) == key]
                  if not type_is_compatible(inarg, arg):
                     raise ValueError

            except ValueError:
               continue
         
         return routine(*args, **kwargs)
            
            
      raise TypeError('No matching routine found in interface %s' % name)
         
   func.__doc__ = '\n'.join(doc)+'\n\n'

   for rname,spec,routine in routines:
       func.__doc__ += routine.__doc__.replace(rname,cls+'.'+name)

   return func


def _update(self):
   logging.debug('updating %s at 0x%x' % (self.__class__, id(self)))
   for name in self._arrays:
      try:
         delattr(self, name)
      except AttributeError:
         pass

   for name, (arrayfunc, doc) in self._arrays.iteritems():
      dtype, shape, loc = arrayfunc(self._p)
      dtype = dtype.strip()
      if shape.any() and loc != 0:
         if dtype in numpy_scalar_types:
             if hasattr(self, '_map_array_shape'):
                 nshape = self._map_array_shape(name, shape)
                 setattr(self, '_'+name, FortranArray(arraydata(shape,dtype,loc), doc))
                 setattr(self, name, getattr(self, '_'+name)[nshape])
             else:
                 setattr(self, name, FortranArray(arraydata(shape,dtype,loc), doc))
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
      p = getfunc(self._p)
      if p != 0:
         savedoc = getattr(self, name).__doc__
         setattr(self, name, FortranDerivedTypes[cls.lower()](p=p,finalise=False))
         getattr(self,name).__doc__ = savedoc

   if hasattr(self,'update_hook'): 
       self.update_hook()
      
         

def call_fortran(fobj, args=(), kwargs={}):
   return fobj(*args, **kwargs)
   
