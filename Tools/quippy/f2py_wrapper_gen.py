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

from f90doc import *
import numpy, re

def wrap_mod(mod, type_map, out=None, kindlines=[], initlines={}, filtertypes=None, prefix='',
             callback_routines={},public_symbols=None,sizeof_fortran_t=12):
   """Given an f90doc.C_module class 'mod', write a F90 wrapper file suitable for f2py to 'out'."""
   spec = {}
   def strip_type(t):
      if t.startswith('type('):
         t = t[t.index('(')+1:t.index(')')]
      return t.lower()

   def println(*args):
      out.write('%s%s\n' % ((' '*indent),' '.join(args)))

   def default_value(type):
      lookup = {'real(dp)': '-1.0_dp',
                'integer': '-1',
                'character*(1)': '"@"',
                'character*(*)': '"@"'}
      if type in lookup: return lookup[type]
      elif type.startswith('character'): return '"@"'
      elif type.startswith('type'): return '-1'

   # Skip things with callbacks (for now...)
   def no_callbacks(sub):
      try:
          types = [x.type for x in sub.arguments]
          attrs = [x.attributes for x in sub.arguments]
      except AttributeError:
         return sub.name.lower() in callback_routines.keys()

      return True

   # Also skip allocatable and pointer array arguments
   # (and return type if sub is a function)
   def no_allocatables_or_pointers(sub):
      if sub.name.endswith('initialise_ptr') or sub.name.endswith('finalise_ptr'):
         return True
      for arg in sub.arguments:
          if not hasattr(arg, 'attributes'): continue
         
          # FIXME: this skips scalar pointer args too

          if ('allocatable' in arg.attributes or 'pointer' in arg.attributes):
             if not 'optional' in arg.attributes:
                return False

          # arrays of derived type are out as well
          dims = filter(lambda x: x.startswith('dimension'), arg.attributes)
          if len(dims) > 0 and arg.type.startswith('type'):
              return False

      if isinstance(sub, C_funct):
          if 'allocatable' in sub.ret_val.attributes or 'pointer' in sub.ret_val.attributes:
              return False

      return True

   def no_complex_scalars(sub):
      for arg in sub.arguments:
          if not hasattr(arg, 'attributes'): continue         
          dims = filter(lambda x: x.startswith('dimension'), arg.attributes)
          if arg.type.startswith('complex') and len(dims) == 0: return False

      return True


   def no_c_ptr(sub):
      for arg in sub.arguments:
          if not hasattr(arg, 'attributes'): continue      
          if arg.type.lower() == 'type(c_ptr)': return False

      return True

   # Only wrap certain types
   def no_bad_types(sub):
      for arg in sub.arguments:
         if not hasattr(arg, 'attributes'): continue      
         if arg.type.startswith('type') and not strip_type(arg.type) in filtertypes:
            print 'omitting routine %s as argument %s of unwrapped type %s' % (sub.name, arg.name, arg.type)
            return False
      return True

   def no_private_symbols(sub):
      subname = sub_to_interface.get(sub.name, sub.name)
      return subname.lower() in [s.lower() for s in public_symbols]
      

   debug(mod.name)
   shortname = mod.name[:mod.name.index('_module')]

   if out is None:
      out = sys.stdout #open('%s.f90' %shortname, 'w')

   # Interfaces
   sub_to_interface = {}
   subts = filter(no_callbacks, mod.subts)
   functs = filter(no_callbacks, mod.functs)
   for intf in mod.interfaces:
      if intf.name.startswith('operator('): continue
      thissubs = filter(no_callbacks, intf.subts)
      thisfuncts = filter(no_callbacks, intf.functs)
      subts.extend(thissubs)
      functs.extend(thisfuncts)
      for sub in thissubs + thisfuncts:
         sub_to_interface[sub.name] = intf.name

   subts = filter(no_allocatables_or_pointers, subts)
   functs = filter(no_allocatables_or_pointers, functs)

   subts = filter(no_complex_scalars, subts)
   functs = filter(no_complex_scalars, functs)

   subts = filter(no_c_ptr, subts)
   functs = filter(no_c_ptr, functs)

   if filtertypes is not None:
      print 'Filtering routines to exclude derived types not in list %s' % filtertypes
      subts  = filter(no_bad_types, subts)
      functs = filter(no_bad_types, functs)

   if public_symbols is not None:
      print 'Filtering routines to exclude routines not in list of public symbols %s' % public_symbols
      subts = filter(no_private_symbols, subts)
      functs = filter(no_private_symbols, functs)

   debug('%s: %d subs' % (mod.name, len(subts+functs)))

   names = [x.name for x in subts+functs]
   if shortname in names: shortname += '_'

   indent = 0

   spec[shortname.lower()] = {'doc': '\n'.join(mod.doc),
                                  'routines': {},
                                  'types': {},
                                  'interfaces': {},
                                  'parameters': []}

   # for some reason setting parameters in Fortran causes an AssertionError
   # when importing module. For now, let's just copy them into python dictionary
   for el in mod.elements:
      if 'parameter' in el.attributes:
          spec[shortname.lower()]['parameters'].append((el.name, el.type, el.attributes, el.value))

   for sub in subts + functs:

      arglines = []
      init = []
      got_types = []
      uses = set()
      type_lines = []

      # Imported routine
      ##uses.append('use %s, only: imported_%s => %s' % (mod.name, sub.name,sub.name))
      uses.add(mod.name.lower())
      
      # Add argument for return value, after last non-optional argument
      if hasattr(sub,'ret_val'):
          ret_val = sub.ret_val
          ret_val.name = 'ret_'+ret_val.name
          ret_val.attributes.append('intent(out)')
          ret_val.is_ret_val = True
          sub.arguments = ([ arg for arg in sub.arguments if not hasattr(arg, 'attributes') or (hasattr(arg, 'attributes') and not 'optional' in arg.attributes) ] + [ret_val] + 
                           [ arg for arg in sub.arguments if hasattr(arg, 'attributes') and 'optional' in arg.attributes
                             and not 'pointer' in arg.attributes
                             and not 'allocatable' in arg.attributes ])
      else:
         sub.arguments = [arg for arg in sub.arguments if not hasattr(arg, 'attributes') or not 'optional' in arg.attributes or
                          (not 'pointer' in arg.attributes and not 'allocatable' in arg.attributes)]


      args = sub.arguments
      old_arg_names = [x.name for x in args]
      for a in args:
         a.name = prefix + a.name
      argnames = [x.name for x in args][:]
      newargnames = argnames[:]

      optionals = []

      n_dummy = 0

      newname = sub.name
      spec[shortname.lower()]['routines'][newname.lower()] = \
          {'doc': '\n'.join(sub.doc),'args':[]}
      thisdoc = spec[shortname.lower()]['routines'][newname.lower()]['args']

      # See if this routine is in any interfaces
      for intf in mod.interfaces:
          subnames = [x.name.lower() for x in intf.subts + intf.functs]
          if sub.name.lower() in subnames:
              if not intf.name.lower() in spec[shortname.lower()]['interfaces']:
                  spec[shortname.lower()]['interfaces'][intf.name.lower()] = \
                      {'doc': intf.doc,
                       'routines': []}

              spec[shortname.lower()]['interfaces'][intf.name.lower()]['routines'].append(sub.name.lower())

      allocates = []
      callbacklines = []
      transfer_in = []
      transfer_out = []

      for argindex, arg in enumerate(args):

          append_argline = True

          if not hasattr(arg, 'attributes') and not hasattr(arg, 'type'):
             callback_spec = callback_routines[sub.name.lower()]
             arglines.extend(callback_spec['arglines'])
             thisdoc.append({'doc': '\n'.join(arg.doc), 'name':arg.name, 'type': 'callback', 'attributes': callback_spec['attributes']})
             callbacklines.extend(['if (.false.) then','call %s(%s)' % (arg.name, callback_spec['call']),'end if'])
             continue

          attributes = arg.attributes
          mytype = arg.type

          if 'optional' in attributes and 'intent(out)' in attributes: 
              attributes.remove('intent(out)')
              attributes.append('intent(inout)')

          if arg.type.startswith('type'):

              tname = strip_type(arg.type)

              if not tname in got_types:
                 ##uses.append('use %s, only: imported_%s => %s' % (type_map[tname], tname, tname))
                 uses.add(type_map[tname].lower())

                 type_lines.extend(['type %s_ptr_type' %  tname,
                                    'type(%s), pointer :: p => NULL()' % tname,
                                    'end type %s_ptr_type' % tname])
                 got_types.append(tname)
                          
              if tname in initlines:
                 use, (exe, exe_optional) = initlines[tname]
                 if use is not None: uses.add(use.lower())
                 if 'optional' in attributes:
                    init.append(exe_optional % {'OLD_ARG':arg.name[len(prefix):], 'ARG':arg.name, 'PTR':arg.name+'_ptr%p'})
                 else:
                    init.append(exe % {'OLD_ARG':arg.name[len(prefix):], 'ARG':arg.name, 'PTR':arg.name+'_ptr%p'})

              ptr_type = 'type(%s_ptr_type)' % tname

              # Preserve original fortran intent
              fintent = [ a for a in attributes if a.startswith('intent')]
              if fintent != []: fintent = fintent[0].replace('intent', 'fintent')

              if ('intent(out)' in attributes or 
                  ((sub.name.lower().find('_initialise') != -1 or sub.name.lower().find('_allocate') != -1) \
                   and len(argnames) > 0 and argnames[0] == prefix+'this' and arg.name == prefix+'this')):
                 if not 'pointer' in attributes:
                    allocates.append(arg.name)
                 intent = 'intent(out)'
                 transfer_out.append(arg.name)
              else:
                 intent = 'intent(in)'
                 transfer_in.append((arg.name, 'optional' in attributes))

              arglines.append('integer, %s%s :: %s(%d)' % (intent,
                                                            ('optional' in attributes and ', optional' or ''),
                                                            arg.name, sizeof_fortran_t))
              
              arglines.append('%s :: %s_ptr' % (ptr_type, arg.name))
              argnames[argindex] = '%s_ptr%%p' % arg.name
              append_argline = False

              
          elif arg.type.startswith('character'):
              # change from '(len=*)' or '(*)' syntax to *(*) syntax
              try:
                  lind = arg.type.index('(')
                  rind = arg.type.rindex(')')
                  mytype = arg.type[:lind]+'*'+arg.type[lind:rind+1].replace('len=','')
                  #mytype = 'character*(*)'

                  # Try to get length of string arguments
                  if not mytype[11:-1] == '*' and not all([x in '0123456789' for x in mytype[11:-1]]):
                      try:
                          mytype = 'character*(%s)' % string_lengths[mytype[11:-1].lower()]
                      except KeyError:
                          mytype = 'character*(%s)' % string_lengths['default']

                  attributes = filter(lambda x: x.startswith('intent') or
                                      x.startswith('dimension') or
                                      x.startswith('optional') or x.startswith('pointer'), attributes)

                  # Default string length for intent(out) strings 
                  if mytype[11:-1] == '*' and 'intent(out)' in attributes:
                      mytype = 'character*(%s)' % string_lengths['default']


              except ValueError:
                  mytype = 'character*(1)'
      
          dims = filter(lambda x: x.startswith('dimension('), attributes)
          if dims != []:
              # replace dimensions with n1,n2
              dim = dims[0][10:-1]
              br = 0
              d = 1
              ds = ['']
              for c in dim:
                  if c != ',': ds[-1] += c
                  if   c == '(': br += 1
                  elif c == ')': br -= 1
                  elif c == ',':
                      if br == 0: ds.append('')
                      else: ds[-1] += ','

              newds = []
              for i,d in enumerate(ds):
                  if valid_dim_re.match(d):
                      #if ',' in d: ds[i] = d.replace('size','shape')
                      for old, new in zip(old_arg_names, argnames):
                         ds[i] = ds[i].replace(old, new)
                      if d.startswith('len'):
                          arglines.append('!f2py %s %s, dimension(%s) :: %s' % \
                                              (arg.type, 
                                               ','.join(filter(lambda a: not a.startswith('dimension'), attributes)), 
                                               ds[i].replace('len','slen'), arg.name))
                      continue
                  ds[i] = (prefix+'n%d' % (n_dummy))
                  newds.append(ds[i])
                  n_dummy += 1

              attributes[attributes.index(dims[0])] = 'dimension(%s)' % \
                                                          ','.join(ds)
              if 'allocatable' in attributes:
                  attributes.remove('allocatable')


              # intent(out) arrays of variable size don't work with f2py
              #if 'intent(out)' in attributes:
              #    attributes.remove('intent(out)')
              #    attributes.append('intent(inout)')

          charflag = None
          if mytype == 'character*(*)' and 'intent(out)' in attributes:
              mytype = 'character*(n%d)' % n_dummy
              charflag = 'n%d' % n_dummy
              n_dummy += 1

          if append_argline:
             if attributes == []:
                arglines.append('%s :: %s' % (mytype, arg.name))
             else:
                arglines.append('%s, %s :: %s' % (mytype, ', '.join(attributes),arg.name))

          f2py_attributes = attributes[:]
          # For types, we want the intent of the f2py 'pointer', rather than the real fortran intent
          # We also store the fortran intent as 'fintent', to help with determining which
          # objects are affected by a call.
          if arg.type.startswith('type'):
              f2py_attributes.append(intent)
              f2py_attributes.append(fintent)
          thisdoc.append({'doc': '\n'.join(arg.doc), 'name':arg.name, 'type': arg.type, 'attributes': f2py_attributes})

          if dims != []:
              for i,d in enumerate(newds):
                  newargnames.append(d)
                  arglines.insert(0,'integer :: %s' % d)
                  if not 'intent(out)' in attributes:
                      arglines.insert(0,'!f2py intent(hide), depend(%s) :: %s = shape(%s,%d)' % (arg.name, d, arg.name, i))
                  else:
                      thisdoc.append({'name': d, 'doc': 'shape(%s,%d)' % (arg.name,i), 'type': 'integer', 'attributes':[]})

          if charflag is not None:
              newargnames.append(charflag)
              arglines.append('integer :: %s' % charflag)
              if not 'intent(out)' in attributes:
                  arglines.append('!f2py intent(hide), depend(%s) :: %s = slen(%s)' % (arg.name, charflag, arg.name))
              else:
                  thisdoc.append({'name':charflag,'doc': 'slen(%s)' % arg.name, 'type': 'integer', 'attributes':[]})


          # Do not mess around with character*(1) multidimensional arrays
          if mytype == 'character*(1)' and dims != [] and 'intent(in)' in attributes:
             arglines.append('!f2py intent(cache) :: %s' % arg.name)


      # prepend each argument with its name so that we can leave out optionals in the middle of argument list
      for argindex, arg in enumerate(args):
          argnames[argindex] = '%s=%s' % (arg.name[len(prefix):], argnames[argindex])

      def join_and_split_lines(args, max_length=80):
          args_str = ''
          line = ''
          args = args[:] # make a copy
          while args:
              line += args.pop(0)
              if args:
                  line += ', '
              if args and len(line) >= max_length:
                  args_str += '%s&\n%s' % (line, ' '*(indent+3))
                  line = ''
          args_str += line
          return args_str

      println('subroutine %s%s(%s)' % (prefix, newname, join_and_split_lines(newargnames)))
      indent += 3
      for line in kindlines:
         if not any([line.startswith('use %s' % umod) for umod in uses]):
            println(line)
      for umod in uses:
         println('use %s' % umod)
      println('implicit none')
      for line in type_lines:
         println(line)
      for line in arglines:
          println(line)
      println()

      for var in allocates: println('allocate(%s_ptr%%p)' % var)

      for line in callbacklines: println(line)

      for var, optional in transfer_in:
         if optional:
            println('if (present(%s)) then' % var)
            indent += 3
            println('%s_ptr = transfer(%s, %s_ptr)' % (var,var,var))
            indent -= 3
            println('else')
            indent += 3
            println('%s_ptr%%p => null()' % var)
            indent -= 3
            println('end if')
         else:
            println('%s_ptr = transfer(%s, %s_ptr)' % (var,var,var))

      for line in init:
         println(line)

      subname = sub_to_interface.get(sub.name, sub.name)

      if subname == 'assignment(=)':
         assert len(argnames) == 2
         to_arg = argnames[0]
         to_arg = to_arg[to_arg.index('=')+1:]
         from_arg = argnames[1]
         from_arg = from_arg[from_arg.index('=')+1:]
         println('%s = %s' % (to_arg, from_arg))
      else:
         if hasattr(sub, 'ret_val'):
             argfilt = [ name for name,arg in zip(argnames,args) if not (hasattr(arg, 'is_ret_val') and arg.is_ret_val) ]
             if sub.ret_val.type.startswith('type'):
                println('%s_ptr%%p = %s(%s)' % (sub.ret_val.name, subname, join_and_split_lines(argfilt)))                          
             else:
                println('%s = %s(%s)' % (sub.ret_val.name, subname, join_and_split_lines(argfilt)))
         else:
             println('call %s(%s)' % (subname, join_and_split_lines(argnames)))

      for var in transfer_out: println('%s = transfer(%s_ptr, %s)' % (var,var,var))


      if sub.name.lower().endswith('finalise') and len(old_arg_names) > 0 and old_arg_names[0] == 'this':
         println('deallocate(%sthis_ptr%%p)' % prefix)

      println()

      indent -= 3
      println('end subroutine %s%s' % (prefix,newname))
      println()


   # add _get_<type> and _set_<type> methods

   numpy_type_map = {'real(dp)':'d',
                     'integer':'i',
                     'logical':'i',
                     'character*(*)':'S',
                     'complex(dp)':'complex',
                     'real(qp)':'float128'}

   T_REAL_A = 6
   T_INTEGER_A = 5
   T_CHAR_A = 10
   T_COMPLEX_A = 7
   
   fortran_type_code = {
      'd': T_REAL_A,
      'i': T_INTEGER_A,
      'S': T_CHAR_A,
      'complex': T_COMPLEX_A
      }

   max_type_len = max(map(len,numpy_type_map.values()))

   subnames = [x.name for x in subts+functs]
   for t in mod.types:

      if filtertypes is not None and not strip_type(t.name) in filtertypes: continue

      spec[shortname.lower()]['types'][t.name] = {'doc': '\n'.join(t.doc), 'elements':{}}
      thisdoc = spec[shortname.lower()]['types'][t.name]['elements']
      for el in t.elements:

          uses = set()
          uses.add(type_map[t.name.lower()].lower())

          if filtertypes is not None and el.type.startswith('type') and not strip_type(el.type) in filtertypes:
             print 'Omitting element %s as of unwrapped type %s' % (el.name, el.type)
             continue

          #f2py misparses arguments that are called "type"
          name = el.name
          if el.name == 'type': name = 'thetype'
          mytype = el.type
          if el.type.startswith('character'):
              # change from '(len=*)' or '(*)' syntax to *(*) syntax
              mytype = 'character*(*)'

          attributes = el.attributes[:]
          dim_list = filter(lambda x: x.startswith('dimension'), attributes)

          if 'pointer' in attributes and dim_list != []: continue
          if mytype.lower() == 'type(c_ptr)': continue

          thisdoc[el.name] = {'doc': '\n'.join(el.doc), 'type': el.type, 'attributes': attributes}

          # If it's a proper array (not a pointer) let's write an __array__ routine 
          # which returns shape and data location (suitable for constructing
          # a numpy array that shares the same data)

          if mytype.startswith('type'):
              typename = strip_type(mytype)
          elif mytype in numpy_type_map:
              typename = numpy_type_map[mytype]
          else:
              typename = mytype

          if  dim_list != []:
              if mytype.startswith('type'):
                 tname = strip_type(t.name)
                 uses.add(type_map[typename].lower())

                 if len(dim_list) != 1: continue

                 println('subroutine %s%s__array_getitem__%s(this, i, the%s)' % (prefix, t.name,name,name))
                 indent += 3
                 #println('use %s, only: imported_%s => %s' % (type_map[t.name.lower()], t.name, t.name))
                 for line in kindlines:
                    if not any([line.startswith('use %s' % umod) for umod in uses]):
                       println(line)
                 for umod in uses:
                    println('use %s' % umod)
                 println('implicit none')
                 println('type %s_ptr_type' %  tname)
                 println('type(%s), pointer :: p' % tname)
                 println('end type %s_ptr_type' % tname)
                 println('integer, intent(in) :: this(%d)' % sizeof_fortran_t)
                 println('type(%s_ptr_type) :: this_ptr' % t.name)

                 println('type %s_ptr_type' %  strip_type(mytype))
                 println('type(%s), pointer :: p' % strip_type(mytype))
                 println('end type %s_ptr_type' % strip_type(mytype))
                 println('integer, intent(in) :: i')
                 println('integer, intent(out) :: the%s(%d)' % (name, sizeof_fortran_t))
                 println('type(%s_ptr_type) :: the%s_ptr' % (typename, name))
                 println()
                 println('this_ptr = transfer(this, this_ptr)')

                 if 'allocatable' in el.attributes:
                     println('if (allocated(this_ptr%%p%%%s)) then' % el.name)
                     indent += 3

                 println('if (i < 1 .or. i > size(this_ptr%%p%%%s)) then' % el.name)
                 indent += 3
                 println('call system_abort("array index out of range")')
                 indent -= 3
                 println('else')
                 indent += 3
                 println('the%s_ptr%%p => this_ptr%%p%%%s(i)' % (name, el.name))
                 println('the%s = transfer(the%s_ptr,the%s)' % (name, name, name))
                 indent -= 3
                 println('endif')

                 if 'allocatable' in el.attributes:
                     indent -= 3
                     println('else')
                     indent += 3
                     println('call system_abort("derived type array not allocated")')
                     indent -= 3
                     println('end if')

                 indent -= 3
                 println('end subroutine %s%s__array_getitem__%s' % (prefix, t.name, name))
                 println()
                 thisdoc[el.name]['array_getitem'] = '%s%s__array_getitem__%s' % (prefix, t.name.lower(), name.lower())

                 println('subroutine %s%s__array_setitem__%s(this, i, the%s)' % (prefix, t.name, name, name))
                 indent += 3
                 for line in kindlines:
                    if not any([line.startswith('use %s' % umod) for umod in uses]):
                       println(line)
                 for umod in uses:
                    println('use %s' % umod)                 

                 println('implicit none')
                 println('type %s_ptr_type' %  tname)
                 println('type(%s), pointer :: p' % tname)
                 println('end type %s_ptr_type' % tname)
                 println('integer, intent(in) :: this(%d)' % sizeof_fortran_t)
                 println('type(%s_ptr_type) :: this_ptr' % t.name)

                 println('type %s_ptr_type' %  strip_type(mytype))
                 println('type(%s), pointer :: p' % strip_type(mytype))
                 println('end type %s_ptr_type' % strip_type(mytype))
                 println('integer, intent(in) :: i')
                 println('integer, intent(in) :: the%s(%d)' % (name, sizeof_fortran_t))
                 println('type(%s_ptr_type) :: the%s_ptr' % (typename, name))
                 println()
                 println('this_ptr = transfer(this, this_ptr)')

                 if 'allocatable' in el.attributes:
                     println('if (allocated(this_ptr%%p%%%s)) then' % el.name)
                     indent += 3

                 println('if (i < 1 .or. i > size(this_ptr%%p%%%s)) then' % el.name)
                 indent += 3
                 println('call system_abort("array index out of range")')
                 indent -= 3
                 println('else')
                 indent += 3
                 println('the%s_ptr = transfer(the%s,the%s_ptr)' % (name, name, name))
                 println('this_ptr%%p%%%s(i) = the%s_ptr%%p' % (name, el.name))
                 indent -= 3
                 println('endif')

                 if 'allocatable' in el.attributes:
                     indent -= 3
                     println('else')
                     indent += 3
                     println('call system_abort("derived type array not allocated")')
                     indent -= 3
                     println('end if')

                 indent -= 3
                 println('end subroutine %s%s__array_setitem__%s' % (prefix, t.name, name))
                 println()
                 thisdoc[el.name]['array_setitem'] = '%s%s__array_setitem__%s' % (prefix, t.name.lower(), name.lower())

                 println('subroutine %s%s__array_len__%s(this, n)' % (prefix, t.name, name))
                 indent += 3
                 for line in kindlines:
                    if not any([line.startswith('use %s' % umod) for umod in uses]):
                       println(line)
                 for umod in uses:
                    println('use %s' % umod)                 

                 println('implicit none')
                 println('type %s_ptr_type' %  tname)
                 println('type(%s), pointer :: p' % tname)
                 println('end type %s_ptr_type' % tname)
                 println('integer, intent(in) :: this(%d)' % sizeof_fortran_t)
                 println('integer, intent(out) :: n')
                 println('type(%s_ptr_type) :: this_ptr' % t.name)
                 println()
                 println('this_ptr = transfer(this, this_ptr)')

                 if 'allocatable' in el.attributes:
                     println('if (allocated(this_ptr%%p%%%s)) then' % el.name)
                     indent += 3

                 println('n = size(this_ptr%%p%%%s)' % el.name)

                 if 'allocatable' in el.attributes:
                     indent -= 3
                     println('else')
                     indent += 3
                     println('n = 0')
                     indent -= 3
                     println('end if')

                 indent -= 3
                 println('end subroutine %s%s__array_len__%s' % (prefix, t.name, name))
                 println()
                 thisdoc[el.name]['array_len'] = '%s%s__array_len__%s' % (prefix, t.name.lower(), name.lower())                 
                 

              else:
                 tname = strip_type(t.name)

                 println('subroutine %s%s__array__%s(this, nd, dtype, dshape, dloc)' % (prefix, t.name,name))
                 indent += 3
                 #println('use %s, only: imported_%s => %s' % (type_map[t.name.lower()], t.name, t.name))
                 for umod in uses:
                    println('use %s' % umod)
                 println('implicit none')
                 println('type %s_ptr_type' %  tname)
                 println('type(%s), pointer :: p' % tname)
                 println('end type %s_ptr_type' % tname)
                 println('integer, intent(in) :: this(%d)' % sizeof_fortran_t)
                 println('type(%s_ptr_type) :: this_ptr' % t.name)
                 println('integer, intent(out) :: nd')
                 println('integer, intent(out) :: dtype')
                 try:
                     rank = dim_list[0].count(',')+1
                     if mytype.startswith('character'): rank += 1
                 except ValueError:
                     rank = 1
                 println('integer, dimension(10), intent(out) :: dshape')
                 println('integer*%d, intent(out) :: dloc' % numpy.dtype('O').itemsize)
                 println()
                 println('nd = %d' % rank)
                 println('dtype = %s' % fortran_type_code[typename])
                 println('this_ptr = transfer(this, this_ptr)')
                 if 'allocatable' in el.attributes:
                     println('if (allocated(this_ptr%%p%%%s)) then' % el.name)
                     indent += 3
                 if mytype.startswith('character'):
                     first = ','.join(['1' for i in range(rank-1)])
                     println('dshape(1:%d) = (/len(this_ptr%%p%%%s(%s)), shape(this_ptr%%p%%%s)/)' % (rank, el.name, first, el.name))
                 else:
                    println('dshape(1:%d) = shape(this_ptr%%p%%%s)' % (rank, el.name))
                 println('dloc = loc(this_ptr%%p%%%s)' % el.name)
                 if 'allocatable' in el.attributes:
                     indent -= 3
                     println('else')
                     indent += 3
                     println('dloc = 0')
                     indent -= 3
                     println('end if')

                 indent -= 3
                 println('end subroutine %s%s__array__%s' % (prefix, t.name, name))
                 thisdoc[el.name]['array'] = '%s%s__array__%s' % (prefix, t.name.lower(), name.lower())
                 println()

          # For scalars write get/set routines
          else:
              if mytype.startswith('type'):
                  typename = strip_type(mytype)
                  uses.add(type_map[typename].lower())
              elif mytype in numpy_type_map:
                  typename = numpy_type_map[mytype]
              else:
                  typename = mytype

              println('subroutine %s%s__get__%s(this, the%s)' % (prefix, t.name, name, name))
              indent += 3
              for line in kindlines:
                 if not any([line.startswith('use %s' % umod) for umod in uses]):
                    println(line)
              for umod in uses:
                 println('use %s' % umod)
              println('implicit none')
              println('type %s_ptr_type' %  t.name)
              println('type(%s), pointer :: p' % t.name)
              println('end type %s_ptr_type' % t.name)
              if mytype.startswith('type') and strip_type(mytype).lower() != t.name.lower():
                 println('type %s_ptr_type' %  strip_type(mytype))
                 println('type(%s), pointer :: p' % strip_type(mytype))
                 println('end type %s_ptr_type' % strip_type(mytype))
              println('integer, intent(in)   :: this(%d)' % sizeof_fortran_t)
              println('type(%s_ptr_type) :: this_ptr' % t.name)


              if el.type.startswith('type'):

                  # For derived types elements, treat as opaque reference
                  println('integer, intent(out) :: the%s(%d)' % (name, sizeof_fortran_t))
                  println('type(%s_ptr_type) :: the%s_ptr' % (typename, name))
                  println()
                  println('this_ptr = transfer(this, this_ptr)')
                  println('the%s_ptr%%p => this_ptr%%p%%%s' % (name, el.name))
                  println('the%s = transfer(the%s_ptr,the%s)' % (name, name, name))

              else:
                  # Return by value
                  if 'pointer' in attributes: attributes.remove('pointer')

                  if el.type.startswith('character'):

                      # change from '(len=*)' or '(*)' syntax to *(*) syntax
                      try:
                          lind = el.type.index('(')
                          rind = el.type.rindex(')')
                          mytype = el.type[:lind]+'*'+el.type[lind:rind+1].replace('len=','')

                          # Try to get length of string arguments
                          if not mytype[11:-1] == '*' and not all([x in '0123456789' for x in mytype[11:-1]]):
                              try:
                                  mytype = 'character*(%s)' % string_lengths[mytype[11:-1].lower()]
                              except KeyError:
                                  mytype = 'character*(%s)' % string_lengths['default']

                          # Default string length for intent(out) strings 
                          if mytype[11:-1] == '*' and 'intent(out)' in attributes:
                              mytype = 'character*(%s)' % string_lengths['default']

                      except ValueError:
                          pass

                  if attributes != []:
                      println('%s, %s, intent(out) :: the%s' % (mytype, ','.join(attributes), name))
                  else:
                      println('%s, intent(out) :: the%s' % (mytype, name))
                  println()
                  println('this_ptr = transfer(this, this_ptr)')
                  println('the%s = this_ptr%%p%%%s' % (name, el.name))

              indent -= 3
              println('end subroutine %s%s__get__%s' % (prefix, t.name, name))
              thisdoc[el.name]['get'] = '%s%s__get__%s' % (prefix, t.name.lower(), name.lower())

              println()

              println('subroutine %s%s__set__%s(this, the%s)' % (prefix, t.name, name, name))
              indent += 3
              for line in kindlines:
                 if not any([line.startswith('use %s' % umod) for umod in uses]):
                    println(line)
              for umod in uses:
                 println('use %s' % umod)                 
              println('implicit none')
              println('type %s_ptr_type' %  t.name)
              println('type(%s), pointer :: p' % t.name)
              println('end type %s_ptr_type' % t.name)
              if mytype.startswith('type') and strip_type(mytype).lower() != t.name.lower():
                 println('type %s_ptr_type' %  strip_type(mytype))
                 println('type(%s), pointer :: p' % strip_type(mytype))
                 println('end type %s_ptr_type' % strip_type(mytype))
              println('integer, intent(in)   :: this(%d)' % sizeof_fortran_t)
              println('type(%s_ptr_type) :: this_ptr' % t.name)
              attributes = el.attributes[:]

              if el.type.startswith('type'):
                  # Set by reference
                  println('integer, intent(in) :: the%s(%d)' % (name, sizeof_fortran_t))
                  println('type(%s_ptr_type) :: the%s_ptr' % (typename,name))
                  println()
                  println('this_ptr = transfer(this, this_ptr)')
                  println('the%s_ptr = transfer(the%s,the%s_ptr)' % (name, name, name))
                  println('this_ptr%%p%%%s = the%s_ptr%%p' % (el.name, name))
                  

              else:
                  # Set by value
                  if attributes != []:
                      println('%s, %s, intent(in) :: the%s' % (mytype, ','.join(attributes), name))
                  else:
                      println('%s, intent(in) :: the%s' % (mytype, name))
                  println()
                  println('this_ptr = transfer(this, this_ptr)')
                  println('this_ptr%%p%%%s = the%s' % (el.name, name))

              indent -= 3
              println('end subroutine %s%s__set__%s' % (prefix, t.name, name))
              thisdoc[el.name]['set'] = '%s%s__set__%s' % (prefix, t.name.lower(), name.lower())
              println()

   println()

   return spec


def filter_source_file(infile, outfile, is_wrapper=False):
   start_sub = re.compile(r'^(double precision\s+)?\s*(recursive\s+)?(subroutine|function)\s+(\w*)')
   end_sub   = re.compile(r'^\s*end\s+(subroutine|function)\s*$')

   in_interface = False
   seen_interfaces = {}
   interface_lines = []

   for line in infile.readlines():
      sline = line.strip()
      start_match = start_sub.match(sline)
      end_match = end_sub.match(sline)
      if start_match:
         last_subroutine = start_match.group(4)
         outfile.write(line)
      elif end_match:
         outfile.write('%s %s\n' % (sline, last_subroutine))
      elif is_wrapper and sline.startswith('interface'):
         in_interface = True
         interface_lines = []
         outfile.write(line)
      elif is_wrapper and (sline.startswith('end interface') or sline.startswith('endinterface')):
         in_interface = False
         outfile.write(line)
         for interface_line in interface_lines:
            if interface_line.startswith('#'):
               outfile.write(interface_line)
            else:
               intfs = [x.strip() for x in interface_line.split(',')]
               for intf in intfs:
                  if intf not in seen_interfaces.keys():
                     outfile.write('public :: %s\n' % intf)
                  seen_interfaces[intf] = 1
      elif is_wrapper and in_interface:
         outfile.write(line)
         if sline.startswith('#'):
            interface_lines.append(line)
         elif sline.startswith('module procedure'):
            sline = sline[len('module procedure'):]
            if '!' in sline: sline = sline[:sline.index('!')]
            interface_lines.append(sline)
      elif is_wrapper and (sline.startswith('private::') or sline.startswith('private ::')):
         outfile.write('!%s' % line)
      else:
         outfile.write(line)

        
def find_public_symbols(file):
   """Return a list of publically exported symbols in the Fortran module with filename 'file'.
      If the module defaults to public access, we return None."""
   data = open(file).read()
   data = data.replace('&\n', '')
   lines = [line.strip() for line in data.split('\n') ]

   if not 'private' in lines:
      print 'module %s defaults to public access' % file
      return None

   public_lines  = [line for line in lines if line.startswith('public')]

   res = []
   for line in public_lines:
      line = line.replace('public', '')
      line = line.replace('::', '')
      line = line.strip()
      res.extend([field.strip() for field in line.split(',')])
   return res
                

