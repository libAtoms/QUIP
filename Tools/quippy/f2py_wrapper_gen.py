from f90doc import *

def wrap_mod(mod, type_map, f2py_docs, out=None, kindlines=[]):
   """Given an f90doc.C_module class 'mod', write a F90 wrapper file suitable for f2py to 'out'."""

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
          return False

      return True

   # Also skip allocatable and pointer array arguments
   # (and return type if sub is a function)
   def no_allocatables_or_pointers(sub):
      for arg in sub.arguments:
          # FIXME: this skips scalar pointer args too
          if 'allocatable' in arg.attributes or 'pointer' in arg.attributes:
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
          dims = filter(lambda x: x.startswith('dimension'), arg.attributes)
          if arg.type.startswith('complex') and len(dims) == 0: return False

      return True


   def no_c_ptr(sub):
      for arg in sub.arguments:
          if arg.type.lower() == 'type(c_ptr)': return False

      return True

   debug(mod.name)
   shortname = mod.name[:mod.name.index('_module')]

   if out is None:
      out = sys.stdout #open('%s.f90' %shortname, 'w')


   # Don't care about interfaces, so put subts and functs back
   # into flat lists. Copy interface doc comment to all members.
   subts = filter(no_callbacks, mod.subts)
   functs = filter(no_callbacks, mod.functs)
   for intf in mod.interfaces:
      thissubs = filter(no_callbacks, intf.subts)
      thisfuncts = filter(no_callbacks, intf.functs)
      subts.extend(thissubs)
      functs.extend(thisfuncts)


   subts = filter(no_allocatables_or_pointers, subts)
   functs = filter(no_allocatables_or_pointers, functs)

   subts = filter(no_complex_scalars, subts)
   functs = filter(no_complex_scalars, functs)

   subts = filter(no_c_ptr, subts)
   functs = filter(no_c_ptr, functs)

   debug('%s: %d subs' % (mod.name, len(subts+functs)))

   names = [x.name for x in subts+functs]
   if shortname in names: shortname += '_'

   indent = 0
   println('module',shortname)
   indent += 3

   use_list = [ 'my_%s => %s' % (x.name,x.name) for x in subts+functs]

   if len(subts+functs) != 0:
      println('use %s, only: &' % mod.name)
      indent += 3
      while use_list:
          take = use_list[:min(1,len(use_list))]
          del use_list[:min(1,len(use_list))]
          term = ', &'
          if use_list == []:
              term = ''
          println(', '.join(take)+term)
      indent -= 3
      println()

   # Add uses clauses for types used in this module
   dep_types = []
   for sub in subts + functs:
      for arg in sub.arguments:
          if arg.type.startswith('type'):
              t = arg.type[arg.type.index('(')+1:arg.type.index(')')].lower()
              if t not in dep_types: dep_types.append(t)

   for t in mod.types:
      tname = t.name.lower()
      if not tname in dep_types: dep_types.append(tname)
      for el in t.elements:
          if el.type.startswith('type'):
              tname = el.type[el.type.index('(')+1:el.type.index(')')].lower()
              if not tname in dep_types: dep_types.append(tname)

   dep_mods = {}
   for dep in dep_types:
      modname = type_map[dep]
      if not modname in dep_mods:
          dep_mods[modname] = []
      dep_mods[modname].append(dep)

   for modname in dep_mods.keys():
      println('use %s, only: &' % modname)
      indent += 3
      println(','.join(['my_%s => %s' % (t, t) for t in dep_mods[modname]]))
      indent -= 3
   println()

   # Kind lines
   for line in kindlines:
      println(line)
   # Must have at least one item in module or f2py gives segfaults
   # Let's put the module name in, could be handy
   println('character*(%d), parameter :: module_name = "%s"' % (len(shortname), shortname))


   if len(subts+functs) > 0:
      println()
      println('contains')
      println()
      indent += 3

   f2py_docs[shortname.lower()] = {'doc': '\n'.join(mod.doc),
                                  'routines': {},
                                  'types': {},
                                  'interfaces': {},
                                  'parameters': []}

   # for some reason setting parameters in Fortran causes an AssertionError
   # when importing module. For now, let's just copy them into python dictionary
   for el in mod.elements:
      if 'parameter' in el.attributes:
          f2py_docs[shortname.lower()]['parameters'].append((el.name, el.type, el.attributes, el.value))

   for sub in subts + functs:

      # Add argument for return value, after last non-optional argument
      if hasattr(sub,'ret_val'):
          ret_val = sub.ret_val
          ret_val.name = 'ret_'+ret_val.name
          ret_val.attributes.append('intent(out)')
          ret_val.is_ret_val = True
          sub.arguments = ([ arg for arg in sub.arguments if not 'optional' in arg.attributes ] + [ret_val] +
                           [ arg for arg in sub.arguments if 'optional' in arg.attributes ])

      args = sub.arguments
      argnames = [x.name for x in args]
      newargnames = argnames[:]

      arglines = []
      optionals = []

      n_dummy = 0

      newname = sub.name
      # Ensure that methods start with the class name followed by an underscore
      #if len(args) > 0 and args[0].type.startswith('type') and args[0].name == 'this':
      #    if not sub.name.lower().startswith(args[0].type[5:-1].lower()+'_'):
      #        if sub.name.startswith('ds_'):
      #            newname = args[0].type[5:-1].lower()+'_'+sub.name[3:]
      #        else:
      #            newname = args[0].type[5:-1].lower()+'_'+sub.name
      #else:

      #br = sub.name.find('_')
      #if br != -1:
      #    oldtype = sub.name[:br]
      #    if len(args) > 0 and oldtype.lower() == args[0].type[5:-1].lower():
      #        basename = sub.name[br+1:]
      #    else:
      #        basename = sub.name
      #else:
      #    basename = sub.name

      #if len(typenames) < 8:
      #    newname = '__'.join(typenames)
      #else:
      #    newname = 'too_many_types_%s' % basename

      f2py_docs[shortname.lower()]['routines'][newname.lower()] = \
          {'doc': '\n'.join(sub.doc),'args':[]}
      thisdoc = f2py_docs[shortname.lower()]['routines'][newname.lower()]['args']

      # See if this routine is in any interfaces
      for intf in mod.interfaces:
          subnames = [x.name.lower() for x in intf.subts + intf.functs]
          if sub.name.lower() in subnames:
              if not intf.name.lower() in f2py_docs[shortname.lower()]['interfaces']:
                  f2py_docs[shortname.lower()]['interfaces'][intf.name.lower()] = \
                      {'doc': intf.doc,
                       'routines': []}

              f2py_docs[shortname.lower()]['interfaces'][intf.name.lower()]['routines'].append(sub.name.lower())

      allocates = []

      for arg in args:

          # Replace all type args with pointers
          if not hasattr(arg, 'attributes') and not hasattr(arg, 'type'):
              arglines.append('external %s' % arg.name)
              continue

          attributes = arg.attributes
          mytype = arg.type

          if 'optional' in attributes and 'intent(out)' in attributes: 
              attributes.remove('intent(out)')
              attributes.append('intent(inout)')

          if arg.type.startswith('type'):
              mytype = 'type(my_%s' % arg.type[arg.type.index('(')+1:]

              # Preserve original fortran intent
              fintent = [ a for a in attributes if a.startswith('intent')]
              if fintent != []: fintent = fintent[0].replace('intent', 'fintent')

              if ('intent(out)' in attributes or 
                  ((sub.name.lower().find('_initialise') != -1 or sub.name.lower().find('_allocate') != -1) \
                       and len(argnames) > 0 and argnames[0] == 'this' and arg.name == 'this')):
                  allocates.append(arg.name)
                  intent = 'intent(out)'
              else:
                  intent = 'intent(in)'

              arglines.append('!f2py integer*SIZEOF_VOID_PTR, %s :: %s' % (intent, arg.name))
              attributes = filter(lambda x:x.startswith('dimension') or
                                  x.startswith('allocatable') or x.startswith('optional'),attributes)
              attributes.append('pointer')
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
                  pass

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
                      if d.startswith('len'):
                          arglines.append('!f2py %s %s, dimension(%s) :: %s' % \
                                              (arg.type, 
                                               ','.join(filter(lambda a: not a.startswith('dimension'), attributes)), 
                                               d.replace('len','slen'), arg.name))
                      continue
                  ds[i] = ('n%d' % (n_dummy))
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
                  arglines.append('integer :: %s' % d)
                  if not 'intent(out)' in attributes:
                      arglines.append('!f2py intent(hide), depend(%s) :: %s = shape(%s,%d)' % (arg.name, d, arg.name, i))
                  else:
                      thisdoc.append({'name': d, 'doc': 'shape(%s,%d)' % (arg.name,i), 'type': 'integer', 'attributes':[]})

          if charflag is not None:
              newargnames.append(charflag)
              arglines.append('integer :: %s' % charflag)
              if not 'intent(out)' in attributes:
                  arglines.append('!f2py intent(hide), depend(%s) :: %s = slen(%s)' % (arg.name, charflag, arg.name))
              else:
                  thisdoc.append({'name':charflag,'doc': 'slen(%s)' % arg.name, 'type': 'integer', 'attributes':[]})


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

      println('subroutine %s(%s)' % (newname, join_and_split_lines(newargnames)))
      indent += 3
      for line in arglines:
          println(line)
      println()

      for var in allocates: println('allocate(%s)' % var)

      if hasattr(sub, 'ret_val'):
          argfilt = [ arg.name for arg in args if not (hasattr(arg, 'is_ret_val') and arg.is_ret_val) ]
          println('%s = my_%s(%s)' % (sub.ret_val.name, sub.name, join_and_split_lines(argfilt)))
      else:
          println('call my_%s(%s)' % (sub.name, join_and_split_lines(argnames)))

      if sub.name.lower().endswith('finalise') and len(argnames) > 0 and argnames[0] == 'this':
          println('deallocate(this)')

      println()

      indent -= 3
      println('end subroutine',newname)
      println()


   # add _get_<type> and _set_<type> methods

   numpy_type_map = {'real(dp)':'d','integer':'i','logical':'i','character*(*)':'S','complex(dp)':'complex'}
   max_type_len = max(map(len,numpy_type_map.values()))

   subnames = [x.name for x in subts+functs]
   for t in mod.types:

      f2py_docs[shortname.lower()]['types'][t.name] = {'doc': '\n'.join(t.doc), 'elements':{}}
      thisdoc = f2py_docs[shortname.lower()]['types'][t.name]['elements']
      for el in t.elements:

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
              typename = mytype[mytype.index('(')+1:mytype.index(')')]
          elif mytype in numpy_type_map:
              typename = numpy_type_map[mytype]
          else:
              typename = mytype

          if  dim_list != []:
              if mytype.startswith('type'): continue

              println('subroutine %s__array__%s(this, dtype, dshape, dloc)' % (t.name,name))
              indent += 3
              println('!f2py integer*SIZEOF_VOID_PTR, intent(in) :: this')
              println('type(my_%s), pointer, intent(in) :: this' % t.name)
              println('character(%d), intent(out) :: dtype' % max_type_len)
              try:
                  rank = dim_list[0].count(',')+1
                  if mytype.startswith('character'): rank += 1
              except ValueError:
                  rank = 1
              println('integer, dimension(%d), intent(out) :: dshape' % rank)
              println('integer*SIZEOF_VOID_PTR, intent(out) :: dloc')
              println()
              println('dtype = "%s"' % typename)
              if 'allocatable' in el.attributes:
                  println('if (allocated(this%%%s)) then' % el.name)
                  indent += 3
              if mytype.startswith('character'):
                  first = ','.join(['1' for i in range(rank-1)])
                  println('dshape = (/len(this%%%s(%s)), shape(this%%%s)/)' % (el.name, first, el.name))
              else:
                  println('dshape = shape(this%%%s)' % el.name)
              println('dloc = loc(this%%%s)' % el.name)
              if 'allocatable' in el.attributes:
                  indent -= 3
                  println('else')
                  indent += 3
                  println('dshape = (/%s/)' % ','.join(['0' for x in range(rank)]))
                  println('dloc   = 0')
                  indent -= 3
                  println('end if')

              indent -= 3
              println('end subroutine %s__array__%s' % (t.name, name))
              thisdoc[el.name]['array'] = '%s__array__%s' % (t.name.lower(), name.lower())
              println()

          # For scalars write get/set routines
          else:
              if mytype.startswith('type'):
                  typename = mytype[mytype.index('(')+1:mytype.index(')')]
              elif mytype in numpy_type_map:
                  typename = numpy_type_map[mytype]
              else:
                  typename = mytype

              println('subroutine %s__get__%s(this, the%s)' % (t.name, name, name))
              indent += 3
              println('!f2py integer*SIZEOF_VOID_PTR, intent(in) :: this')
              println('type(my_%s), pointer, intent(in) :: this' % t.name)


              if el.type.startswith('type'):
                  # For derived types elements, just treat as a pointer
                  println('!f2py integer*SIZEOF_VOID_PTR, intent(out) :: the%s' % name)
                  #if dim_list != []:
                  #    #println('type(my_%s, pointer, %s, intent(out) :: the%s' % (el.type[5:], dim_list[0], name))
                  #    println('integer, %s, intent(out) :: the%s' % (el.type[5:], dim_list[0], name))
                  #else:
                      #println('type(my_%s, pointer, intent(out) :: the%s' % (el.type[5:], name))
                  println('integer*SIZEOF_VOID_PTR, intent(out) :: the%s' % name)
                  println()
                  println('the%s = loc(this%%%s)' % (name, el.name))

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
                  println('the%s = this%%%s' % (name, el.name))

              indent -= 3
              println('end subroutine %s__get__%s' % (t.name, name))
              thisdoc[el.name]['get'] = '%s__get__%s' % (t.name.lower(), name.lower())

              println()

              println('subroutine %s__set__%s(this, the%s)' % (t.name, name, name))
              indent += 3
              println('!f2py integer*SIZEOF_VOID_PTR, intent(in) :: this')
              println('type(my_%s), pointer, intent(inout) :: this' % t.name)
              attributes = el.attributes[:]

              if el.type.startswith('type'):
                  # Set by reference
                  println('!f2py integer*SIZEOF_VOID_PTR, intent(in) :: the%s' % name)
                  println('type(my_%s, pointer, intent(in) :: the%s' % (el.type[el.type.index('(')+1:], name))
                  println()
                  println('this%%%s = the%s' % (el.name, name))

              else:
                  # Set by value
                  if attributes != []:
                      println('%s, %s, intent(in) :: the%s' % (mytype, ','.join(attributes), name))
                  else:
                      println('%s, intent(in) :: the%s' % (mytype, name))
                  println()
                  println('this%%%s = the%s' % (el.name, name))

              indent -= 3
              println('end subroutine %s__set__%s' % (t.name, name))
              thisdoc[el.name]['set'] = '%s__set__%s' % (t.name.lower(), name.lower())
              println()

   indent -= 6
   println('end module',shortname)
   println()

        


