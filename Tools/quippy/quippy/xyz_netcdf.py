from quippy import (atoms_reader, AtomsReaders, AtomsWriters, TABLE_STRING_LENGTH,
                    PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL,
                    T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                    T_CHAR, T_LOGICAL, T_INTEGER_A,
                    T_REAL_A, T_COMPLEX_A, T_CHAR_A, T_LOGICAL_A)

from atomslist import *
from ordereddict import *
from quippy import Dictionary
from farray import *
import logging, StringIO
from math import pi

def make_lattice(a, b=None, c=None, alpha=pi/2.0, beta=pi/2.0, gamma=pi/2.0):
   "Form 3x3 lattice from cell lengths (a,b,c) and angles (alpha,beta,gamma) in radians"

   if b is None: b = a
   if c is None: c = a
   
   lattice = fzeros((3,3),'d')

   cos_alpha = numpy.cos(alpha); cos2_alpha = cos_alpha*cos_alpha
   cos_beta  = numpy.cos(beta);  cos2_beta  = cos_beta *cos_beta
   cos_gamma = numpy.cos(gamma); cos2_gamma = cos_gamma*cos_gamma
   sin_gamma = numpy.sin(gamma); sin2_gamma = sin_gamma*sin_gamma
   
   lattice[1,1] = a

   lattice[1,2] = b * cos_gamma
   lattice[2,2] = b * sin_gamma

   lattice[1,3] = c * cos_beta
   lattice[2,3] = c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma
   lattice[3,3] = c * numpy.sqrt(1.0 - (cos2_alpha + cos2_beta - 2.0*cos_alpha*cos_beta*cos_gamma)/ sin2_gamma)

   return lattice

def get_lattice_params(lattice):
   """Return 6-tuple of cell lengths and angles a,b,c,alpha,beta,gamma"""
   from math import acos
   from numpy import dot
   a_1 = lattice[:,1]
   a_2 = lattice[:,2]
   a_3 = lattice[:,3]

   return (a_1.norm(), a_2.norm(), a_3.norm(),
           acos(dot(a_2,a_3)/(a_2.norm()*a_3.norm())),
           acos(dot(a_3,a_1)/(a_3.norm()*a_1.norm())),
           acos(dot(a_1,a_2)/(a_1.norm()*a_2.norm())))


def _props_dtype(props):
   "Return a record array dtype for the specified properties (default all)"

   result = []
   fmt_map = {PROPERTY_REAL:'d',
              PROPERTY_INT: 'i',
              PROPERTY_STR: 'S' + str(TABLE_STRING_LENGTH),
              PROPERTY_LOGICAL: 'bool'}

   for prop in props:
      ptype, col_start, col_stop = props[prop]
      if col_start == col_stop:
         result.append((prop,fmt_map[ptype]))
      else:
         for c in range(col_stop-col_start+1):
            result.append((prop+str(c),fmt_map[ptype]))

   return numpy.dtype(result)


@atoms_reader('pupyxyz', False)
def PuPyXYZReader(xyz):
   "Read from extended XYZ filename, string or open file."
   from quippy import Table

   def parse_properties(prop_str):
      """Parse a property description string in the format
      name:ptype:cols:...  and return a 5-tuple of an OrderedDict
      mapping property names to (ptype,col_start,col_end) tuples and the number of
      int, real, string and logical properties"""

      from quippy import (PROPERTY_INT,
                          PROPERTY_REAL,
                          PROPERTY_STR,
                          PROPERTY_LOGICAL)

      pmap = {'I': PROPERTY_INT,
              'R': PROPERTY_REAL,
              'S': PROPERTY_STR,
              'L': PROPERTY_LOGICAL}

      properties = OrderedDict()
      if prop_str.startswith('pos:R:3'):
         prop_str = 'species:S:1:'+prop_str

      fields = prop_str.split(':')

      n_real = 0
      n_int  = 0
      n_str = 0
      n_logical = 0
      for name, ptype, cols in zip(fields[::3], fields[1::3], map(int,fields[2::3])):
         if ptype == 'R':
            start_col = n_real
            n_real = n_real + cols
            end_col = n_real
         elif ptype == 'I':
            start_col = n_int
            n_int = n_int + cols
            end_col = n_int
         elif ptype == 'S':
            start_col = n_str
            n_str = n_str + cols
            end_col = n_str
         elif ptype == 'L':
            start_col = n_logical
            n_logical = n_logical + cols
            end_col = n_logical
         else:
            raise ValueError('Unknown property type: '+ptype)

         properties[name] = (pmap[ptype],start_col+1,end_col)

      return properties, n_int, n_real, n_str, n_logical


   def table_update_from_recarray(self, props, data):
      for prop in props:
         ptype, col_start, col_stop = props[prop]
         if ptype == PROPERTY_REAL:
            if col_start == col_stop:
               self.real[col_start,:] = data[prop]
            else:
               for c in range(col_stop-col_start+1):
                  self.real[col_start+c,:] = data[prop+str(c)]
         elif ptype == PROPERTY_INT:
            if col_start == col_stop:
               self.int[col_start,:] = data[prop]
            else:
               for c in range(col_stop-col_start+1):
                  self.int[col_start+c,:] = data[prop+str(c)]
         elif ptype == PROPERTY_STR:
            if col_start == col_stop:
               self.str[:,col_start,:] = padded_str_array(data[prop], TABLE_STRING_LENGTH)
            else:
               for c in range(col_stop-col_start+1):
                  self.str[:,col_start+c,:] = padded_str_array(data[prop+str(c)], TABLE_STRING_LENGTH)
         elif ptype == PROPERTY_LOGICAL:
            if col_start == col_stop:
               self.logical[col_start,:] = data[prop]
            else:
               for c in range(col_stop-col_start+1):
                  self.logical[col_start+c,:] = data[prop+str(c)]
         else: 
            raise ValueError('Bad property type : '+str(ptype))

   def _getconv(dtype):
       typ = dtype.type
       if issubclass(typ, numpy.bool_):
          return lambda x: {'T' : True, 'F' : False, 'True' : True, 'False' : False}.get(x)
       if issubclass(typ, numpy.integer):
           return int
       elif issubclass(typ, numpy.floating):
           return float
       elif issubclass(typ, numpy.complex):
           return complex
       else:
           return str

   opened = False
   if type(xyz) == type(''):
      if '\n' in xyz and xyz[:xyz.index('\n')].isdigit():
         # string looks like an embedded XYZ file
         xyz = iter(xyz.split('\n'))
      else:
         xyz = open(xyz,'r')
         opened = True
   else:
      xyz = iter(xyz)

   while True:
      line = xyz.next()
      if not line: raise StopIteration

      n = int(line.strip())
      comment = (xyz.next()).strip()

      # Parse comment line
      params = Dictionary(comment)

      if not 'Properties' in params:
         raise ValueError('Properties missing from comment line')

      properties, n_int, n_real, n_str, n_logical = parse_properties(params['Properties'])
      del params['Properties']

      # Get lattice
      if not 'Lattice' in params:
         raise ValueError('No lattice found in xyz file')

      lattice = numpy.reshape(params['Lattice'], (3,3), order='F')
      del params['Lattice']

      props_dtype = _props_dtype(properties)

      converters = [_getconv(props_dtype.fields[name][0]) \
                    for name in props_dtype.names]

      X = []
      for i,line in enumerate(xyz):
         vals = line.split()
         row = tuple([converters[j](val) for j,val in enumerate(vals)])
         X.append(row)
         if i == n-1: break # Only read self.n lines

      try:
         data = numpy.array(X,props_dtype)
      except TypeError:
         raise IOError('Badly formatted data, or end of file reached before end of frame')

      tab = Table(n_int,n_real,n_str,n_logical)
      tab.append(blank_rows=n)
      table_update_from_recarray(tab, properties, data)

      yield Atoms(n=n,lattice=lattice, properties=properties, params=params, data=tab)

   if opened: xyz.close()


class PuPyXYZWriter(object):
   "Write atoms in extended XYZ format. xyz can be a filename or open file"

   def __init__(self, xyz):
      self.opened = False
      self.string = False
      if type(xyz) == type(''):
         if xyz == 'stdout':
            self.xyz = sys.stdout
         elif xyz == 'string':
            self.xyz = StringIO.StringIO()
            self.string = True
         else:
            self.xyz = open(xyz, 'w')
            self.opened = True
      else:
         self.xyz = xyz

   def write(self, at, properties=None):

      def _getfmt(dtype):
         typ = dtype.type
         if issubclass(typ, numpy.bool_):
            return ' %.1s '
         if issubclass(typ, numpy.integer):
            return '%8d '
         elif issubclass(typ, numpy.floating):
            return '%16.8f'
         elif issubclass(typ, numpy.complex):
            return '(%f,%f)'
         else:
            return '%s '

      def comment(self, properties=None):
         if properties is None:
            props = self.properties.keys()
         else:
            props = properties

         pmap = { PROPERTY_INT: 'I',
                  PROPERTY_REAL: 'R',
                  PROPERTY_STR: 'S',
                  PROPERTY_LOGICAL: 'L'}

         lattice_str = 'Lattice="' + ' '.join(map(str, numpy.reshape(self.lattice,9,order='F'))) + '"'

         props_str = 'Properties=' + ':'.join(map(':'.join,
                                                  zip(props,
                                                      [pmap[self.properties[k][1]] for k in props],
                                                      [str(self.properties[k][3]-self.properties[k][2]+1) for k in props])))

         return lattice_str+' '+props_str+' '+str(self.params)


      def table_to_recarray(self, props):
         # Create empty record array with correct dtype
         data = numpy.zeros(self.n, _props_dtype(props))

         # Copy cols from self.real and self.int into data recarray
         for prop in props:
            ptype, col_start, col_stop = props[prop]
            if ptype == PROPERTY_REAL:
               if col_start == col_stop:
                  data[prop] = self.real[col_start,:]
               else:
                  for c in range(col_stop-col_start+1):
                     data[prop+str(c)] = self.real[col_start+c,:]
            elif ptype == PROPERTY_INT:
               if col_start == col_stop:
                  data[prop] = self.int[col_start,:]
               else:
                  for c in range(col_stop-col_start+1):
                     data[prop+str(c)] = self.int[col_start+c,:]
            elif ptype == PROPERTY_STR:
               if col_start == col_stop:
                  data[prop] = self.str[col_start,:].stripstrings()
               else:
                  for c in range(col_stop-col_start+1):
                     data[prop+str(c)] = self.str[col_start+c,:].stripstrings()
            elif ptype == PROPERTY_LOGICAL:
               if col_start == col_stop:
                  data[prop] = self.logical[col_start,:]
               else:
                  for c in range(col_stop-col_start+1):
                     data[prop+str(c)] = self.logical[col_start+c,:]
            else:
               raise ValueError('Bad property type : '+str(ptype))

         return data

      if properties is None:
         props = at.properties.keys()

      if 'species' in props:
         i = props.index('species')
         props[0], props[i] = props[i], props[0]

      if 'pos' in props:
         i = props.index('pos')
         props[1], props[i] = props[i], props[1]

      species = getattr(at, props[0].lower())
      if species.shape != (TABLE_STRING_LENGTH, at.n) or species.dtype.kind != 'S':
         raise ValueError('First property must be species like')

      pos = getattr(at, props[1].lower())
      if pos.shape != (3, at.n) or pos.dtype.kind != 'f':
         raise ValueError('Second property must be position like')

      subprops = OrderedDict.frompairs([(k, at.properties[k]) for k in props])
      data = table_to_recarray(at.data, subprops)
      format = ''.join([_getfmt(data.dtype.fields[name][0]) for name in data.dtype.names])+'\n'


      self.xyz.write('%d\n' % at.n)
      self.xyz.write('%s\n' % comment(at, props))
      for i in range(at.n):
         self.xyz.write(format % tuple(data[i]))

      if self.string: return self.xyz.getvalue()

   def close(self):
      if self.opened: self.xyz.close()

try:
   from netCDF4 import Dataset
   netcdf_file = Dataset
except ImportError:
   from pupynere import netcdf_file
   logging.warning('netCDF4 not found. falling back on (slower) pupynere.')

def netcdf_dimlen(obj, name):
   """Return length of dimension 'name'. Works for both netCDF4 and pupynere."""
   n = obj.dimensions[name]
   try:
      return len(n)
   except TypeError:
      return n

try:
   from quippy import CInOutput, INPUT, OUTPUT, INOUT

   class CInOutputReader(object):
      """Class to read atoms from a CInOutput. Supports generator and random access via indexing."""

      lazy = True

      def __init__(self, source, frame=None, start=0, stop=None, step=1, zero=False):
         if isinstance(source, str):
            self.opened = True
            self.source = CInOutput(source, action=INPUT, append=False, zero=zero)
            try:
               self.netcdf_file = netcdf_file(source)
            except (RuntimeError, AssertionError):
               self.netcdf_file = None
         else:
            self.opened = False
            self.source = source
            self.netcdf_file = None

         if frame is not None:
            self.start = frame
            self.stop = frame+1
            self.step = 1
         else:
            self.start = start
            if stop is None:
               stop = len(self.source)
            self.stop = stop
            self.step = step

      def __iter__(self):
         for frame in range(self.start,self.stop,self.step):
            yield self.source[frame]

      def close(self):
         if self.opened: self.source.close()

      def __len__(self):
         return len(self.source)

      def __getitem__(self, idx):
         return self.source[idx]

      def __getattr__(self, name):
         if self.netcdf_file is not None:
            try:
               return self.netcdf_file.__getattr__(name)
            except AttributeError:
               try:
                  return farray(self.netcdf_file.variables[name][:])
               except KeyError:
                  raise AttributeError('Attribute %s not found' % name)
         else:
            raise AttributeError('Attribute %s not found' % name)


   AtomsReaders['xyz'] = AtomsReaders['nc'] = AtomsReaders[CInOutput] = CInOutputReader

   @atoms_reader('stdin', True)
   def CInOutputStdinReader(source='stdin'):
      assert source == 'stdin'
      source = CInOutput(source, action=INPUT)
      while True:
         try:
            yield source.read()
         except RuntimeError:
            break         

   class CInOutputWriter(object):
      """Class to write atoms sequentially to a CInOutput stream"""

      def __init__(self, dest, append=False, netcdf4=True):
         self.opened = False
         if isinstance(dest, str):
            self.opened = True
            self.dest = CInOutput(dest, action=OUTPUT, append=append, netcdf4=netcdf4)
         else:
            self.dest = dest

      def write(self, at, frame=None, properties=None):
         self.dest.write(at, frame=frame, properties=properties)

      def close(self):
         self.dest.close()

   AtomsWriters['xyz'] = AtomsWriters['nc'] = AtomsWriters[CInOutput] = CInOutputWriter

   
except ImportError:
   logging.warning('CInOutput not found - falling back on (slower) pure python I/O')
   AtomsReaders['xyz'] = PuPyXYZReader
   AtomsWrtiters['xyz'] = PuPyXYZWriter

AtomsReaders['pupyxyz'] = AtomsReaders['string'] = PuPyXYZReader
AtomsWriters['pupyxyz'] = AtomsWriters['string'] = AtomsWriters['stdout'] = PuPyXYZWriter

@atoms_reader(netcdf_file, True)
@atoms_reader('nc', True)
def NetCDFReader(source, frame=None, start=0, stop=None, step=1):

   opened = False
   if isinstance(source, str):
      opened = True
      source = netcdf_file(source)

   from quippy import Atoms

   DEG_TO_RAD = pi/180.0

   remap_names = {'coordinates': 'pos',
                  'velocities': 'velo',
                  'cell_lengths': None,
                  'cell_angles': None}

   prop_type_to_value = {PROPERTY_INT: 0,
                         PROPERTY_REAL: 0.0,
                         PROPERTY_STR: "",
                         PROPERTY_LOGICAL: False}

   prop_dim_to_ncols = {('frame','atom','spatial'): 3,
                        ('frame','atom','label'): 1,
                        ('frame', 'atom'): 1}

   if frame is not None:
      start = frame
      stop = frame+1
      step = 1
   else:
      if stop is None:
         stop = source.variables['cell_lengths'].shape[0]
            
   for frame in range(start, stop, step):
      cl = source.variables['cell_lengths'][frame]
      ca = source.variables['cell_angles'][frame]
      lattice = make_lattice(cl[0],cl[1],cl[2],ca[0]*DEG_TO_RAD,ca[1]*DEG_TO_RAD,ca[2]*DEG_TO_RAD)

      at = Atoms(n=netcdf_dimlen(source, 'atom'), lattice=lattice)

      for name, var in source.variables.iteritems():
         name = remap_names.get(name, name)

         if name is None:
            continue

         if 'frame' in var.dimensions:
            if 'atom' in var.dimensions:
               # It's a property
               at.add_property(name, prop_type_to_value[var.type],
                              n_cols=prop_dim_to_ncols[var.dimensions])
               getattr(at, name.lower())[...] = var[frame].T
            else:
               # It's a param
               if var.dimensions == ('frame','string'):
                  # if it's a single string, join it and strip it
                  at.params[name] = ''.join(var[frame]).strip()
               else:
                  at.params[name] = var[frame].T

      yield at

   def close(self):
      if self.opened: source.close()

class NetCDFWriter(object):

   def __init__(self, dest, frame=None, append=False, netcdf4=True):
      self.opened = False
      self.frame = frame
      if isinstance(dest, str):
         self.opened = True
         try:
            FORMAT = {True:  'NETCDF4',
                      False: 'NETCDF3_CLASSIC'}
            MODE = {True:  'a',
                    False: 'w'}

            self.dest = netcdf_file(dest, MODE[append], format=FORMAT[netcdf4])
         except ValueError:
            self.dest = netcdf_file(dest, 'w')
      else:
         self.dest = dest

   def write(self, at):
      remap_names_rev = {'pos': 'coordinates',
                         'velocities': 'velo'}


      prop_type_ncols_to_dtype_dim = {(PROPERTY_INT,3):       ('i', ('frame','atom','spatial')),
                                      (PROPERTY_REAL,3):     ('d', ('frame','atom','spatial')),
                                      (PROPERTY_LOGICAL,3):  ('d', ('frame','atom','spatial')),
                                      (PROPERTY_INT,1):      ('i', ('frame', 'atom')),
                                      (PROPERTY_REAL, 1):    ('d', ('frame', 'atom')),
                                      (PROPERTY_LOGICAL, 1): ('i', ('frame', 'atom')),
                                      (PROPERTY_STR,1):      ('S1', ('frame','atom','label'))}


      param_type_to_dtype_dim = {T_INTEGER:   ('i', ('frame',)),
                                 T_REAL:      ('d', ('frame',)),
                                 T_CHAR:      ('S1', ('frame', 'string')),
                                 T_LOGICAL:   ('i', ('frame',)),
                                 T_INTEGER_A: ('i', ('frame', 'spatial')),
                                 T_REAL_A:    ('d', ('frame', 'spatial')),
                                 T_LOGICAL_A: ('i', ('frame', 'spatial'))}


      if self.dest.dimensions == {}:
         self.dest.createDimension('frame', None)
         self.dest.createDimension('spatial', 3)
         self.dest.createDimension('atom', at.n)
         self.dest.createDimension('cell_spatial', 3)
         self.dest.createDimension('cell_angular', 3)
         self.dest.createDimension('label', 10)
         self.dest.createDimension('string', 1024)

         self.dest.createVariable('spatial', 'S1', ('spatial',))
         self.dest.variables['spatial'][:] = list('xyz')

         self.dest.createVariable('cell_spatial', 'S1', ('cell_spatial',))
         self.dest.variables['cell_spatial'][:] = list('abc')

         self.dest.createVariable('cell_angular', 'S1', ('cell_angular', 'label'))
         self.dest.variables['cell_angular'][0,:] = list('alpha     ')
         self.dest.variables['cell_angular'][1,:] = list('beta      ')
         self.dest.variables['cell_angular'][2,:] = list('gamma     ')

         self.dest.createVariable('cell_lengths', 'd', ('frame', 'spatial'))
         self.dest.createVariable('cell_angles', 'd', ('frame', 'spatial'))

      if self.frame is None:
         self.frame = netcdf_dimlen(self.dest, 'frame')
         if self.frame is None: self.frame = 0


      assert at.n == netcdf_dimlen(self.dest, 'atom')

      a, b, c, alpha, beta, gamma = get_lattice_params(at.lattice)

      RAD_TO_DEG = 180.0/pi

      self.dest.variables['cell_lengths'][self.frame] = (a, b, c)
      self.dest.variables['cell_angles'][self.frame] = (alpha*RAD_TO_DEG, beta*RAD_TO_DEG, gamma*RAD_TO_DEG)

      for origname, (ptype, start, stop) in at.properties.iteritems():
         name = remap_names_rev.get(origname, origname)
         ncols = stop-start+1
         dtype, dims = prop_type_ncols_to_dtype_dim[(ptype, ncols)]

         if not name in self.dest.variables:
            self.dest.createVariable(name, dtype, dims)
            self.dest.variables[name].type = ptype

         assert self.dest.variables[name].dimensions == dims
         self.dest.variables[name][self.frame] = getattr(at, origname.lower()).T

      for name in at.params.keys():
         t, s = at.params.get_type_and_size()
         dtype, dims = param_type_to_dype_dims[(name, t)]

         if not name in self.dest.variables:
            self.dest.createVariable(name, dtype, dims)
            self.dest.variables[name].type = t

         assert self.dest.variables.type == t
         assert self.dest.variables.dimensions == dims

         self.dest.variables[name][self.frame] = at.params[name]

      self.frame += 1


   def close(self):
      if self.opened: self.dest.close()


AtomsWriters[netcdf_file] = NetCDFWriter
if not 'nc' in AtomsWriters: AtomsWriters['nc'] = NetCDFWriter
