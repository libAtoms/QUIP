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

from quippy import (ElementName, atoms_reader, AtomsReaders, AtomsWriters, TABLE_STRING_LENGTH,
                    PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL,
                    T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                    T_CHAR, T_LOGICAL, T_INTEGER_A,
                    T_REAL_A, T_COMPLEX_A, T_CHAR_A, T_LOGICAL_A, T_INTEGER_A2, T_REAL_A2)


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
   fmt_map = {'R':'d',
              'I': 'i',
              'S': 'S' + str(TABLE_STRING_LENGTH),
              'L': 'bool'}

   for prop in props:
      ptype, cols = props[prop]
      if cols == 1:
         result.append((prop,fmt_map[ptype]))
      else:
         for c in range(cols):
            result.append((prop+str(c),fmt_map[ptype]))

   return numpy.dtype(result)

def parse_properties(prop_str):
   properties = OrderedDict()
   if prop_str.startswith('pos:R:3'):
      prop_str = 'species:S:1:'+prop_str

   fields = prop_str.split(':')

   for name, ptype, cols in zip(fields[::3], fields[1::3], map(int,fields[2::3])):
      if ptype not in ('R', 'I', 'S', 'L'):
         raise ValueError('Unknown property type: '+ptype)
      properties[name] = (ptype, cols)

   return properties


@atoms_reader('pupyxyz', False)
def PuPyXYZReader(xyz):
   "Read from extended XYZ filename, string or open file."
   from quippy import Table

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

      properties = parse_properties(params['Properties'])
      del params['Properties']

      # Get lattice
      if not 'Lattice' in params:
         raise ValueError('No lattice found in xyz file')

      lattice = params.get_value('Lattice') # make a copy
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

      # Empty dictionary passed for properties to avoid creating pos, species, z.
      at = Atoms(n=n,lattice=lattice,params=params,properties={})

      for p in properties:
         ptype, cols = properties[p]
         if cols == 1:
            value = data[p]
            if value.dtype.kind == 'S':
               value = s2a(value, TABLE_STRING_LENGTH).T
         else:
            value = numpy.vstack([data[p+str(c)] for c in range(cols) ])
         at.add_property(p, value, overwrite=True)

      if not at.has_property('Z')  and not at.has_property('species'):
         raise ValueError('Atoms read from XYZ has neither Z nor species')
      elif at.has_property('Z') and not at.has_property('species'):
         at.add_property('species', ' '*TABLE_STRING_LENGTH)
         at.set_atoms(at.z)
      elif at.has_property('species') and not at.has_property('z'):
         at.add_property('Z', 0)
         at.z[:] = [ElementName.index(sp) for sp in at.species.stripstrings()]

      yield at

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

      def ncols(props, name):
         t, l, (l1,l2) = props.get_type_and_size(name)
         if t in (T_INTEGER_A, T_REAL_A, T_LOGICAL_A, T_CHAR_A):
            return 1 
         elif t in (T_INTEGER_A2, T_REAL_A2):
            return l1
         else:
            raise TypeError('bad property type %d' % t)

      def properties_comment(self, properties=None):
         if properties is None:
            props = self.properties.keys()
         else:
            props = properties

         pmap = { T_INTEGER_A: 'I',
                  T_REAL_A: 'R',
                  T_CHAR_A: 'S',
                  T_LOGICAL_A: 'L',
                  T_INTEGER_A2: 'I',
                  T_REAL_A2: 'R'}

         lattice_str = 'Lattice="' + ' '.join(map(str, numpy.reshape(self.lattice,9,order='F'))) + '"'

         props_str =  ':'.join(map(':'.join,
                                   zip(props,
                                       [pmap[self.properties.get_type_and_size(k)[0]] for k in props],
                                       [str(ncols(self.properties, k)) for k in props])))

         return props_str, lattice_str+' Properties='+props_str+' '+str(self.params)


      if properties is None:
         props = at.properties.keys()
      else:
         props = properties

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

      props_str, comment = properties_comment(at, props)

      subprops = OrderedDict.frompairs([(k, at.properties[k]) for k in props])
      data = numpy.zeros(at.n, _props_dtype(parse_properties(props_str)))

      for prop in props:
         value = getattr(at, prop)
         if at.properties.get_type_and_size(prop)[0] in (T_INTEGER_A2, T_REAL_A2):
            for c in range(value.shape[0]):
               data[prop+str(c)] = value[c+1,:]
         else:
            if value.dtype.kind == 'S':
               value = a2s(value)
            data[prop] = value
      
      format = ''.join([_getfmt(data.dtype.fields[name][0]) for name in data.dtype.names])+'\n'

      self.xyz.write('%d\n' % at.n)
      self.xyz.write('%s\n' % comment)
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

      def write(self, at, frame=None, properties=None, prefix=None):
         self.dest.write(at, frame=frame, properties=properties, prefix=prefix)

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

      at = Atoms(n=netcdf_dimlen(source, 'atom'), lattice=lattice, properties={})

      for name, var in source.variables.iteritems():
         name = remap_names.get(name, name)

         if name is None:
            continue

         if 'frame' in var.dimensions:
            if 'atom' in var.dimensions:
               # It's a property
               value = var[frame]
               if value.dtype.kind != 'S': value = value.T
               at.add_property(name, value)
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


      prop_type_ncols_to_dtype_dim = {(T_INTEGER_A2,3):       ('i', ('frame','atom','spatial')),
                                      (T_REAL_A2,3):     ('d', ('frame','atom','spatial')),
                                      (T_LOGICAL_A,3):  ('d', ('frame','atom','spatial')),
                                      (T_INTEGER_A,1):      ('i', ('frame', 'atom')),
                                      (T_REAL_A, 1):    ('d', ('frame', 'atom')),
                                      (T_LOGICAL_A, 1): ('i', ('frame', 'atom')),
                                      (T_CHAR_A,1):      ('S1', ('frame','atom','label'))}


      param_type_to_dtype_dim = {T_INTEGER:   ('i', ('frame',)),
                                 T_REAL:      ('d', ('frame',)),
                                 T_CHAR:      ('S1', ('frame', 'string')),
                                 T_LOGICAL:   ('i', ('frame',)),
                                 T_INTEGER_A: ('i', ('frame', 'spatial')),
                                 T_REAL_A:    ('d', ('frame', 'spatial')),
                                 T_LOGICAL_A: ('i', ('frame', 'spatial'))}

      # Temporary hack to retain compatibility with old NetCDF property codes
      convert_ptype = {
         T_INTEGER_A2: PROPERTY_INT,
         T_REAL_A2: PROPERTY_REAL,
         T_LOGICAL_A: PROPERTY_LOGICAL,
         T_INTEGER_A: PROPERTY_INT,
         T_REAL_A:   PROPERTY_REAL,
         T_LOGICAL_A: PROPERTY_LOGICAL,
         T_CHAR_A : PROPERTY_STR
         }

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

      for origname in at.properties:
         name = remap_names_rev.get(origname, origname)
         ptype, s1, (s2, s3) = at.properties.get_type_and_size(origname)
         if ptype in (T_INTEGER_A, T_REAL_A, T_LOGICAL_A, T_CHAR_A):
            ncols = 1
         elif ptype in (T_INTEGER_A2, T_REAL_A2):
            ncols = s2
         else:
            raise TypeError('bad property type %d' % ptype)
         dtype, dims = prop_type_ncols_to_dtype_dim[(ptype, ncols)]

         if not name in self.dest.variables:
            self.dest.createVariable(name, dtype, dims)
            self.dest.variables[name].type = convert_ptype[ptype]

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
