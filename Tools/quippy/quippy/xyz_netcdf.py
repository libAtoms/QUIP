from quippy import atoms_reader, AtomsReaders, AtomsWriters
from atomslist import *
from farray import *
import logging
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

   lattice[2,1] = b * cos_gamma
   lattice[2,2] = b * sin_gamma

   lattice[3,1] = c * cos_beta
   lattice[3,2] = c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma
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
           acos(dot(a_2,a_3))/a_2.norm()*a_3.norm(),
           acos(dot(a_3,a_1))/a_3.norm()*a_1.norm(),
           acos(dot(a_1,a_2))/a_1.norm()*a_2.norm())

try:
   from quippy import CInOutput, INPUT, OUTPUT, INOUT

   class CInOutputReader(object):
      """Generator to read atoms sequentially from a CInOutput"""

      def __init__(self, source, frame=None, start=0, stop=None, step=1):
         if isinstance(source, str):
            self.opened = True
            self.source = CInOutput(source, action=INPUT, append=False)
         else:
            self.opened = False
            self.source = source

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

   AtomsReaders['xyz'] = AtomsReaders['nc'] = AtomsReaders[CInOutput] = CInOutputReader

   @atoms_reader('stdin')
   def CInOutputStdinReader(source):
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

      def write(self, at):
         self.dest.write(at)

      def close(self):
         self.dest.close()

   AtomsWriters['xyz'] = AtomsWriters['nc'] = AtomsWriters[CInOutput] = AtomsWriters['stdout'] = CInOutputWriter

   
except ImportError:
   logging.warning('CInOutput not found - falling back on (slower) pure python I/O')


try:
   from netCDF4 import Dataset
   netcdf_file = Dataset
except ImportError:
   from pupynere import netcdf_file
   logging.warning('netCDF4 not found. falling back on (slower) pupynere.')

def netcdf_dimlen(obj, name):
   """Return length of dimension 'name'. Works for both netCD4 and pupynere."""
   n = obj.dimensions[name]
   try:
      return len(n)
   except TypeError:
      return n


@atoms_reader(netcdf_file)
@atoms_reader('nc')
def NetCDFReader(source, frame=None, start=0, stop=None, step=1):

   opened = False
   if isinstance(source, str):
      opened = True
      source = netcdf_file(source)

   from quippy import Atoms
   from quippy import (PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL,
                       T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                       T_CHAR, T_LOGICAL, T_INTEGER_A,
                       T_REAL_A, T_COMPLEX_A, T_CHAR_A, T_LOGICAL_A)

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
                  at.params[name] = var[frame]

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
      from quippy import (PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL,
                          T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                          T_CHAR, T_LOGICAL, T_INTEGER_A,
                          T_REAL_A, T_COMPLEX_A, T_CHAR_A, T_LOGICAL_A)

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

class NetCDFAtomsList(netcdf_file, AtomsList):
   def __init__(self, source):
      AtomsList.__init__(self, source)
      netcdf_file.__init__(self, source)

   def __getattr__(self, name):
      return netcdf_file.__getattr__(self,name)

   def __setattr__(self, name, value):
      if name.startswith('_'):
         self.__dict__[name] = value
      else:
         netcdf_file.__setattr__(self, name, value)
      
