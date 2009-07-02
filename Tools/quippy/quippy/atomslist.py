import sys, os, itertools
from quippy import Atoms, AtomsReaders, AtomsWriters


class AtomsList(object):
   """Class for dealing with sequences of Atoms objects"""

   def __init__(self, source, format=None, lazy=True, *args, **kwargs):

      if format is None:
         if isinstance(source, str):
            if source in AtomsReaders:
               format = source
            else:
               base, ext = os.path.splitext(source)
               format = ext[1:]
         else:
            format = source.__class__

      if format in AtomsReaders:
         self._source = iter(AtomsReaders[format](source, *args, **kwargs))
      else:
         if not hasattr(source, '__iter__'):
            self._source = [source]
         else:
            self._source = source
         
      self._randomaccess = hasattr(self._source, '__getitem__') and hasattr(self._source, '__len__')
      if lazy:
         if self._randomaccess:
            self._list = [None for a in range(len(self._source)) ]
         else:
            self._list = []
      else:
         self._list = list(self._source)

   def __repr__(self):
      return '%s(source=%r)' % (self.__class__.__name__, self._source)

   def __str__(self):
      return '%s(%r)' % (self.__class__.__name__, self._list)

   def __len__(self):
      return len(self._list)

   def __delitem__(self, frame):
      a = self._list[frame]
      self._list[frame] = None
      del a

   def iteratoms(self):
      for frame in range(len(self._list)):
         yield self[frame]
      if not self._randomaccess:
         for at in self._source:
            self._list.append(at)
            yield at
      self._randomaccess = True

   def __iter__(self):
      return self.iteratoms()

   def __reversed__(self):
      if not self._randomaccess:
         raise ValueError('Cannot reverse AtomsList with generator source. Call loadall() then try again')
      for frame in reversed(range(len(self._list))):
         yield self[frame]

   def __getslice__(self, first, last):
      return self.__getitem__(slice(first,last,None))

   def __getitem__(self, frame):
      if isinstance(frame, int):
         if frame < 0: frame = frame + len(self._list)
         if not self._randomaccess and frame >= len(self._list):
            self._list.extend(list(itertools.islice(self._source, frame-len(self._list)+1)))
         if self._list[frame] is None:
            self._list[frame] = self._source[frame]
         return self._list[frame]

      elif isinstance(frame, slice):
         if not self._randomaccess and frame.stop >= len(self._list):
            self._list.extend(list(itertools.islice(self._source, frame.stop != sys.maxint and frame.stop-len(self._list)+1 or None)))
         allatoms = range(len(self._list))
         subatoms = [ allatoms[i] for i in range(*frame.indices(len(allatoms))) ]
         res = []
         for f in subatoms:
            if self._list[f] is None:
               self._list[f] = self._source[f]
            res.append(self._list[f])
         return AtomsList(res, lazy=False)
      else:
         raise TypeError('frame should be either an integer or a slice')

   def __setslice__(self, i, j, value):
      self._list[i:j] = value

   def __setitem__(self, i, value):
      self._list[i] = value

   def loadall(self, progress=False, progress_width=80, update_interval=None, show_value=True):
      if progress:
         if self._randomaccess:
            from progbar import ProgressBar
            pb = ProgressBar(0,len(self),progress_width,showValue=show_value)
            update_interval = update_interval or max(1, len(self)/progress_width)
         else:
            update_interval = update_interval or 100
      
      for i, a in enumerate(self.iteratoms()):
         if progress and i % update_interval == 0:
            if self._randomaccess:
               pb(i)
            else:
               print i

   def show(self, property=None, frame=0):
      try:
         import atomeye
         atomeye.show(self, property, frame)
      except ImportError:
         raise RuntimeError('AtomEye not available')

   def write(self, dest, format=None, progress=False, progress_width=80, update_interval=None,
             show_value=True, *args, **kwargs):
      if format is None:
         if isinstance(dest, str):
            if dest in AtomsWriters:
               format = dest
            else:
               base, ext = os.path.splitext(dest)
               format = ext[1:]
         else:
            format = dest.__class__
      
      if progress:
         if self._randomaccess:
            from progbar import ProgressBar
            pb = ProgressBar(0,len(self),progress_width,showValue=show_value)
            update_interval = update_interval or max(1, len(self)/progress_width)
         else:
            update_interval = update_interval or 100

      if format in AtomsWriters:
         dest = iter(AtomsWriters[format](dest, *args, **kwargs))
         dest.next()

      res = []
      for i, a in enumerate(self.iteratoms()):
         res.append(a.write(dest))
         if progress and i % update_interval == 0:
            if self._randomaccess:
               pb(i)
            else:
               print i
      return res


try:
   from quippy import CInOutput, INPUT, OUTPUT, INOUT

   class CInOutputReader(object):
      """Generator to read atoms sequentially from a CInOutput"""

      def __init__(self, source, frame=None, start=0, stop=None, step=1):
         if isinstance(source, str):
            self.source = CInOutput(source, action=INPUT, append=False)
         else:
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

         if self.source.got_index:
            self.__len__ = lambda self: len(self.source)
            self.__getitem__ = lambda self, idx: self.source[idx]

      def __iter__(self):
         if not self.source.got_index:
            while True:
               try:
                  yield self.source.read()
               except RuntimeError:
                  break
         else:
            for frame in range(self.start,self.stop,self.step):
               yield self.source[frame]


   def CInOutputWriter(dest, append=False):
      """Generator to write atoms sequentially to a CInOutput"""

      if isinstance(dest, str):
         dest = CInOutput(dest, action=OUTPUT, append=append)

      at = yield None
      while True:
         dest.write(at)
         at = yield None

   AtomsReaders['xyz'] = AtomsReaders['nc'] = AtomsReaders[CInOutput] = CInOutputReader
   AtomsWriters['xyz'] = AtomsWriters['nc'] = AtomsWriters[CInOutput] = CInOutputWriter
   AtomsReaders['stdin'] = CInOutputReader
   AtomsWriters['stdout'] = CInOutputWriter

   
except ImportError:
   print 'Warning: CInOutput not found - Atoms read/write disabled'


try:
   from netCDF4 import Dataset

   def netcdf4_read(source, frame=0):

      if isinstance(source, str):
         source = Dataset(source)
      
      from quippy import Atoms, make_lattice
      from math import pi
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


      cl = source.variables['cell_lengths'][frame]
      ca = source.variables['cell_angles'][frame]
      lattice = make_lattice(cl[0],cl[1],cl[2],ca[0]*DEG_TO_RAD,ca[1]*DEG_TO_RAD,ca[2]*DEG_TO_RAD)

      at = Atoms(n=len(source.dimensions['atom']), lattice=lattice)

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

   AtomsReaders[Dataset] = netcdf4_read
   if not 'nc' in AtomsReaders:
      AtomsReaders['nc'] = netcdf4_read
   

   class NetCDFAtomsList(Dataset, AtomsList):
      def __init__(self, source):
         AtomsList.__init__(self, source)
         Dataset.__init__(self, source)

      def __getattr__(self, name):
         return Dataset.__getattr__(self,name)

      def __setattr__(self, name, value):
         if name.startswith('_'):
            self.__dict__[name] = value
         else:
            Dataset.__setattr__(self, name, value)

except ImportError:
   print 'Warning: netCDF4 not found. NetCDFAtomsReader disabled.'
      
