class AtomsList(object):
   def __init__(self, seq):
      if not hasattr(seq, '__iter__'):
         seq = [seq]
      self._list = list(seq)
   
   def __repr__(self):
      return '%s(%r)' % (self.__class__.__name__, self._list)

   def show(self, property=None, frame=0):
      try:
         import atomeye
         atomeye.show(self, property, frame)
      except ImportError:
         raise RuntimeError('AtomEye not available')

   def __getitem__(self, idx):
      return self._list[idx]

   def __iter__(self):
      return self._list.__iter__()

   def __len__(self):
      return len(self._list)

   def write(self, dest):
      for a in self:
         dest.write(self)

       
class GenericFrameReader(AtomsList):
   """Read-only access to an XYZ or NetCDF trajectory. The file is opened
   and then read lazily as frames are asked for. Supports list-like interface:

   fr = FrameReader('foo.xyz')
   at1 = fr[0]      # First frame
   at2 = fr[-1]     # Last frame
   ats = fr[0:10:3] # Every third frame between 0 and 10
   ats = [ a for a in fr if a.n == 100 ]  # Only frames with exactly 100 atoms
   """
   def __init__(self, source, start=0, stop=-1, step=None, count=None):

      self._init(source)
      if count is not None:
         self._frames = slice(start,start+count,1)
      else:
         self._frames = slice(start,stop,step)

      self._list = [None for a in range(*(self._frames.indices(self._nframe()+1)))]
      self._iterframes = None


   def __del__(self):
      self._close()

   def __len__(self):
      return len(range(*self._frames.indices(self._nframe()+1)))

   def __getitem__(self, frame):
      start, stop, step = self._frames.indices(self._nframe()+1)

      if isinstance(frame, int):
         if frame < 0: frame = frame + (stop - start)
         if start + frame >= stop:
            raise ValueError("frame %d out of range %d" % (start + frame, stop))
         if self._list[frame] is None:
            self._list[frame] = self._getframe(start+frame)
         return self._list[frame]

      elif isinstance(frame, slice):
         allframes = range(start, stop, step)
         subframes = [ allframes[i] for i in range(*frame.indices(len(allframes))) ]
         res = []
         for f in subframes:
            if self._list[f] is None:
               self._list[f] = self._getframe(f)
            res.append(self._list[f])
         return res
      else:
         raise TypeError('frame should be either an integer or a slice')

   def __delitem__(self, frame):
      a = self._list[frame]
      self._list[frame] = None
      del a

   def __getslice__(self, first, last):
      return self.__getitem__(slice(first,last,None))

   def iterframes(self):
      for frame in range(*self._frames.indices(self._nframe()+1)):
         yield self[frame]
      raise StopIteration

   def __iter__(self):
      if self._iterframes is None:
         self._iterframes = self.iterframes()
      return self._iterframes

   def next(self):
      if self._iterframes is None:
         self._iterframes = self.iterframes()
      return self._iterframes.next()

   def __reversed__(self):
      for frame in reversed(range(*self._frames.indices(self._nframe()+1))):
         yield self[frame]
      raise StopIteration

   def loadall(self, progress=True, progress_width=80, show_value=True):
      if progress:
         from progbar import ProgressBar
         pb = ProgressBar(0,len(self),progress_width,showValue=show_value)
         update_interval = max(1, len(self)/progress_width)
      
      for i in range(len(self)):
         if progress and i % update_interval == 0:
            pb(i)
         self[i]

   def clearall(self):
      for i in range(len(self)):
         del self[i]


try:
   from quippy import CInOutput
   
   class CInOutputFrameReader(GenericFrameReader):
      def _init(self, source):
         self._cio = CInOutput(source)
         self._cio.query()

      def _close(self):
         pass

      def _getframe(self, frame):
         return self._cio.read(frame) 

      def _nframe(self):
         return self._cio.n_frame
   
   FrameReader = CInOutputFrameReader
except ImportError:
   print 'Cannot import CInOutput. CInoutputFrameReader disabled.'

try:
   from netCDF4 import Dataset

   class NetCDFFrameReader(CInOutputFrameReader,Dataset):

      def _init(self, source):
         Dataset.__init__(self, source)
         CInOutputFrameReader._init(self, source)

      def _close(self):
         self.close()
         CInOutputFrameReader._close()

      def __getattr__(self, name):
         return Dataset.__getattr__(self,name)

      def __setattr__(self, name, value):
         if name[0] == '_':
            self.__dict__[name] = value
         else:
            Dataset.__setattr__(self, name, value)
         

##       def _getframe(self, frame):
##          from quippy import Atoms, make_lattice
##          from math import pi
##          from quippy import (PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL,
##                              T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
##                              T_CHAR, T_LOGICAL, T_INTEGER_A,
##                              T_REAL_A, T_COMPLEX_A, T_CHAR_A, T_LOGICAL_A)

##          DEG_TO_RAD = pi/180.0

##          remap_names = {'coordinates': 'pos',
##                         'velocities': 'velo',
##                         'cell_lengths': None,
##                         'cell_angles': None}
         
##          prop_type_to_value = {PROPERTY_INT: 0,
##                                PROPERTY_REAL: 0.0,
##                                PROPERTY_STR: "",
##                                PROPERTY_LOGICAL: False}

##          prop_dim_to_ncols = {('frame','atom','spatial'): 3,
##                               ('frame','atom','label'): 1,
##                               ('frame', 'atom'): 1}


##          cl = self.variables['cell_lengths'][frame]
##          ca = self.variables['cell_angles'][frame]
##          lattice = make_lattice(cl[0],cl[1],cl[2],ca[0]*DEG_TO_RAD,ca[1]*DEG_TO_RAD,ca[2]*DEG_TO_RAD)

##          a = Atoms(len(self.dimensions['atom']),lattice)

##          for name, var in self.variables.iteritems():
##             name = remap_names.get(name, name)

##             if name is None:
##                continue

##             if 'frame' in var.dimensions:
##                if 'atom' in var.dimensions:
##                   # It's a property
##                   a.add_property(name, prop_type_to_value[var.type],
##                                  n_cols=prop_dim_to_ncols[var.dimensions])
##                   getattr(a,name.lower())[...] = var[frame].T
##                else:
##                   # It's a param
##                   if var.dimensions == ('frame','string'):
##                      # if it's a single string, join it and strip it
##                      a.params[name] = ''.join(var[frame]).strip()
##                   else:
##                      a.params[name] = var[frame]
                  
##          return a

##       def _nframe(self):
##          return len(self.dimensions['frame'])

except ImportError:
   print 'Cannot import netCDF4. NetCDFFrameReader disabled.'


