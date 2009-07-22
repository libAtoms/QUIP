import sys, os, itertools
from quippy import Atoms, AtomsReaders, AtomsWriters

def AtomsReader(source, format=None, *args, **kwargs):
   """Generator to read sucessive frames from source"""

   if format is None:
      if isinstance(source, str):
         if source in AtomsReaders:
            format = source
         else:
            base, ext = os.path.splitext(source)
            format = ext[1:]
      else:
         format = source.__class__

   opened = False
   if format in AtomsReaders:
      source = iter(AtomsReaders[format](source, *args, **kwargs))
      opened = True

   for at in source:
      yield at

   if opened and hasattr(source, 'close'):
      source.close()
   

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
         self._source = AtomsReaders[format](source, *args, **kwargs)
         self._itersource = iter(self._source)
      else:
         if not hasattr(source, '__iter__'):
            self._source = [source]
            self._itersource = iter(self._source)
         else:
            self._source = source
            self._itersource = iter(self._source)
         
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
         for at in self._itersource:
            self._list.append(at)
            yield at
      if hasattr(self._source, 'close'):
         self._source.close()
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
            self._list.extend(list(itertools.islice(self._itersource, frame-len(self._list)+1)))
         if self._list[frame] is None:
            self._list[frame] = self._source[frame]
         return self._list[frame]

      elif isinstance(frame, slice):
         if not self._randomaccess and frame.stop >= len(self._list):
            self._list.extend(list(itertools.islice(self._itersource, frame.stop != sys.maxint and frame.stop-len(self._list)+1 or None)))
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

      opened = False
      if format in AtomsWriters:
         dest = AtomsWriters[format](dest, *args, **kwargs)
         opened = True
         
      res = []
      for i, a in enumerate(self.iteratoms()):
         res.append(dest.write(a))
         if progress and i % update_interval == 0:
            if self._randomaccess:
               pb(i)
            else:
               print i

      if opened:
         dest.close()
      return res

# Decorator to add a new reader

def atoms_reader(source):
   def decorate(func):
      from quippy import AtomsReaders
      if not source in AtomsReaders:
         AtomsReaders[source] = func
      return func
   return decorate


