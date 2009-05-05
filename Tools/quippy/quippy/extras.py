"""Various mixin classes to add functionality to automatically generated classes."""

import numpy
from farray import *

class AtomsExtras(object):

   @classmethod
   def _init_extra(cls):
      cls._select = cls.select
      del cls.select

   def write(self, target, append=True):
      "Write atoms object to an XYZ or NetCDF file."
      cio = CInOutput(target, OUTPUT, append)
      try:
         cio.write(self)
      finally:
         cio.close()

   @classmethod
   def read(cls, source, frame=0):
      """Read a single frame from an XYZ or NetCDF file, then close the file.
      This is a classmethod so should be called as at = Atoms.read(source)."""
      from atomslist import FrameReader
      at = FrameReader(source, start=frame)[0]
      return at

   def show(self, property=None):
      """Show this atoms object in AtomEye."""
      try:
         import atomeye
         atomeye.show(self,property=property)
      except ImportError:
         raise RuntimeError('AtomEye not available')


   def select(self, mask=None, list=None):
      """Select a subset of the atoms in an atoms object, either using a logical
      mask array or list of atom indices to include.

      small_at = at.select([mask, list])
      
      """
      from quippy import Atoms
      if mask is not None:
         out = Atoms(mask.count(),self.lattice)
         out._select(self, mask=mask)
      elif list is not None:
         out = Atoms(list.size(), self.lattice)
         out._select(self, list=list)
      else:
         raise ValueError('Either mask or list must be present.')
      return out
    
   def update_hook(self):
      # Remove existing pointers
      if hasattr(self, '_props'):
         for prop in self._props:
            try:
               delattr(self, prop)
            except AttributeError:
               pass

      from quippy import PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL

      type_names = {PROPERTY_INT: "integer",
                    PROPERTY_REAL: "real",
                    PROPERTY_STR: "string",
                    PROPERTY_LOGICAL: "logical"}

      self._props = {}
      for prop,(ptype,col_start,col_end) in self.properties.iteritems():
         prop = prop.lower()
         self._props[prop] = (ptype,col_start,col_end)
         doc = "%s : %s property with %d %s" % (prop, type_names[ptype], col_end-col_start+1, 
                                                {True:"column",False:"columns"}[col_start==col_end])
         if ptype == PROPERTY_REAL:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.real[col_start,1:self.n],doc))
            else:
               setattr(self, prop, FortranArray(self.data.real[col_start:col_end,1:self.n],doc))
         elif ptype == PROPERTY_INT:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.int[col_start,1:self.n],doc))
            else:
               setattr(self, prop, FortranArray(self.data.int[col_start:col_end,1:self.n],doc))
         elif ptype == PROPERTY_STR:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.str[:,col_start,1:self.n],doc))
            else:
               setattr(self, prop, FortranArray(self.data.str[:,col_start:col_end,1:self.n],doc))
         elif ptype == PROPERTY_LOGICAL:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.logical[col_start,1:self.n],doc))
            else:
               setattr(self, prop, FortranArray(self.data.logical[col_start:col_end,1:self.n],doc))
         else:
            raise ValueError('Bad property type :'+str(self.properties[prop]))

   def copy(self):
      from quippy import Atoms
      other = Atoms(self.n, self.lattice, self.data, 
                    self.properties, self.params)
      return other

   def calc_bonds(self):
      if hasattr(self, '_bonds'): return
      
      if not self.connect.initialised:
         self.calc_connect()
         
      self._bonds = []
      from matplotlib.lines import Line2D
      import pylab
      fig = pylab.gcf()
      ax = pylab.gca()
      
      for i in frange(self.n):
         for n in frange(self.n_neighbours(i)):
            if self.is_nearest_neighbour(i, n):
               j = self.neighbour(i, n)
               self._bonds.append(Line2D((self.pos[1,i],self.pos[1,j]), (self.pos[2,i],self.pos[2,j]),zorder=-100))

   def plot(self, bonds=False, **kwargs):
      import pylab

      ax = pylab.gca()

      print kwargs
      ax.scatter(self.pos[1,:],self.pos[2,:], **kwargs)

      if bonds:
         self.calc_bonds()
         for line in self._bonds:
            ax.add_line(line)

      pylab.gcf().show()

   def plot_arrows(self, value):

      import pylab
      pylab.quiver(numpy.array(self.pos[1,:]),numpy.array(self.pos[2,:]),
                   numpy.array(value[1,:]),numpy.array(value[2,:]))



from dictmixin import DictMixin
class DictionaryExtras(DictMixin):

   @classmethod
   def _init_extra(cls):
      cls._arrays['_keys'] = cls._arrays['keys']
      del cls._arrays['keys']
      del cls.keys

   def keys(self):
      return[''.join(self._keys[:,i]).strip() for i in frange(self.n)]

   def __getitem__(self, k):
      i = self.lookup_entry_i(k)
      if i == -1: 
         raise KeyError('Key "%s" not found ' % k)

      t, s = self.get_type_and_size(k)

      from quippy import (T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                          T_CHAR, T_LOGICAL, T_INTEGER_A,
                          T_REAL_A, T_COMPLEX_A, T_CHAR_A,
                          T_LOGICAL_A)

      if t == T_NONE:
         raise ValueError('Key %s has no associated value' % k)
      elif t == T_INTEGER:
         v,p = self.get_value_i(k)
      elif t == T_REAL:
         v,p = self.get_value_r(k)
      elif t == T_COMPLEX:
         v,p = self.get_value_c(k)
      elif t == T_CHAR:
         v,p = self.get_value_s(k)
         v = v.strip()
      elif t == T_LOGICAL:
         v,p = self.get_value_l(k)
      elif t == T_INTEGER_A:
         v,p = self.get_value_i_a(k,s)
      elif t == T_REAL_A:
         v,p = self.get_value_r_a(k,s)
      elif t == T_COMPLEX_A:
         v,p = self.get_value_c_a(k,s)
      elif t == T_CHAR_A:
         a,p = self.get_value_s_a(k,s)
         v = [''.join(line).strip() for line in a]
      elif t == T_LOGICAL_A:
         v,p = self_get_value_l_a(k,s)
      else:
         raise ValueError('Unsupported dictionary entry type %d' % t)

      return v

   def __setitem__(self, k, v):
      if type(v) == type(0):
         self.set_value_i(k,v)
      elif type(v) == type(0.0):
         self.set_value_r(k,v)
      elif type(v) == type(0j):
         self.set_value_c(k,v)
      elif type(v) == type(''):
         self.set_value_s(k,v)
      elif type(v) == type(True):
         self.set_value_l(k,v)
      elif hasattr(v, '__iter__'):
         v0 = v[0]
         if type(v0) == type(0):
            self.set_value_i_a(k,v)
         elif type(v0) == type(0.0):
            self.set_value_r_a(k,v)
         elif type(v0) == type(0j):
            self.set_value_c_a(k,v)
         elif type(v0) == type(''):
            self.set_value_s_a(k,numpy.array(v))
         elif type(v0) == type(True):
            self.set_value_l_a(k,v)
      else:
         raise ValueError('Unsupported dictionary entry type %s' % v)

   def __repr__(self):
      return 'Dictionary(%s)' % DictMixin.__repr__(self)


class TableExtras(object):

   @classmethod
   def _init_extra(cls):
      pass
   
   def copy(self):
      from quippy import Table
      t = Table(self.intsize, self.realsize, self.strsize, self.logicalsize, self.n)
      t.n = self.n
      t._update()
      t.int[...] = self.int[...]
      t.real[...] = self.real[...]
      t.str[...] = self.str[...]
      t.logical[...] = self.logical[...]
      return t

   def equal(self, other):
      for t1, t2 in zip((self.n,  self.intsize,  self.realsize,  self.strsize,  self.logicalsize),
                        (other.n, other.intsize, other.realsize, other.strsize, other.logicalsize)):
         if t1 != t2: return False

      for a1, a2, in zip((self.int,  self.real,  self.str,  self.logical),
                         (other.int, other.real, other.str, other.logical)):
         if (a1 != a2).any(): return False

      return True


   def map_array_shape(self, name, shape):
       if name in ('int','real','logical'):
           return (slice(None),slice(1,self.n))
       elif name == 'str':
           return (slice(None),slice(None),slice(1,self.n))


class DynamicalSystemExtras(object):

   @classmethod
   def _init_extra(cls):
      pass

   def run(self, pot, dt=1.0, n_steps=10, save_interval=1, connect_interval=10):
      from atomslist import AtomsList, Trajectory
      traj = Trajectory(self, pot, dt, n_steps, save_interval, connect_interval)
      return AtomsList(traj)



extras_map = {'Atoms': AtomsExtras,
              'Dictionary': DictionaryExtras,
              'Table': TableExtras,
              'DynamicalSystem': DynamicalSystemExtras}

