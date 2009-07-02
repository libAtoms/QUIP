"""Various subclasses to add functionality to automatically generated classes."""

import numpy, hashlib, os
from farray import *

from quippy import FortranAtoms, FortranDictionary, FortranTable, FortranDynamicalSystem, FortranCInOutput
from quippy import AtomsReaders, AtomsWriters

class Atoms(FortranAtoms):
   """
   Thin pythonic wrapper over auto-generated FortranAtoms class
   
   An atoms object contains atomic numbers, all dynamical variables
   and connectivity information for all the atoms in the simulation cell.
   It is initialised like this:
   >         call initialise(MyAtoms,N,lattice)
   where 'N' is the number of atoms to allocate space for and 'lattice' is a $3\times3$
   matrix of lattice vectors given as column vectors.
   
   Atoms also contains a Connection object, which stores distance information about
   the atom neghbours after 'calc_connect' has been called. Rather than using a minimum
   image convention, all neighbours are stored up to a radius of 'cutoff', including images
   """

   def __new__(cls, source=None, n=0, lattice=fidentity(3), fpointer=None, finalise=True,
                data=None, properties=None, params=None, *readargs, **readkwargs):
      """
      Initialise an Atoms object.

      The arugments, n, lattice, fpointer, finalise, data, properties and params are passed on
      to the FortranAtoms constructor.

      If source is given, Atoms.read(source, frame=frame) is called.
      """

      if source is None:
         self = object.__new__(cls)
      else:
         self = cls.read(source, *readargs, **readkwargs)

      return self

   def __init__(self, source=None, n=0, lattice=fidentity(3), fpointer=None, finalise=True,
                data=None, properties=None, params=None, *readargs, **readkwargs):
      if source is None:
         FortranAtoms.__init__(self, n=n, lattice=lattice, fpointer=fpointer, finalise=finalise,
                               data=data, properties=properties, params=params)

   @classmethod
   def read(cls, source, format=None, *args, **kwargs):
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
         source = iter(AtomsReaders[format](source, *args, **kwargs))

      at = source.next()
      if not isinstance(at, cls):
         raise ValueError('Object %r read from source %r is not Atoms instance' % (at, source))
      return at
         

   def write(self, dest, format=None, *args, **kwargs):
      if format is None:
         if isinstance(dest, str):
            if dest in AtomsWriters:
               format = dest
            else:
               base, ext = os.path.splitext(dest)
               format = ext[1:]
         else:
            format = dest.__class__

      if format in AtomsWriters:
         dest = iter(AtomsWriters[format](dest, *args, **kwargs))
         dest.next()

      return dest.send(self)

   def show(self, property=None):
      """Show this Atoms object in AtomEye."""
      try:
         import atomeye
         atomeye.show(self,property=property)
      except ImportError:
         raise RuntimeError('AtomEye not available')

   def select(self, mask=None, list=None):
      """Select a subset of the atoms in an atoms object

      Use either a logical mask array or a list of atom indices to include.

      small_at = at.select([mask, list])
      
      """
      if mask is not None:
         out = Atoms(n=mask.count(),lattice=self.lattice)
         FortranAtoms.select(out, self, mask=mask)
      elif list is not None:
         out = Atoms(n=list.size(), lattice=self.lattice)
         FortranAtoms.select(out, self, list=list)
      else:
         raise ValueError('Either mask or list must be present.')
      return out

   def set_species(self, i, species):
      self.species[i] = [x for x in species] + [' ' for a in range(self.species.shape[0]-len(species))]
    
   def _update_hook(self):
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
               #setattr(self, prop, self.data.real[col_start-1,0:self.n])
            else:
               setattr(self, prop, FortranArray(self.data.real[col_start:col_end,1:self.n],doc,transpose_on_print=True))
               #setattr(self, prop, self.data.real[col_start-1:col_end,0:self.n])
         elif ptype == PROPERTY_INT:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.int[col_start,1:self.n],doc))
               #setattr(self, prop, self.data.int[col_start-1,0:self.n])
            else:
               setattr(self, prop, FortranArray(self.data.int[col_start:col_end,1:self.n],doc,transpose_on_print=True))
               #setattr(self, prop, self.data.int[col_start-1:col_end,0:self.n])
         elif ptype == PROPERTY_STR:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.str[:,col_start,1:self.n],doc))
               #setattr(self, prop, self.data.str[:,col_start-1,0:self.n])
            else:
               setattr(self, prop, FortranArray(self.data.str[:,col_start:col_end,1:self.n],doc,transpose_on_print=True))
               #setattr(self, prop, self.data.str[:,col_start-1:col_end,0:self.n])
         elif ptype == PROPERTY_LOGICAL:
            if col_end == col_start:
               setattr(self, prop, FortranArray(self.data.logical[col_start,1:self.n],doc))
               #setattr(self, prop, self.data.logical[col_start-1,0:self.n])
            else:
               setattr(self, prop, FortranArray(self.data.logical[col_start:col_end,1:self.n],doc,transpose_on_print=True))
               #setattr(self, prop, self.data.logical[col_start-1:col_end,0:self.n])
         else:
            raise ValueError('Bad property type :'+str(self.properties[prop]))

   def copy(self):
      other = Atoms(n=self.n, lattice=self.lattice, data=self.data, 
                    properties=self.properties, params=self.params)
      return other

   def calc_bonds(self):
      if hasattr(self, '_bonds'): return
      
      if not self.connect.initialised:
         self.calc_connect()
         
      self._bonds = []
      from matplotlib.lines import Line2D
      for i in frange(self.n):
         for n in frange(self.n_neighbours(i)):
            if self.is_nearest_neighbour(i, n):
               j = self.neighbour(i, n)
               self._bonds.append(Line2D((self.pos[1,i],self.pos[1,j]), (self.pos[2,i],self.pos[2,j]),zorder=-100))

   def plot(self, bonds=False, **kwargs):
      import pylab

      ax = pylab.gca()
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

   def __len__(self):
      return self.n

   def __getattr__(self, name):
      if isinstance(self.params, FortranDictionary) and name in self.params:
         return self.params[name]
      raise AttributeError('Unknown Atoms parameter %s' % name)

   def __getitem__(self, i):
      if i < 1 or i > self.n:
         raise ValueError('Atoms index should be in range 1..self.n(%d)' % self.n)
      res = {}
      for k in self.properties.keys():
         v = getattr(self,k.lower())[...,i]
         if isinstance(v,FortranArray) and v.dtype.kind == 'S':
             v = ''.join(v).strip()
         res[k] = v
      return res

   def __eq__(self, other):
      tol = 1e-10
      
      if self.n != other.n: return False
      if sorted(self.properties.keys()) != sorted(other.properties.keys()): return False
      for key in self.properties:
         v1, v2 = self.properties[key], other.properties[key]
         if v1[1] != v2[1]: return False # type
         if v1[3]-v1[2] != v2[3]-v2[2]: return False # ncols
         
      if self.params != other.params: return False
      if abs(self.lattice - other.lattice).max() > tol: return False
      if self.data != other.data: return False
      return True

   def __ne__(self, other):
      return not self.__eq__(other)
      

from dictmixin import DictMixin
class Dictionary(DictMixin, FortranDictionary):

   def __init__(self, *args, **kwargs):
      FortranDictionary.__init__(self, *args, **kwargs)
      if 'keys' in self._arrays:
         self._arrays['_keys'] = self._arrays['keys']
         del self._arrays['keys']
         del self.keys
         self._update()

   def keys(self):
      new_hash = hashlib.md5(self._keys).hexdigest()
      if new_hash != self.__dict__.get('_keys_hash', None):
         self._cache_misses = self.__dict__.get('_cache_misses',0) + 1
         self._keys_cache = [''.join(self._keys[:,i]).strip() for i in frange(self.n)]
      else:
         self._cache_hits = self.__dict__.get('_cache_hits',0) + 1
      self._keys_hash = new_hash
      return self._keys_cache

   def slowkeys(self):
      return [''.join(self._keys[:,i]).strip() for i in frange(self.n)]

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
         v,p = self._get_value_i(k)
      elif t == T_REAL:
         v,p = self._get_value_r(k)
      elif t == T_COMPLEX:
         v,p = self._get_value_c(k)
      elif t == T_CHAR:
         v,p = self._get_value_s(k)
         v = v.strip()
      elif t == T_LOGICAL:
         v,p = self._get_value_l(k)
      elif t == T_INTEGER_A:
         v,p = self._get_value_i_a(k,s)
         v = farray(v)
      elif t == T_REAL_A:
         v,p = self._get_value_r_a(k,s)
         v = farray(v)
      elif t == T_COMPLEX_A:
         v,p = self._get_value_c_a(k,s)
         v = farray(v)
      elif t == T_CHAR_A:
         a,p = self._get_value_s_a(k,s)
         v = [''.join(line).strip() for line in a]
      elif t == T_LOGICAL_A:
         v,p = self._get_value_l_a(k,s)
         v = farray(v)
      else:
         raise ValueError('Unsupported dictionary entry type %d' % t)

      return v

   def __setitem__(self, k, v):
      self.set_value(k, v)

   def __repr__(self):
      return 'Dictionary(%s)' % DictMixin.__repr__(self)

   def __eq__(self, other):
      if self.keys() != other.keys(): return False

      for key in self:
         v1, v2 = self[key], other[key]
         if isinstance(v1, FortranArray):
            if (v1 != v2).any(): return False
         else:
            if v1 != v2: return False
      return True
         

   def __ne__(self, other):
      return not self.__eq__(other)


class Table(FortranTable):

   def copy(self):
      t = Table(self.intsize, self.realsize, self.strsize, self.logicalsize, self.n)
      t.n = self.n
      t._update()
      t.int[...] = self.int[...]
      t.real[...] = self.real[...]
      t.str[...] = self.str[...]
      t.logical[...] = self.logical[...]
      return t

   def __eq__(self, other):
      return self.equal(other)

   def __ne__(self, other):
      return not self.equal(other)

   def equal(self, other):
      tol = 1e-10
      
      for t1, t2 in zip((self.n,  self.intsize,  self.realsize,  self.strsize,  self.logicalsize),
                        (other.n, other.intsize, other.realsize, other.strsize, other.logicalsize)):
         if t1 != t2: return False

      for n, a1, a2, in zip((self.intsize, self.realsize, self.strsize, self.logicalsize),
                            (self.int,  self.real,  self.str,  self.logical),
                            (other.int, other.real, other.str, other.logical)):

         if n == 0: continue
         try:
            if abs(a1 - a2).max() > tol: return False
         except TypeError:
            if (a1 != a2).any(): return False

      return True


   def _get_array_shape(self, name):
       if name in ('int','real','logical'):
           return (slice(None),slice(1,self.n))
       elif name == 'str':
           return (slice(None),slice(None),slice(1,self.n))
       else:
          return None


class Trajectory(object):
   def __init__(self, ds, pot, dt, n_steps, save_interval, connect_interval, args_str):
      self.ds = ds
      self.pot = pot
      self.dt = dt
      self.n_steps = n_steps
      self.save_interval = save_interval
      self.connect_interval = connect_interval
      self.args_str = args_str

      self.ds.atoms.add_property('force', 0.0, n_cols=3)
      self.f = fzeros((3,ds.atoms.n))
      self.e = farray(0.0)
      self.pot.calc(ds.atoms, f=self.f, e=self.e, args_str=self.args_str)
      self.ds.atoms.force[:] = self.f[:]
      self.ds.atoms.params['energy'] = self.e
      
   def __iter__(self):
      for n in range(self.n_steps):
         self.ds.advance_verlet1(self.dt, self.f)
         self.pot.calc(self.ds.atoms, e=self.e, f=self.f, args_str=self.args_str)
         self.ds.atoms.force[:] = self.f[:]
         self.ds.atoms.params['energy'] = self.e
         self.ds.advance_verlet2(self.dt, self.f)
         self.ds.print_status(epot=self.e)
         if n % self.connect_interval == 0:
            self.ds.atoms.calc_connect()
         if n % self.save_interval == 0:
            yield self.ds.atoms.copy()
      raise StopIteration


class DynamicalSystem(FortranDynamicalSystem):
   
   def run(self, pot, dt=1.0, n_steps=10, save_interval=1, connect_interval=10, args_str=""):
      return Trajectory(self, pot, dt, n_steps, save_interval, connect_interval, args_str)


from quippy import INPUT, OUTPUT, INOUT


class CInOutput(FortranCInOutput):

   def __init__(self, filename, action=INPUT, append=False):
      FortranCInOutput.__init__(self, filename, action, append)

   def __len__(self):
      return int(self.n_frame)

   def __getitem__(self, index):
      return self.read(frame=index)

   def __setitem__(self, index, value):
      self.write(value, frame=index)


