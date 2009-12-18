"""Various subclasses to add functionality to automatically generated classes."""

import sys

major, minor = sys.version_info[0:2]

assert (major, minor) >= (2, 4)

if (major, minor) < (2, 5):
   import md5
   got_hashlib = False
else:
   import hashlib
   got_hashlib = True

import numpy, os, weakref
from farray import *

from quippy import FortranAtoms, FortranDictionary, FortranTable, FortranDynamicalSystem, FortranCInOutput
from quippy import AtomsReaders, AtomsWriters

class NeighbourInfo(object):
   __slots__ = ('j','distance','diff','cosines','shift')

   def __init__(self, j, distance, diff, cosines, shift):
      self.j = j
      self.distance = float(distance)
      self.diff = diff.copy()
      self.cosines = cosines.copy()
      self.shift = shift.copy()

   def __int__(self):
      return self.j

   def __repr__(self):
      return 'NeighbourInfo(j=%d, distance=%f, diff=%s, cosines=%s, shift=%s)' % (self.j, self.distance, self.diff, self.cosines, self.shift)

class Neighbours(object):

   def __init__(self, at):
      self.atref = weakref.ref(at)

   def __iter__(self):
      return self.iterneighbours()

   def iterneighbours(self):
      at = self.atref()

      for i in frange(at.n):
         yield self[i]

   def __getitem__(self, i):

      at = self.atref()
      if not at.connect.initialised:
         at.calc_connect()
      
      distance = farray(0.0)
      diff = fzeros(3)
      cosines = fzeros(3)
      shift = fzeros(3,dtype=int)

      res = []
      for n in frange(at.n_neighbours(i)):
         j = at.neighbour(i, n, distance, diff, cosines, shift)
         res.append(NeighbourInfo(j,distance,diff,cosines,shift))

      return farray(res) # to give 1-based indexing

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

      The arguments, n, lattice, fpointer, finalise, data, properties and params are passed on
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
         if params is not None and not isinstance(params, Dictionary):
            params = Dictionary(params)
         if properties is not None and not isinstance(properties, Dictionary):
            properties = Dictionary(properties)
         FortranAtoms.__init__(self, n=n, lattice=lattice, fpointer=fpointer, finalise=finalise,
                               data=data, properties=properties, params=params)
         self.neighbours = Neighbours(self)

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

      opened = False
      if format in AtomsReaders:
         source = AtomsReaders[format](source, *args, **kwargs)
         opened = True

      if not hasattr(source, '__iter__'):
         raise ValueError('Cannot read from source %r - not an iterator' % source)
      at = iter(source).next()
      if not isinstance(at, cls):
         raise ValueError('Object %r read from source %r is not Atoms instance' % (at, source))
      if opened and hasattr(source, 'close'):
         source.close()
      return at
         

   def write(self, dest, format=None, *args, **kwargs):
      opened = False
      if format is None:
         if isinstance(dest, str):
            opened = True
            if dest in AtomsWriters:
               format = dest
            else:
               base, ext = os.path.splitext(dest)
               format = ext[1:]
         else:
            format = dest.__class__

      if format in AtomsWriters:
         dest = AtomsWriters[format](dest, *args, **kwargs)

      res = dest.write(self)
      if opened and hasattr(dest, 'close'):
         dest.close()
      return res


   def show(self, property=None, highlight=None, arrows=None, *arrowargs, **arrowkwargs):
      """Show this Atoms object in AtomEye."""
      try:
         import atomeye
         atomeye.show(self,property=property, highlight=highlight, arrows=arrows, *arrowargs, **arrowkwargs)
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
         list = farray(list)
         out = Atoms(n=len(list), lattice=self.lattice)
         FortranAtoms.select(out, self, list=list)
      else:
         raise ValueError('Either mask or list must be present.')
      return out

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
      other.use_uniform_cutoff = self.use_uniform_cutoff
      other.cutoff = self.cutoff
      other.cutoff_break = self.cutoff_break
      other.nneightol = self.nneightol
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


   def density(self):
      from quippy import ElementMass, N_A, MASSCONVERT
      
      """Density in units of :math:`g/m^3`. If `mass` property exists, use that, otherwise we use `z`"""
      if hasattr(self, 'mass'):
         mass = sum(self.mass)/MASSCONVERT/1.0e3
      else:
         mass = sum(ElementMass[z] for z in self.z)/1.0e3

      return mass/(N_A*self.cell_volume()*1.0e-30)/1.0e3


   def add_property(self, name, value, n_cols=1):
      """
      Add a new property to this Atoms object.

      `name` is the name of the new property and `value` should be
      either a scalar or an array representing the value, which should
      be either integer, real, logical or string.

      If a scalar is given for `value` it is copied to every element
      in the new property.  `n_cols` can be specified to create a 2D
      property from a scalar initial value - the default is 1 which
      creates a 1D property.

      If an array is given for `value` it should either have shape
      (self.n,) for a 1D property or (n_cols,self.n) for a 2D
      property.  In this case `n_cols` is inferred from the shape of
      the `value` and shouldn't be passed as an argument.

      If property with the same type is already present then no error
      occurs. A warning is printed if the verbosity level is VERBOSE
      or higher. The value will be overwritten with that given in
      `value`.
      """
      
      if hasattr(value, 'shape'):
         # some kind of array:
         if len(value.shape) == 1:
            if value.shape[0] != self.n:
               raise ValueError('Bad array length for "value" - len(value.shape[0])=%d != self.n=%d'
                                % (value.shape[0], self.n))
            n_cols = 1
            value_ref = value[1]
         elif len(value.shape) == 2:
            if value.shape[1] != self.n:
               raise ValueError('Bad array length for "value" - len(value.shape[1])=%d != self.n=%d'
                                % (value.shape[1], self.n))
            value_ref = value[1,1]
            if value.dtype.kind == 'S':
               n_cols = 1
            else:
               n_cols = value.shape[0]
         else:
            raise ValueError('Bad array shape for "value" - should be either 1D or 2D')
      else:
         # some kind of scalar
         value_ref = value

      FortranAtoms.add_property(self, name, value_ref, n_cols)
      getattr(self, name.lower())[:] = value
            


from dictmixin import DictMixin, ParamReaderMixin
class Dictionary(DictMixin, ParamReaderMixin, FortranDictionary):

   def __init__(self, D=None, *args, **kwargs):
      FortranDictionary.__init__(self, *args, **kwargs)
      if 'keys' in self._arrays:
         self._arrays['_keys'] = self._arrays['keys']
         del self._arrays['keys']
         del self.keys
         self._update()
      if D is not None:
         self.read(D) # copy from D

   def keys(self):
      return [''.join(self._keys[:,i]).strip() for i in frange(self.n)]

   def __getitem__(self, k):
      i = self.lookup_entry_i(k)
      if i == -1: 
         raise KeyError('Key "%s" not found ' % k)

      t, s, s2 = self.get_type_and_size(k)

      from quippy import (T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                          T_CHAR, T_LOGICAL, T_INTEGER_A,
                          T_REAL_A, T_COMPLEX_A, T_CHAR_A,
                          T_LOGICAL_A, T_INTEGER_A2, T_REAL_A2)

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
         v = bool(v)
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
         v = farray(v, dtype=bool)
      elif t == T_INTEGER_A2:
         v,p = self._get_value_i_a2(k, s2[1], s2[2])
         v = farray(v)
      elif t == T_REAL_A2:
         v,p = self._get_value_r_a2(k, s2[1], s2[2])
         v = farray(v)
      else:
         raise ValueError('Unsupported dictionary entry type %d' % t)

      return v

   def __setitem__(self, k, v):
      self.set_value(k, v)

   def __delitem__(self, k):
      i = self.lookup_entry_i(k)
      if i == -1:
         raise KeyError('Key %s not found in Dictionary' % k)
      self.remove_entry(i)

   def __repr__(self):
      return ParamReaderMixin.__repr__(self)

   def __eq__(self, other):
      if sorted(self.keys()) != sorted(other.keys()): return False

      for key in self:
         v1, v2 = self[key], other[key]
         if isinstance(v1, FortranArray):
            if (v1 != v2).any(): return False
         else:
            if v1 != v2: return False
      return True
         
   def __ne__(self, other):
      return not self.__eq__(other)

   def __str__(self):
      return ParamReaderMixin.__str__(self)

   def copy(self):
      return Dictionary(self)


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


class DynamicalSystem(FortranDynamicalSystem):
   
   def run(self, pot, dt=1.0, n_steps=10, save_interval=1, connect_interval=10, args_str="",
           out=None, write_interval=10):
      self.atoms.calc_connect()
      pot.calc(self.atoms, args_str=args_str, calc_force=True, calc_energy=True)
      for n in range(n_steps):
         self.advance_verlet1(dt, self.atoms.force)
         pot.calc(self.atoms, args_str=args_str, calc_force=True, calc_energy=True)
         self.advance_verlet2(dt, self.atoms.force)
         self.print_status(epot=self.atoms.energy)
         self.atoms.params['time'] = self.t
         if connect_interval is not None and n % connect_interval == 0:
            self.atoms.calc_connect()
         if out is not None and n % write_interval == 0:
            out.write(self.atoms)
         if save_interval is not None and n % save_interval == 0:
            yield self.atoms.copy()
      raise StopIteration


from quippy import INPUT, OUTPUT, INOUT


class CInOutput(FortranCInOutput):

   def __init__(self, filename, action=INPUT, append=False, netcdf4=True, zero=False):
      FortranCInOutput.__init__(self, filename, action, append, netcdf4)
      self.zero = zero

   def __len__(self):
      return int(self.n_frame)

   def __getitem__(self, index):
      return self.read(frame=index, zero=self.zero)

   def __setitem__(self, index, value):
      self.write(value, frame=index)

   def write(self, at, properties=None, int_format=None, real_format=None, frame=None,
             shuffle=None, deflate=None, deflate_level=None, status=None):

      if properties is not None and (not hasattr(properties, 'dtype') or properties.dtype != dtype('S1')):
         properties = padded_str_array(properties, max([len(x) for x in properties])).T
      
      FortranCInOutput.write(self, at, properties, int_format, real_format, frame,
                             shuffle, deflate, deflate_level, status)
                             
      
             
             


