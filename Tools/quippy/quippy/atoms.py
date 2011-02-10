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

import numpy, os, weakref
from farray import *
from quippy import FortranAtoms
from quippy import AtomsReaders, AtomsWriters
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

   def __init__(self, at, hysteretic=False):
      self.atref = weakref.ref(at)
      self.hysteretic = hysteretic

   def __iter__(self):
      return self.iterneighbours()

   def iterneighbours(self):
      at = self.atref()

      for i in frange(at.n):
         yield self[i]

   def __getitem__(self, i):

      at = self.atref()
      if self.hysteretic:
         connect = at.hysteretic_connect
      else:
         connect = at.connect

      if not connect.initialised:
         if self.hysteretic:
            at.calc_connect_hysteretic(connect)
         else:
            at.calc_connect(connect)

      distance = farray(0.0)
      diff = fzeros(3)
      cosines = fzeros(3)
      shift = fzeros(3,dtype=numpy.int32)

      res = []
      for n in frange(at.n_neighbours(i,alt_connect=connect)):
         j = at.neighbour(i, n, distance, diff, cosines, shift, alt_connect=connect)
         res.append(NeighbourInfo(j,distance,diff,cosines,shift))

      return farray(res) # to give 1-based indexing      

   def distances(self, Z1=None, Z2=None):
      """Distances between pairs of neighbours, optionally filtered by species (Z1,Z2)"""
      at = self.atref()
      for i in frange(at.n):
         for neighb in self[i]:
            if neighb.j > i: continue
            if Z1 is not None and Z2 is not None:
               if sorted((at.z[i],at.z[neighb.j])) == sorted((Z1,Z2)):
                  yield neighb.distance
            else:
               yield neighb.distance

class Atom(dict):
   _cmp_tol = 1e-8
   
   def __init__(self, *args, **kwargs):
      dict.__init__(self, *args, **kwargs)

   def __getattr__(self, key):
      return self[key]

   def __setattr__(self, key, value):
      raise ValueError('Atom instances are read-only')
      
   def __repr__(self):
      return 'Atom(%s)' % ', '.join(['%s=%r' % (k,v) for (k, v) in self.iteritems()])

   def __str__(self):
      title = 'Atom %d' % self['i']
      title = title + '\n' + '-'*len(title)+'\n\n'
      return title+'\n'.join(['%-15s =  %s' % (k,self[k]) for k in sorted(self.keys()) if k != 'i'])

   def __eq__(self, other):
      self_ = dict([(k.lower(),v) for (k,v) in self.iteritems() if k != 'i'])
      other_ = dict([(k.lower(),v) for (k,v) in other.iteritems() if k != 'i'])

      if sorted(self_.keys()) != sorted(other_.keys()): return False
      
      for k in self_.keys():
         v1, v2 = self_[k], other_[k]
         if hasattr(v1, '__iter__') and hasattr(v2, '__iter__'):
            v1 = farray(v1)
            v2 = farray(v2)
            if v1.dtype.kind != 'f':
               if (v1 != v2).any(): return False
            else:
               if abs(v1 - v2).max() > self._cmp_tol: return False
         else:
            if v1 != v2:
               return False
      return True

   def __ne__(self, other):
      return not self.__eq__(other)


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
   _cmp_skip_fields = ['own_this', 'ref_count', 'domain']

   def __new__(cls, source=None, n=0, lattice=fidentity(3), fpointer=None, finalise=True,
                properties=None, params=None, fixed_size=False, *readargs, **readkwargs):
      """
      Initialise an Atoms object.

      The arguments, n, lattice, fpointer, finalise, data, properties and params are passed on
      to the FortranAtoms constructor.

      If source is given, Atoms.read(source, frame=frame) is called.
      """

      if source is None:
         self = object.__new__(cls)
      elif isinstance(source, cls):
         self = cls.copy(source)
      else:
         self = cls.read(source, *readargs, **readkwargs)

      return self

   def __init__(self, source=None, n=0, lattice=fidentity(3), fpointer=None, finalise=True,
                properties=None, params=None, fixed_size=False, *readargs, **readkwargs):
      from quippy import Dictionary
      if source is None:
         if params is not None and not isinstance(params, Dictionary):
            params = Dictionary(params)
         if properties is not None and not isinstance(properties, Dictionary):
            properties = Dictionary(properties)
         FortranAtoms.__init__(self, n=n, lattice=lattice, fpointer=fpointer, finalise=finalise,
                               properties=properties, params=params, fixed_size=fixed_size)
         self.neighbours = Neighbours(self)
         self.hysteretic_neighbours = Neighbours(self,hysteretic=True)

   @classmethod
   def read(cls, source, format=None, *args, **kwargs):
      if format is None:
         if isinstance(source, str):
            if source in AtomsReaders:
               format = source
            else:
               source = os.path.expanduser(source)
               base, ext = os.path.splitext(source)
               format = ext[1:]
         else:
            format = source.__class__

      opened = False
      if format in AtomsReaders:
         source = AtomsReaders[format](source, *args, **kwargs)
         opened = True

      if isinstance(source, str):
         raise IOError("Don't know how to read Atoms from file '%s'" % source)
      if not hasattr(source, '__iter__'):
         raise IOError('Cannot read from source %r - not an iterator' % source)
      
      at = iter(source).next()
      if not isinstance(at, cls):
         raise ValueError('Object %r read from source %r is not Atoms instance' % (at, source))
      if opened and hasattr(source, 'close'):
         source.close()
      return at
         

   def write(self, dest, format=None, properties=None, prefix=None, *args, **kwargs):
      opened = False
      if format is None:
         if isinstance(dest, str):
            opened = True
            if dest in AtomsWriters:
               format = dest
            else:
               dest = os.path.expanduser(dest)
               base, ext = os.path.splitext(dest)
               format = ext[1:]
         else:
            format = dest.__class__

      if format in AtomsWriters:
         dest = AtomsWriters[format](dest, *args, **kwargs)

      if not hasattr(dest, 'write'):
         raise ValueError("Don't know how to write to destination \"%s\" in format \"%s\"" % (dest, format))

      write_kwargs = {}
      if properties is not None: write_kwargs['properties'] = properties
      if prefix is not None: write_kwargs['prefix'] = prefix
      try:
         res = dest.write(self, **write_kwargs)
      except TypeError:
         raise ValueError('Atoms.write destination %r does not support arguments %r' % (dest, write_kwargs))

      if opened and hasattr(dest, 'close'):
         dest.close()
      return res


   def show(self, *args, **kwargs):
      """Show this Atoms object in AtomEye."""
      try:
         import atomeye
         atomeye.show(self, *args, **kwargs)
      except ImportError:
         raise RuntimeError('AtomEye not available')

   def select(self, mask=None, list=None, orig_index=None):
      """Select a subset of the atoms in an atoms object

      Use either a logical mask array or a list of atom indices to include.

      small_at = at.select([mask, list])
      
      """
      if mask is not None:
         mask = farray(mask)
         out = Atoms(n=mask.count(),lattice=self.lattice)
         FortranAtoms.select(out, self, mask=mask, orig_index=orig_index)
      elif list is not None:
         list = farray(list)
         out = Atoms(n=len(list), lattice=self.lattice)
         FortranAtoms.select(out, self, list=list, orig_index=orig_index)
      else:
         raise ValueError('Either mask or list must be present.')
      return out

   def copy(self):
      other = Atoms(n=self.n, lattice=self.lattice, 
                    properties=self.properties, params=self.params)
      other.use_uniform_cutoff = self.use_uniform_cutoff
      other.cutoff = self.cutoff
      other.cutoff_break = self.cutoff_break
      other.nneightol = self.nneightol
      return other

   def __len__(self):
      return self.n

   def __getattr__(self, name):
      try:
         return self.properties[name]
      except KeyError:
         try:
            return self.params[name]
         except KeyError:
            raise AttributeError('Unknown Atoms attribute %s' % name)

   def __getitem__(self, i):
      if i < 1 or i > self.n:
         raise ValueError('Atoms index should be in range 1..self.n(%d)' % self.n)
      res = {}
      res['i'] = i
      for k in self.properties.keys():
         v = getattr(self,k.lower())[...,i]
         if isinstance(v,FortranArray):
            if v.dtype.kind == 'S':
               v = ''.join(v).strip()
            elif v.shape == ():
               v = v.item()
         res[k.lower()] = v

      return Atom(**res)

   def __setitem__(self, i, atom):
      if i < 1 or i > self.n:
         raise ValueError('Atoms index should be in range 1..self.n(%d)' % self.n)
      for k, v in atom.iteratoms():
         setattr(self, k.lower())[...,i] = v

      
   def density(self):
      from quippy import ElementMass, N_A, MASSCONVERT
      
      """Density in units of :math:`g/m^3`. If `mass` property exists, use that, otherwise we use `z`"""
      if hasattr(self, 'mass'):
         mass = sum(self.mass)/MASSCONVERT/1.0e3
      else:
         mass = sum(ElementMass[z] for z in self.z)/1.0e3

      return mass/(N_A*self.cell_volume()*1.0e-30)/1.0e3


   def add_property(self, name, value, n_cols=None, overwrite=None, property_type=None):
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

      If `property_type` is present, then no attempt is made to
      infer the type from `value`. This is necessary to resolve
      ambiguity between integer and logical types.

      If property with the same type is already present then no error
      occurs.If `overwrite` is true, the value will be overwritten with
      that given in `value`, otherwise the old value is retained.
      """

      kwargs = {}
      if n_cols is not None: kwargs['n_cols'] = n_cols
      if overwrite is not None: kwargs['overwrite'] = overwrite

      if property_type is None:
         FortranAtoms.add_property(self, name, value, **kwargs)

      else:
         # override value_ref if property_type is specified
         
         from quippy import (T_INTEGER_A, T_REAL_A, T_LOGICAL_A, T_CHAR_A,
                             T_INTEGER_A2, T_REAL_A2, TABLE_STRING_LENGTH)

         new_property = not self.has_property(name)
         
         type_to_value_ref = {
            T_INTEGER_A  : 0,
            T_REAL_A : 0.0,
            T_CHAR_A  : " "*TABLE_STRING_LENGTH,
            T_LOGICAL_A : False,
            T_INTEGER_A2 : 0,
            T_REAL_A2: 0.0
            }
         try:
            value_ref = type_to_value_ref[property_type]
         except KeyError:
            raise ValueError('Unknown property_type %d' % property_type)

         if hasattr(value, 'shape') and len(value.shape) == 2 and property_type != T_CHAR_A and n_cols is None:
            kwargs['n_cols'] = value.shape[0]

         FortranAtoms.add_property(self, name, value_ref, **kwargs)
         if new_property or overwrite:
            getattr(self, name.lower())[:] = value            



