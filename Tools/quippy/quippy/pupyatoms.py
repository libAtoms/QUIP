from ordereddict import OrderedDict
from dictmixin import ParamReaderMixin
from farray import *

TABLE_STRING_LENGTH = 10
PROPERTY_INT = 1
PROPERTY_REAL = 2
PROPERTY_STR = 3
PROPERTY_LOGICAL = 4

class Dictionary(OrderedDict, ParamReaderMixin):
   """Subclass of OrderedDict for reading key/value pairs from strings or files.
      The original order of items is maintained. Values that looks like floats or ints
      or lists of floats or ints are automatically converted on reading."""
   
   def __init__(self, source=None):
      OrderedDict.__init__(self)
      if source is not None:
         self.read(source)

   def __repr__(self):
      return ParamReaderMixin.__repr__(self)

   def __str__(self):
      return ParamReaderMixin.__str__(self)

   def copy(self):
      return Dictionary(OrderedDict.copy(self))


class Table:
   """Pure python Table class"""

   def __init__(self, nint=0, nreal=0, nstr=0, nlogical=0, max_length=0):
      self.n = 0
      self.intsize = 0
      self.realsize = 0
      self.strsize = 0
      self.logicalsize = 0
      self.int = fzeros((0,0),int)
      self.real = fzeros((0,0),'d')
      self.str = fzeros((0,0,0),'S')
      self.logical = fzeros((0,0),bool)
      self.realloc(nint, nreal, nstr, nlogical, max_length)

   def allocate(self, nint=None, nreal=None, nstr=None, nlogical=None, max_length=None):

      if nint is not None:       self.intsize = nint
      if nreal is not None:      self.realsize = nreal
      if nstr is not None:       self.strsize = nstr
      if nlogical is not None:   self.logicalsize = nlogical
      if max_length is not None: self.max_length = max_length
      
      oldint = self.int[:]
      del self.int
      self.int = fzeros((self.intsize, self.n), int)
      if all([oldint.shape[i] > 0 for i in range(2)]):
         self.int[1:min(self.int.shape[0], oldint.shape[0]),
                  1:min(self.int.shape[1], oldint.shape[1])] = \
                  oldint[1:min(self.int.shape[0], oldint.shape[0]),
                         1:min(self.int.shape[1], oldint.shape[1])]

      oldreal = self.real[:]
      del self.real
      self.real = fzeros((self.realsize, self.n), 'd')
      if all([oldreal.shape[i] > 0 for i in range(2)]):
         self.real[1:min(self.real.shape[0], oldreal.shape[0]),
                   1:min(self.real.shape[1], oldreal.shape[1])] = \
                   oldreal[1:min(self.real.shape[0], oldreal.shape[0]),
                   1:min(self.real.shape[1], oldreal.shape[1])]

      oldstr = self.str[:]
      del self.str
      self.str = fzeros((TABLE_STRING_LENGTH, self.strsize, self.n), 'S')
      if all([oldstr.shape[i] > 0 for i in range(3)]):
         self.str[:,
                  1:min(self.str.shape[1], oldstr.shape[1]),
                  1:min(self.str.shape[2], oldstr.shape[2])] = \
                  oldstr[:,
                         1:min(self.str.shape[1], oldstr.shape[1]),
                         1:min(self.str.shape[2], oldstr.shape[2])]

      oldlogical = self.logical[:]
      del self.logical
      self.logical = fzeros((self.logicalsize, self.n), bool)
      if all([oldlogical.shape[i] > 0 for i in range(2)]):
         self.logical[1:min(self.logical.shape[0], oldlogical.shape[0]),
                      1:min(self.logical.shape[1], oldlogical.shape[1])] = \
                      oldlogical[1:min(self.logical.shape[0], oldlogical.shape[0]),
                                 1:min(self.logical.shape[1], oldlogical.shape[1])]
      


class Atoms:
   """Pure python Atoms class"""

   def __init__(self, filename=None, *allocargs, **allockwargs):

      self._atomsptr = None
      self.alloc(*allocargs, **allockwargs)

      if filename is not None:
         self.read(filename)

   def alloc(self, n=0, n_int=0, n_real=3, n_str=1, n_logical=0, use_libatoms=False, atomsptr=None, properties=None, \
                lattice=None, params=None,element='Si'):

      if use_libatoms or atomsptr is not None:
         if atomsptr is None:
            self.attach(libatoms.atoms_initialise(n,lattice))
         else:
            self.attach(atomsptr)
      else:
         self.n = n
         self.lattice = lattice
         self.g = numpy.linalg.inv(self.lattice)
         self.params = params

         # Create single property for atomic positions
         self.real    = numpy.zeros((self.n,n_real),dtype=float)
         self.int     = numpy.zeros((self.n,n_int),dtype=int)
         self.str     = numpy.zeros((self.n,n_str),dtype='S10')
         self.logical = numpy.zeros((self.n,n_logical),dtype=bool)

         if properties is None:
            self.properties = Dictionary({'species':('S',slice(0,1)),'pos':('R',slice(0,3))})
         else:
            self.properties = properties
            
         self.repoint()

   def __repr__(self):
      return 'Atoms(n=%d, properties=%s, params=%s, lattice=%s)' % \
            (self.n, repr(self.properties), repr(self.params), repr(self.lattice))


   def __cmp__(self, other):
       if other is None:
          return 1

       # Quick checks
       if (self.n != other.n) or (self.comment() != other.comment()):
          return 1

       # Check if arrays match one by one
       for this, that in \
           (self.lattice, other.lattice), \
           (self.real, other.real), (self.int, other.int), \
           (self.str, other.str), (self.logical, other.logical):

          if (not numpy.all(this == that)):
             return 1

       return 0

   def update(self, other):
      "Overwrite contents of this Atoms object with a copy of an other"
      self.n = other.n
      self.lattice = other.lattice.copy()
      self.g = other.g.copy()
      self.params = other.params.copy()
      self.properties = other.properties.copy()

      self.real    = other.real[:]
      self.int     = other.int[:]
      self.str     = other.str[:]
      self.logical = other.logical[:]

      self.repoint()


   def add_property(self, name, value, n_cols=1):
      "Add a new property to this Atoms object. Value can be a scalar int or float, or an array."

      # Scalar int or list of all ints
      if (type(value) == type(0)) or \
             ((type(value) == type([])) and numpy.all(numpy.array(map(type,value)) == type(0))):
         n_int = self.int.shape[1]
         intcopy = self.int.copy()
         self.int = numpy.zeros((self.n,n_int+n_cols),dtype=int)
         self.int[:,:n_int] = intcopy
         if n_cols == 1:
            self.int[:,n_int] = value
         else:
            self.int[:,n_int:n_int+n_cols] = value
         self.properties[name] = ('I',slice(n_int,n_int+n_cols))
         self.repoint()

      # Scalar real or list of all reals
      elif (type(value) == type(0.0)) or \
               (type(value) == type([]) and numpy.all(numpy.array(map(type,value)) == type(0.0))):
         n_real = self.real.shape[1]
         realcopy = self.real.copy()
         self.real = numpy.zeros((self.n,n_real+n_cols),dtype=float)
         self.real[:,:n_real] = realcopy
         if n_cols == 1:
            self.real[:,n_real] = value
         else:
            self.real[:,n_real:n_real+n_cols] = value
         self.properties[name] = ('R',slice(n_real,n_real+n_cols))
         self.repoint()

      # Scalar string or list of strings
      elif (type(value) == type('')) or \
             ((type(value) == type([])) and numpy.all(numpy.array(map(type,value)) == type(''))):
         n_str = self.str.shape[1]
         strcopy = self.str.copy()
         self.str = numpy.zeros((self.n,n_str+n_cols),dtype='S10')
         self.str[:,:n_str] = strcopy
         if n_cols == 1:
            self.str[:,n_str] = value
         else:
            self.str[:,n_str:n_str+n_cols] = value
         self.properties[name] = ('S',slice(n_str,n_str+n_cols))
         self.repoint()

      # Scalar logical or list of logicals
      elif (type(value) == type(False)) or \
             ((type(value) == type([])) and numpy.all(numpy.array(map(type,value)) == type(False))):
         n_logical = self.logical.shape[1]
         logicalcopy = self.logical.copy()
         self.logical = numpy.zeros((self.n,n_logical+n_cols),dtype=bool)
         self.logical[:,:n_logical] = logicalcopy
         if n_cols == 1:
            self.logical[:,n_logical] = value
         else:
            self.logical[:,n_logical:n_logical+n_cols] = value
         self.properties[name] = ('L',slice(n_logical,n_logical+n_cols))
         self.repoint()

      # Array type
      elif type(value) == type(numpy.array([])):
         if value.shape[0] != self.n:
            raise ValueError('length of value array (%d) != number of atoms (%d)' % \
                             (value.shape[0],self.n))

         if value.dtype.kind == 'f':
            try:
               n_cols = value.shape[1]
            except IndexError:
               n_cols = 1
            n_real = self.real.shape[1]
            realcopy = self.real.copy()
            self.real = numpy.zeros((self.n,n_real+n_cols),dtype=float)
            self.real[:,:n_real] = realcopy
            if n_cols == 1:
               self.real[:,n_real] = value.copy()
            else:
               self.real[:,n_real:n_real+n_cols] = value.copy()
            self.properties[name] = ('R',slice(n_real,n_real+n_cols))
            self.repoint()
         elif value.dtype.kind == 'i':
            try:
               n_cols = value.shape[1]
            except IndexError:
               n_cols = 1
            n_int = self.int.shape[1]
            intcopy = self.int.copy()
            self.int = numpy.zeros((self.n,n_int+n_cols),dtype=int)
            self.int[:,:n_int] = intcopy

            if n_cols == 1:
               self.int[:,n_int] = value.copy()
            else:
               self.int[:,n_int:n_int+n_cols] = value.copy()
            self.properties[name] = ('I',slice(n_int,n_int+n_cols))
            self.repoint()

         elif value.dtype.kind == 'S':
            try:
               n_cols = value.shape[1]
            except IndexError:
               n_cols = 1
            n_str = self.str.shape[1]
            strcopy = self.str.copy()
            self.str = numpy.zeros((self.n,n_str+n_cols),dtype='S10')
            self.str[:,:n_str] = strcopy

            if n_cols == 1:
               self.str[:,n_str] = value.copy()
            else:
               self.str[:,n_str:n_str+n_cols] = value.copy()
            self.properties[name] = ('S',slice(n_str,n_str+n_cols))
            self.repoint()

         elif value.dtype == numpy.dtype('bool'):
            try:
               n_cols = value.shape[1]
            except IndexError:
               n_cols = 1
            n_logical = self.logical.shape[1]
            logicalcopy = self.logical.copy()
            self.logical = numpy.zeros((self.n,n_logical+n_cols),dtype=numpy.dtype('bool'))
            self.logical[:,:n_logical] = logicalcopy

            if n_cols == 1:
               self.logical[:,n_logical] = value.copy()
            else:
               self.logical[:,n_logical:n_logical+n_cols] = value.copy()
            self.properties[name] = ('S',slice(n_logical,n_logical+n_cols))
            self.repoint()


         else:
            raise ValueError("Don't know how to add array property of type %r" % value.dtype)
         
      else:
         raise ValueError("Don't know how to add property of type %r" % type(value))
         

   def repoint(self):
      "Make pointers to columns in real and int"

      for prop, (ptype, cols) in self.properties.items():
         if ptype == 'R':
            if cols.stop-cols.start == 1:
               setattr(self, prop, self.real[:,cols.start])
            else:
               setattr(self, prop, self.real[:,cols])
         elif ptype == 'I':
            if cols.stop-cols.start == 1:
               setattr(self, prop, self.int[:,cols.start])
            else:
               setattr(self, prop, self.int[:,cols])
         elif ptype == 'S':
            if cols.stop-cols.start == 1:
               setattr(self, prop, self.str[:,cols.start])
            else:
               setattr(self, prop, self.str[:,cols])
         elif ptype == 'L':
            if cols.stop-cols.start == 1:
               setattr(self, prop, self.logical[:,cols.start])
            else:
               setattr(self, prop, self.logical[:,cols])
         else:
            raise ValueError('Bad property type :'+str(self.properties[prop]))


   def filter(self, mask):
      "Return smaller Atoms with only the elements where mask is true"

      other = Atoms()

      if mask is None:
         mask = numpy.zeros((self.n,),numpy.bool)
         mask[:] = True

      other.n = count(mask)
      other.lattice = self.lattice.copy()
      other.g = self.g.copy()
      other.params = self.params.copy()
      other.properties = self.properties.copy()

      other.real = self.real[mask]
      other.int  = self.int[mask]
      other.str  = self.str[mask]
      other.logical = self.logical[mask]

      other.repoint()

      return other

   def copy(self):
      if self.n == 0:
         return Atoms()
      else:
         return self.filter(mask=None)

   def add(self, newpos, newspecies):

      if type(newpos) == type([]):
         newpos = numpy.array(newpos)

      if len(newpos.shape) == 1:
         n_new = 1
      else:
         n_new = newpos.shape[0]

      oldn = self.n
      self.n = self.n + n_new

      self.real    = numpy.resize(self.real,    (self.n,self.real.shape[1]))
      self.int     = numpy.resize(self.int,     (self.n,self.int.shape[1]))
      self.str     = numpy.resize(self.str,     (self.n,self.str.shape[1]))
      self.logical = numpy.resize(self.logical, (self.n,self.logical.shape[1]))

      self.repoint()

      self.pos[oldn:self.n] = newpos
      self.species[oldn:self.n] = newspecies


   def remove(self, discard):
      keep = [ i for i in range(self.n) if not i in discard ]
               
      self.n = len(keep)
      self.real = self.real[keep]
      self.int  = self.int[keep]
      self.str  = self.str[keep]
      self.logical = self.logical[keep]
      self.repoint()

def diamond(a, z):
   "Bulk cube of element with lattice constant a"
   at = Atoms(n=8)

   at.z[:] = z
   at.species[:] = species
   at.lattice = numpy.array([[a,0.,0.],[0.,a,0.],[0.,0.,a]])
   at.g = numpy.linalg.inv(at.lattice)
   at.pos[0,:] = a*numpy.array([0.00, 0.00, 0.00])
   at.pos[1,:] = a*numpy.array([0.25, 0.25, 0.25])
   at.pos[2,:] = a*numpy.array([0.50, 0.50, 0.00])
   at.pos[3,:] = a*numpy.array([0.75, 0.75, 0.25])
   at.pos[4,:] = a*numpy.array([0.50, 0.00, 0.50])
   at.pos[5,:] = a*numpy.array([0.75, 0.25, 0.75])
   at.pos[6,:] = a*numpy.array([0.00, 0.50, 0.50])
   at.pos[7,:] = a*numpy.array([0.25, 0.75, 0.75])

   return at


def supercell(self, n1, n2, n3):

   other = Atoms(n=self.n*n1*n2*n3,lattice=self.lattice,
                 properties=self.properties.copy())

   other.lattice[1,:] = self.lattice[1,:]*n1
   other.lattice[2,:] = self.lattice[2,:]*n2
   other.lattice[3,:] = self.lattice[3,:]*n3
   other.g = numpy.linalg.inv(other.lattice)

   for i in range(n1):
      for j in range(n2):
         for k in range(n3):
            p = numpy.dot(self.lattice, numpy.array([i,j,k]))
            for n in range(self.n):
               nn = ((i*n2+j)*n3+k)*self.n + n
               other.int[nn,:] = self.int[n,:]
               other.real[nn,:] = self.real[n,:]
               other.logical[nn,:] = self.logical[n,:]
               other.str[nn,:] = self.str[n,:]
               other.pos[nn,:] = self.pos[n,:] + p

   other.repoint()
   return other



