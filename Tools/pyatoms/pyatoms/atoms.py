# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import sys, string, operator, os
from ordereddict import OrderedDict
from paramreader import ParamReader
import numpy
import castep

#from f2pyatoms import libatoms, arraydata
#libatoms.system_initialise()

BOHR_TO_ANG = 0.529177249
HARTREE_TO_EV = 27.2107

# List of elements, in order of increasing atomic number
ElementName = map(string.strip,\
[ "xx", "H  ","He ","Li ","Be ","B  ","C  ","N  ","O  ","F  ","Ne ","Na ","Mg ","Al ","Si ","P  ","S  ", \
   "Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ","Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ","Ge ",\
   "As ","Se ","Br ","Kr ","Rb ","Sr ","Y  ","Zr ","Nb ","Mo ","Tc ","Ru ","Rh ","Pd ","Ag ","Cd ",\
   "In ","Sn ","Sb ","Te ","I  ","Xe ","Cs ","Ba ","La ","Ce ","Pr ","Nd ","Pm ","Sm ","Eu ","Gd ",\
   "Tb ","Dy ","Ho ","Er ","Tm ","Yb ","Lu ","Hf ","Ta ","W  ","Re ","Os ","Ir ","Pt ","Au ","Hg ",\
   "Tl ","Pb ","Bi ","Po ","At ","Rn ","Fr ","Ra ","Ac ","Th ","Pa ","U  ","Np ","Pu ","Am ","Cm ",\
   "Bk ","Cf ","Es ","Fm ","Md ","No ","Lr ","Rf ","Db ","Sg ","Bh ","Hs ","Mt ","Ds ","Rg ","Uub",\
   "Uut","Uuq","Uup","Uuh"])

# Temp array to intialise ElementMass dictionary
mass = \
[ 0.0, 1.00794, 4.00260, 6.941, 9.012187, 10.811, 12.0107, 14.00674, 15.9994, 18.99840, 20.1797, 22.98977, \
24.3050, 26.98154, 28.0855, 30.97376, 32.066, 35.4527, 39.948, 39.0983, 40.078, 44.95591, 47.867,     \
50.9415, 51.9961, 54.93805, 55.845, 58.93320, 58.6934, 63.546, 65.39, 69.723, 72.61, 74.92160, 78.96, \
79.904, 83.80, 85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550, 106.42,    \
107.8682, 112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.29, 132.90545, 137.327, 138.9055,\
140.116, 140.90765, 144.24, 145.0, 150.36, 151.964, 157.25, 158.92534, 162.50, 164.93032, 167.26,     \
168.93421, 173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,   \
200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0381, 231.03588,    \
238.0289, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0, 261.0, 262.0,  \
263.0, 264.0, 265.0, 268.0, 271.0, 272.0, 285.0, 284.0, 289.0, 288.0, 292.0 ]

# Mapping of element name to mass
ElementMass = dict(zip(ElementName, mass))
del mass

class Atoms:
   "Class to deal with a single frame of an xyz movie"

   def __init__(self, filename=None, *allocargs, **allockwargs):

      self._atomsptr = None
      self.alloc(*allocargs, **allockwargs)

      if filename is not None:
         self.read(filename)

   def alloc(self, n=0, n_int=0, n_real=3, n_str=1, n_logical=0, use_libatoms=False, atomsptr=None, properties=None, \
                lattice=numpy.array([[100.,0.,0.],[0.,100.,0.],[0.,0.,100.]]), \
                params=ParamReader(),element='Si'):

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
            self.properties = OrderedDict({'species':('S',slice(0,1)),'pos':('R',slice(0,3))})
         else:
            self.properties = properties
            
         self.repoint()

   def attach(self, atomsptr):
      self.finalise()
      self._atomsptr = atomsptr
      
      self.n, n_int, n_real, n_str, n_logical, iloc, rloc, sloc, lloc, latticeloc, gloc = \
                 libatoms.atoms_get_data(self._atomsptr)

      self.int     = arraydata((self.n,n_int),int,iloc)
      self.real    = arraydata((self.n,n_real),float,rloc)
      self.str     = arraydata((self.n,n_str),'S10',sloc)
      self.logical = arraydata((self.n,n_logical),bool,sloc)

      self.lattice = arraydata((3,3),float,latticeloc)
      self.g       = arraydata((3,3),float,gloc)
         
      self.params  = {}

      property_code_map = {1:'I', 2: 'R', 3:'S', 4:'L'}
      self.properties = OrderedDict()
      for i in range(libatoms.atoms_n_properties(self._atomsptr)):
         key, (code, startcol, stopcol) = libatoms.atoms_nth_property(self._atomsptr, i+1)
         self.properties[key.strip()] = (property_code_map[code],slice(startcol-1,stopcol))

      self.repoint()


   def finalise(self):
      if self._atomsptr is not None:
         libatoms.atoms_finalise(self._atomsptr)
         self._atomsptr = None


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


   def add_property(self, name, value, ncols=1):
      "Add a new property to this Atoms object. Value can be a scalar int or float, or an array."

      # Scalar int or list of all ints
      if (type(value) == type(0)) or \
             ((type(value) == type([])) and numpy.all(numpy.array(map(type,value)) == type(0))):
         n_int = self.int.shape[1]
         intcopy = self.int.copy()
         self.int = numpy.zeros((self.n,n_int+ncols),dtype=int)
         self.int[:,:n_int] = intcopy
         if ncols == 1:
            self.int[:,n_int] = value
         else:
            self.int[:,n_int:n_int+ncols] = value
         self.properties[name] = ('I',slice(n_int,n_int+ncols))
         self.repoint()

      # Scalar real or list of all reals
      elif (type(value) == type(0.0)) or \
               (type(value) == type([]) and numpy.all(numpy.array(map(type,value)) == type(0.0))):
         n_real = self.real.shape[1]
         realcopy = self.real.copy()
         self.real = numpy.zeros((self.n,n_real+ncols),dtype=float)
         self.real[:,:n_real] = realcopy
         if ncols == 1:
            self.real[:,n_real] = value
         else:
            self.real[:,n_real:n_real+ncols] = value
         self.properties[name] = ('R',slice(n_real,n_real+ncols))
         self.repoint()

      # Scalar string or list of strings
      elif (type(value) == type('')) or \
             ((type(value) == type([])) and numpy.all(numpy.array(map(type,value)) == type(''))):
         n_str = self.str.shape[1]
         strcopy = self.str.copy()
         self.str = numpy.zeros((self.n,n_str+ncols),dtype='S10')
         self.str[:,:n_str] = strcopy
         if ncols == 1:
            self.str[:,n_str] = value
         else:
            self.str[:,n_str:n_str+ncols] = value
         self.properties[name] = ('S',slice(n_str,n_str+ncols))
         self.repoint()

      # Scalar logical or list of logicals
      elif (type(value) == type(False)) or \
             ((type(value) == type([])) and numpy.all(numpy.array(map(type,value)) == type(False))):
         n_logical = self.logical.shape[1]
         logicalcopy = self.logical.copy()
         self.logical = numpy.zeros((self.n,n_logical+ncols),dtype=bool)
         self.logical[:,:n_logical] = logicalcopy
         if ncols == 1:
            self.logical[:,n_logical] = value
         else:
            self.logical[:,n_logical:n_logical+ncols] = value
         self.properties[name] = ('L',slice(n_logical,n_logical+ncols))
         self.repoint()

      # Array type
      elif type(value) == type(numpy.array([])):
         if value.shape[0] != self.n:
            raise ValueError('length of value array (%d) != number of atoms (%d)' % \
                             (value.shape[0],self.n))

         if value.dtype.kind == 'f':
            try:
               ncols = value.shape[1]
            except IndexError:
               ncols = 1
            n_real = self.real.shape[1]
            realcopy = self.real.copy()
            self.real = numpy.zeros((self.n,n_real+ncols),dtype=float)
            self.real[:,:n_real] = realcopy
            if ncols == 1:
               self.real[:,n_real] = value.copy()
            else:
               self.real[:,n_real:n_real+ncols] = value.copy()
            self.properties[name] = ('R',slice(n_real,n_real+ncols))
            self.repoint()
         elif value.dtype.kind == 'i':
            try:
               ncols = value.shape[1]
            except IndexError:
               ncols = 1
            n_int = self.int.shape[1]
            intcopy = self.int.copy()
            self.int = numpy.zeros((self.n,n_int+ncols),dtype=int)
            self.int[:,:n_int] = intcopy

            if ncols == 1:
               self.int[:,n_int] = value.copy()
            else:
               self.int[:,n_int:n_int+ncols] = value.copy()
            self.properties[name] = ('I',slice(n_int,n_int+ncols))
            self.repoint()

         elif value.dtype.kind == 'S':
            try:
               ncols = value.shape[1]
            except IndexError:
               ncols = 1
            n_str = self.str.shape[1]
            strcopy = self.str.copy()
            self.str = numpy.zeros((self.n,n_str+ncols),dtype='S10')
            self.str[:,:n_str] = strcopy

            if ncols == 1:
               self.str[:,n_str] = value.copy()
            else:
               self.str[:,n_str:n_str+ncols] = value.copy()
            self.properties[name] = ('S',slice(n_str,n_str+ncols))
            self.repoint()

         elif value.dtype == numpy.dtype('bool'):
            try:
               ncols = value.shape[1]
            except IndexError:
               ncols = 1
            n_logical = self.logical.shape[1]
            logicalcopy = self.logical.copy()
            self.logical = numpy.zeros((self.n,n_logical+ncols),dtype=numpy.dtype('bool'))
            self.logical[:,:n_logical] = logicalcopy

            if ncols == 1:
               self.logical[:,n_logical] = value.copy()
            else:
               self.logical[:,n_logical:n_logical+ncols] = value.copy()
            self.properties[name] = ('S',slice(n_logical,n_logical+ncols))
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


   def comment(self, properties=None):
      "Return the comment line for this Atoms object"

      if properties is None:
         props = self.properties.keys()
      else:
         props = properties
         
      lattice_str = 'Lattice="' + ' '.join(map(str, numpy.reshape(self.lattice,9))) + '"'
      
      props_str = 'Properties=' + ':'.join(map(':'.join, \
               zip(props, \
                   [self.properties[k][0] for k in props], \
                   [str(self.properties[k][1].stop-self.properties[k][1].start) for k in props])))

      return lattice_str+' '+props_str+' '+str(self.params)


   def _props_dtype(self, props=None):
      "Return a record array dtype for the specified properties (default all)"

      if props is None:
         props = self.properties.keys()

      result = []
      fmt_map = {'R':'d','I':'i','S':'S10','L':'bool'}

      for prop in props:
         ptype, cols = self.properties[prop]
         if cols.start == cols.stop-1:
            result.append((prop,fmt_map[ptype]))
         else:
            for c in range(cols.stop-cols.start):
               result.append((prop+str(c),fmt_map[ptype]))

      return numpy.dtype(result)


   def to_recarray(self, props=None):
      "Return a record array contains specified properties in order (defaults to all properties)"

      if props is None:
         props = self.properties.keys()

      # Create empty record array with correct dtype
      data = numpy.zeros(self.n,self._props_dtype(props))

      # Copy cols from self.real and self.int into data recarray
      for prop in props:
         ptype, cols = self.properties[prop]
         if ptype == 'R':
            if cols.start == cols.stop-1:
               data[prop] = self.real[:,cols.start]
            else:
               for c in range(cols.stop-cols.start):
                  data[prop+str(c)] = self.real[:,cols.start+c]
         elif ptype == 'I':
            if cols.start == cols.stop-1:
               data[prop] = self.int[:,cols.start]
            else:
               for c in range(cols.stop-cols.start):
                  data[prop+str(c)] = self.int[:,cols.start+c]
         elif ptype == 'S':
            if cols.start == cols.stop-1:
               data[prop] = self.str[:,cols.start]
            else:
               for c in range(cols.stop-cols.start):
                  data[prop+str(c)] = self.str[:,cols.start+c]
         elif ptype == 'L':
            if cols.start == cols.stop-1:
               data[prop] = self.logical[:,cols.start]
            else:
               for c in range(cols.stop-cols.start):
                  data[prop+str(c)] = self.logical[:,cols.start+c]
         else:
            raise ValueError('Bad property type :'+str(self.properties[prop][1]))

      return data


   def update_from_recarray(self, data, props=None):
      """Update Atoms data from a record array. By default all properties
      are updated; use the props argument to update only a subset"""

      if props is None:
         props = self.properties.keys()

      if data.dtype != self._props_dtype(props) or data.shape != (self.n,):
         raise ValueError('Data shape is incorrect')

      # Copy cols from data recarray into self.real and self.int
      for prop in props:
         ptype, cols = self.properties[prop]
         if ptype == 'R':
            if cols.start == cols.stop-1:
               self.real[:,cols.start] = data[prop]
            else:
               for c in range(cols.stop-cols.start):
                  self.real[:,cols.start+c] = data[prop+str(c)]
         elif ptype == 'I':
            if cols.start == cols.stop-1:
               self.int[:,cols.start] = data[prop]
            else:
               for c in range(cols.stop-cols.start):
                  self.int[:,cols.start+c] = data[prop+str(c)]
         elif ptype == 'S':
            if cols.start == cols.stop-1:
               self.str[:,cols.start] = data[prop]
            else:
               for c in range(cols.stop-cols.start):
                  self.str[:,cols.start+c] = data[prop+str(c)]
         elif ptype == 'L':
            if cols.start == cols.stop-1:
               self.logical[:,cols.start] = data[prop]
            else:
               for c in range(cols.stop-cols.start):
                  self.logical[:,cols.start+c] = data[prop+str(c)]
         else:
            raise ValueError('Bad property type :'+str(self.properties[prop][1]))


   def read_xyz(self, xyz):
      "Read from extended XYZ filename or open file."

      opened = False
      if type(xyz) == type(''):
         xyz = open(xyz,'r')
         opened = True

      line = xyz.next()
      if not line: return False

      n = int(line.strip())
      comment = (xyz.next()).strip()

      # Parse comment line
      params = ParamReader(comment)

      if not 'Properties' in params:
         raise ValueError('Properties missing from comment line')

      properties, n_int, n_real, n_str, n_logical = _parse_properties(params['Properties'])
      del params['Properties']

      # Get lattice
      if not 'Lattice' in params:
         raise ValueError('No lattice found in xyz file')

      lattice = numpy.reshape(params['Lattice'], (3,3))
      del params['Lattice']

      self.alloc(n=n,lattice=lattice,properties=properties,params=params,\
                 n_int=n_int,n_real=n_real,n_str=n_str,n_logical=n_logical)

      props_dtype = self._props_dtype()
      
      converters = [_getconv(props_dtype.fields[name][0]) \
                    for name in props_dtype.names]

      X = []
      for i,line in enumerate(xyz):
         vals = line.split()
         row = tuple([converters[j](val) for j,val in enumerate(vals)])
         X.append(row)
         if i == self.n-1: break # Only read self.n lines

      try:
         data = numpy.array(X,props_dtype)
      except TypeError:
         raise IOError('End of file reached before end of frame')

      if opened: xyz.close()

      try:
         self.update_from_recarray(data)
      except ValueError:
         # got a partial frame, must be end of file
         return False
      else:
         return True

   def read_netcdf(self, fname, frame=0):
      from pupynere import netcdf_file

      nc = netcdf_file(fname)
      
      self.n = nc.dimensions['atom']
      self.lattice = make_lattice(nc.variables['cell_lengths'][frame], 
                                  nc.variables['cell_angles'][frame])
      self.g = numpy.linalg.inv(self.lattice)
      self.params = OrderedDict()
      self.properties = OrderedDict()

      self.real    = numpy.zeros((self.n,0),dtype=float)
      self.int     = numpy.zeros((self.n,0),dtype=int)
      self.str     = numpy.zeros((self.n,0),dtype='S10')
      self.logical = numpy.zeros((self.n,0),dtype=bool)

      vars = nc.variables.keys()
      vars = filter(lambda v: not v in ('cell_angles', 'cell_lengths'), vars)

      # ensure first var is species and second positions
      sp = vars.index('species')
      if sp != 0:
         vars[sp], vars[0] = vars[0], vars[sp]
      pos = vars.index('coordinates')
      if pos != 1:
         vars[pos], vars[1] = vars[1], vars[pos]

      for v in vars:
         d = nc.variables[v].dimensions

         if d[0] != 'frame': continue

         value = nc.variables[v][frame]
         if value.dtype == numpy.dtype('|S1'):
            value = [''.join(x).strip() for x in value]

         if len(d) == 1 or (len(d) == 2 and d[1] in ('label','string')):
            if (len(d) == 2 and d[1] in ('label','string')):
               value = ''.join(value)
            self.params[v] = value
         else:
            # Name mangling
            if v == 'coordinates': 
               p = 'pos'
            elif v == 'velocities': 
               p = 'velo'
            else:
               p = v
            value = nc.variables[v][frame]
            if value.dtype == numpy.dtype('|S1'):
               value = [''.join(x).strip() for x in value]
            self.add_property(p, value)


   def write_xyz(self, xyz=sys.stdout, properties=None):
      "Write atoms in extended XYZ format. xyz can be a filename or open file"

      if properties is None:
         # Sort by original order
         props = self.properties.keys()
      else:
         props = properties

      species = getattr(self,props[0])
      if len(species.shape) != 1 or species.dtype.kind != 'S':
         raise ValueError('First property must be species like')

      pos = getattr(self,props[1])
      if pos.shape[1] != 3 or pos.dtype.kind != 'f':
         raise ValueError('Second property must be position like')

      data = self.to_recarray(props)
      format = ''.join([_getfmt(data.dtype.fields[name][0]) for name in data.dtype.names])+'\n'

      opened = False
      if type(xyz) == type(''):
         xyz = open(xyz, 'w')
         opened = True
      
      xyz.write('%d\n' % self.n)
      xyz.write(self.comment(properties)+'\n')
      for i in range(self.n):
         xyz.write(format % tuple(data[i]))

      if opened: xyz.close()


   def read_cell(self, cell):
      "Read atoms from a CastepCell object or file"
      
      if hasattr(cell, 'next'): # looks like a file
         cell = castep.CastepCell(cell)

      self.update(cell.to_atoms())


   def write_cell(self, fname):
      "Write Atoms to a cell file"

      cell = castep.CastepCell()
      cell.update_from_atoms(self)
      cell.write(fname)


   def read_geom(self, geom):
      "Read from a CASTEP .geom file"
      self.update(castep.read_geom(geom))
      

   def read_castep(self, castepfile):
      "Read from a .castep output file"

      if self.n != 0:
         self.update(castep.read_castep_output(castepfile, self, abort=False))
      else:
         self.update(castep.read_castep_output(castepfile, abort=False))


   def read(self, fname, filetype=None):
      "Attempt to guess type of file from extension and call appropriate read method"

      opened = False
      if type(fname) == type(''):
         if fname.endswith('.gz'):
            import gzip
            fh = gzip.open(fname)
            fname = fname[:-3] # remove .gz
         elif fname.endswith('.nc'):
            fh = fname
         else:
            fh = open(fname, 'r')
            opened = True

         # Guess file type from extension
         if filetype is None:
            root, filetype = os.path.splitext(fname)
            filetype = filetype[1:] # remove '.'
      else:
         fh = fname

      # Default to xyz format
      if not filetype in ['cell','geom','xyz','castep','nc']:
         filetype = 'xyz'

      if filetype == 'xyz':
         self.read_xyz(fh)
      elif filetype == 'cell':
         self.read_cell(fh)
      elif filetype == 'geom':
         self.read_geom(fh)
      elif filetype == 'castep':
         self.read_castep(fh)
      elif filetype == 'nc':
         self.read_netcdf(fh)

      if opened: fh.close()

      
   def write(self, fname, filetype=None):
      opened = False
      if type(fname) == type(''):
         if fname.endswith('.gz'):
            import gzip
            fh = gzip.open(fname,'w')
            fname = fname[:-3] # remove .gz
         else:
            fh = open(fname, 'w')

         # Guess file type from extension
         if filetype is None:
            root, filetype = os.path.splitext(fname)
            filetype = filetype[1:] # remove '.'
         opened = True
      else:
         fh = fname

      # Default to xyz format
      if not filetype in ['xyz','cfg','cell']:
         filetype = 'xyz'

      if filetype == 'xyz':
         self.write_xyz(fh)
      elif filetype == 'cfg':
         self.write_cfg(fh)
      elif filetype == 'cell':
         self.write_cell(fh)

      if opened: fh.close()


   def write_cfg(self, cfg=sys.stdout, shift=numpy.array([0.,0.,0.]), properties=None):
      """Write atoms in AtomEye extended CFG format. Returns a list of auxiliary properties
      actually written to CFG file, which may be abbreviated compared to those requested since
      AtomEye has a maximum of 32 aux props."""

      opened = False
      if type(cfg) == type(''):
         cfg = open(cfg, 'w')
         opened = True

      if properties is None:
         properties = self.properties.keys()

      # Header line
      cfg.write('Number of particles = %d\n' % self.n)
      cfg.write('# '+self.comment(properties)+'\n')

      # Lattice vectors
      for i in 0,1,2:
         for j in 0,1,2:
            cfg.write('H0(%d,%d) = %16.8f\n' % (i+1, j+1, self.lattice[i,j]))
   
      cfg.write('.NO_VELOCITY.\n')

      # Check first property is position-like
      species = getattr(self,properties[0])
      if len(species.shape) != 1 or species.dtype.kind != 'S':
         raise ValueError('First property must be species like')

      pos = getattr(self,properties[1])
      if pos.shape[1] != 3 or pos.dtype.kind != 'f':
         raise ValueError('Second property must be position like')

      if not self.properties.has_key('frac_pos'):
         self.add_property('frac_pos',0.0,ncols=3)
      self.frac_pos[:] = numpy.array([ numpy.dot(pos[i,:],self.g) + shift for i in range(self.n) ])

      if not self.properties.has_key('mass'):
         self.add_property('mass', map(ElementMass.get, self.species))

      properties = filter(lambda p: p not in ('pos','frac_pos','mass','species'), properties)

      # AtomEye can handle a maximum of 32 columns, so we might have to throw away
      # some of the less interesting propeeties
      
      def count_cols():
         n_aux = 0
         for p in properties:
            s = getattr(self,p).shape
            if len(s) == 1: n_aux += 1
            else:           n_aux += s[1]
         return n_aux

      boring_properties = ['travel','avgpos','oldpos','acc','velo']
      while count_cols() > 32:
         if len(boring_properties) == 0:
            raise ValueError('No boring properties left!')
         try:
            next_most_boring = boring_properties.pop(0)
            del properties[properties.index(next_most_boring)]
         except IndexError:
            pass # this boring property isn't in the list: move on to next

      properties = ['species','mass','frac_pos'] + properties
      data = self.to_recarray(properties)

      cfg.write('entry_count = %d\n' % (len(data.dtype.names)-2))

      # 3 lines per atom: element name, mass and other data
      format = '%s\n%12.4f\n'
      for i,name in enumerate(data.dtype.names[2:]):
         if i > 2: cfg.write('auxiliary[%d] = %s\n' % (i-3,name))
         format = format + _getfmt(data.dtype.fields[name][0])
      format = format + '\n'

      for i in range(self.n):
         cfg.write(format % tuple(data[i]))

      if opened: cfg.close()

      # Return column names as a list
      return list(data.dtype.names)


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


   def supercell(self, n1, n2, n3):

      other = Atoms(n=self.n*n1*n2*n3,n_int=self.int.shape[1],\
                    n_real=self.real.shape[1], \
                    properties=self.properties.copy())

      other.lattice[0,:] = self.lattice[0,:]*n1
      other.lattice[1,:] = self.lattice[1,:]*n2
      other.lattice[2,:] = self.lattice[2,:]*n3
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


   def cell_volume(self):
      return abs(numpy.dot(numpy.cross(self.lattice[0,:], self.lattice[1,:]),self.lattice[2,:]))


def diamond(species, a):
   "Bulk cube of element with lattice constant a"
   at = Atoms(n=8)

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


def norm2(a):
   "Nx1 array of the norm**2 of each vector in a Nx3 array"
   n, m = a.shape
   if m != 3:
      raise ValueError('second array dimension should be of size 3')
   out = numpy.zeros(n, dtype=float)
   for i in range(n):
      out[i] = numpy.dot(a[i,:],a[i,:])
   return out

def norm(a):
   "Return sqrt(norm2(a))"
   return numpy.sqrt(norm2(a))
                  

def count(condition):
   return numpy.where(condition)[0].size

def rms_diff(at1, at2):
   """Root mean square position difference between two Atoms objects.
    Raises ValueError if atom numbers don't match. """

   if (at1.n != at2.n):
      raise ValueError('at1.n (%d) != at2.n (%d)' % (at1.n, at2.n))
   
   return numpy.sqrt(sum(norm2(at1.pos - at2.pos))/at1.n)


def _getfmt(dtype):
   typ = dtype.type
   if issubclass(typ, numpy.bool_):
      return '%s'
   if issubclass(typ, numpy.integer):
      return '%8d'
   elif issubclass(typ, numpy.floating):
      return '%16.8f'
   elif issubclass(typ, numpy.complex):
      return '(%f,%f)'
   else:
      return '%s '


def _getconv(dtype):
    typ = dtype.type
    if issubclass(typ, numpy.bool_):
        return lambda x: bool(x)
    if issubclass(typ, numpy.integer):
        return int
    elif issubclass(typ, numpy.floating):
        return float
    elif issubclass(typ, numpy.complex):
        return complex
    else:
        return str

def _parse_properties(prop_str):
   """Parse a property description string in the format
   name:ptype:cols:...  and return a 5-tuple of an OrderedDict
   mapping property names to (ptype,cols) tuples and the number of
   int, real, string and logical properties"""

   properties = OrderedDict()
   if prop_str.startswith('pos:R:3'):
      prop_str = 'species:S:1:'+prop_str
      
   fields = prop_str.split(':')

   n_real = 0
   n_int  = 0
   n_str = 0
   n_logical = 0
   for name, ptype, cols in zip(fields[::3], fields[1::3], map(int,fields[2::3])):
      if ptype == 'R':
         start_col = n_real
         n_real = n_real + cols
         end_col = n_real
      elif ptype == 'I':
         start_col = n_int
         n_int = n_int + cols
         end_col = n_int
      elif ptype == 'S':
         start_col = n_str
         n_str = n_str + cols
         end_col = n_str
      elif ptype == 'L':
         start_col = n_logical
         n_logical = n_logical + cols
         end_col = n_logical
      else:
         raise ValueError('Unknown property type: '+ptype)

      properties[name] = (ptype,slice(start_col,end_col))

   return (properties, n_int, n_real, n_str, n_logical)



def frame_reader(fname,filetype=None):
   "Return an iterator which reads successive frames from a file"
   opened = False
   if type(fname) == type(''):
      if fname.endswith('.gz'):
         import gzip
         fh = gzip.open(fname)
         fname = fname[:-3] # remove .gz
      elif fname.endswith('.nc'):
         fh = fname
      else:
         fh = open(fname, 'r')
         opened = True

      # Guess file type from extension
      if filetype is None:
         root, filetype = os.path.splitext(fname)
         filetype = filetype[1:] # remove '.'

   else:
      fh = fname

   # Default to xyz format
   if filetype is None:
      filetype = 'xyz'

   if not filetype in ['cell','geom','xyz','castep','nc']:
      raise ValueError('Unrecognised filetype %s' % filetype)

   while fh:
      try:
         a = Atoms()
         if filetype == 'xyz':
            a.read_xyz(fh)
         elif filetype == 'cell':
            a.read_cell(fh)
         elif filetype == 'geom':
            a.read_geom(fh)
         elif filetype == 'castep':
            a.read_castep(fh)
         elif filetype == 'nc':
            a.read_netcdf(fh)
         yield a
      except IOError:
         if opened: fh.close()
         raise StopIteration

   if opened: fh.close()


def make_lattice(lengths, angles):
   "Form 3x3 lattice from cell lengths (a,b,c) and angles (alpha,beta,gamma)"

   a, b, c = lengths
   alpha, beta, gamma = angles

   lattice = numpy.zeros((3,3),'d')

   if alpha == 90.0 and beta == 90.0 and gamma == 90.0:
      lattice[0,0] = a
      lattice[1,1] = b
      lattice[2,2] = c
      return lattice

   alpha *= numpy.pi/180.0;  beta *= numpy.pi/180.0; gamma *= numpy.pi/180.0
   

   cos_alpha = numpy.cos(alpha); cos2_alpha = cos_alpha*cos_alpha
   cos_beta  = numpy.cos(beta);  cos2_beta  = cos_beta *cos_beta
   cos_gamma = numpy.cos(gamma); cos2_gamma = cos_gamma*cos_gamma
   sin_gamma = numpy.sin(gamma); sin2_gamma = sin_gamma*sin_gamma

   
   lattice[0,0] = a

   lattice[1,0] = b * cos_gamma
   lattice[1,1] = b * sin_gamma

   lattice[2,0] = c * cos_beta
   lattice[2,1] = c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma
   lattice[2,2] = c * numpy.sqrt(1.0 - (cos2_alpha + cos2_beta - 2.0*cos_alpha*cos_beta*cos_gamma)/ sin2_gamma)

   return lattice
