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

from quippy import Atoms, AtomsReaders, AtomsWriters, BOHR
from farray import *
import sys

class CubeWriter(object):

   def __init__(self, f, comment=None, data=None, origin=None, comment2=None):
      if type(f) == type(''):
         if f == 'stdout':
            self.f = sys.stdout
            self.opened = False
         else:
            self.opened = True
            self.f = open(f,'w')
      else:
         self.opened = False
         self.f = f

      self.comment1 = comment
      self.comment2 = comment2
      self.data = data
      self.origin = origin

   def write(self, at):
       data = self.data
       if data is None and hasattr(at, 'data'):
          data = at.data
       if data is None:
          raise ValueError("Cannot write .cube file without any volumetric data")

       comment1 = self.comment1
       if comment1 is None and 'comment1' in at.params:
          comment1 = at.params['comment1']
       if comment1 is None: comment1 = ''

       comment2 = self.comment2
       if comment2 is None and 'comment2' in at.params:
          comment2 = at.params['comment2']
       if comment2 is None: comment2 = ''

       origin = self.origin
       if self.origin is None and 'origin' in at.params:
          origin = at.params['origin']
       if origin is None: origin = (0., 0., 0.)

       self.f.write(comment1.strip()+'\n')
       self.f.write(comment2.strip()+'\n')
       
       origin_x, origin_y, origin_z = origin
       self.f.write('%d %f %f %f\n' % (at.n, origin_x, origin_y, origin_z))

       for i in (1,2,3):
           n = data.shape[i-1]
           voxel_x, voxel_y, voxel_z = at.lattice[:,i]/BOHR/n
           self.f.write('%d %f %f %f\n' % (n, voxel_x, voxel_y, voxel_z))

       for i in frange(at.n):
           self.f.write('%d 0.0 %f %f %f\n' % (at.z[i], at.pos[1,i]/BOHR, at.pos[2,i]/BOHR, at.pos[3,i]/BOHR))

       for d in data.flat:
           self.f.write('%f\n' % d)

   def close(self):
       if self.opened:  self.f.close()


def CubeReader(f):

   def convert_line(line, *fmts):
      return (f(s) for f,s in zip(fmts, line.split()))
   
   if type(f) == type(''):
      f = open(f)
      opened = True

   # First two lines are comments
   comment1 = f.readline()
   comment2 = f.readline()

   # Now number of atoms and origin
   n_atom, origin_x, origin_y, origin_z = convert_line(f.readline(), int, float, float, float)

   # Next three lines define number of voxels and shape of each element
   shape = [0,0,0]
   voxel = fzeros((3,3))
   for i in (1,2,3):
      shape[i-1], voxel[1,i], voxel[2,i], voxel[3,i] = convert_line(f.readline(), int, float, float, float)

   at = Atoms(n=n_atom, lattice=voxel*BOHR*shape)
   at.params['comment1'] = comment1
   at.params['comment2'] = comment2
   at.params['origin'] = (origin_x, origin_y, origin_z)

   # Now there's one line per atom
   for i in frange(at.n):
      at.z[i], unknown, at.pos[1,i], at.pos[2,i], at.pos[3,i] = convert_line(f.readline(), int, float, float, float, float)
      at.pos[:,i] *= BOHR
   at.set_atoms(at.z)

   # Rest of file is data, with one data point per line
   data = numpy.loadtxt(f)
   if data.size != shape[0]*shape[1]*shape[2]:
      raise IOError("Bad array length - expected shape %r, but got size %d" % (shape, data.size))

   # Save volumetric data in at.data
   at.data = farray(data.reshape(shape))

   if opened:
      f.close()

   yield at
   


AtomsReaders['cube'] = CubeReader
AtomsWriters['cube'] = CubeWriter
