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

from quippy import AtomsWriters, BOHR
from farray import *
import sys

class CubeWriter(object):

   def __init__(self, f, comment='', data=None, origin=(0., 0., 0.)):
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

      self.comment = comment
      self.data = data
      self.origin = origin

   def write(self, at):
       self.f.write(self.comment.strip()+'\n\n')
       
       origin_x, origin_y, origin_z = self.origin
       self.f.write('%d %f %f %f\n' % (at.n, origin_x, origin_y, origin_z))

       for i in (1,2,3):
           n = self.data.shape[i-1]
           voxel_x, voxel_y, voxel_z = at.lattice[:,i]/BOHR/n
           self.f.write('%d %f %f %f\n' % (n, voxel_x, voxel_y, voxel_z))

       for i in frange(at.n):
           self.f.write('%d 0.0 %f %f %f\n' % (at.z[i], at.pos[1,i]/BOHR, at.pos[2,i]/BOHR, at.pos[3,i]/BOHR))

       for d in self.data.flat:
           self.f.write('%f\n' % d)

       

   def close(self):
       if self.opened:  self.f.close()



AtomsWriters['cube'] = CubeWriter
