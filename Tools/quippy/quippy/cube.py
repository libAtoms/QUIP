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

from quippy.atoms import Atoms
from quippy.io import AtomsReaders, AtomsWriters
from quippy.units import BOHR
from quippy.farray import *
import sys
import numpy as np

__all__ = ['CubeWriter', 'CubeReader']

class CubeWriter(object):

    def __init__(self, f, comment=None, data=None, origin=None, extent=None, comment2=None):
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
        self.extent = extent
        self.comment2 = comment2

    def write(self, at):

        comment  = self.comment
        data     = self.data
        origin   = self.origin
        extent   = self.extent
        comment2 = self.comment2
        
        if data is None and hasattr(at, 'data'):
            data = at.data
        if data is None:
            raise ValueError("Cannot write .cube file without any volumetric data")

        if comment is None and 'comment' in at.params:
            comment = at.params['comment1']
        if comment is None: comment = ''

        comment2 = comment2
        if comment2 is None and 'comment2' in at.params:
            comment2 = at.params['comment2']
        if comment2 is None: comment2 = ''

        origin = origin
        if origin is None and 'origin' in at.params:
            origin = at.params['origin']
        if origin is None: origin = [0., 0., 0.]
        origin = farray(origin)

        self.f.write(comment.strip()+'\n')
        self.f.write(comment2.strip()+'\n')

        origin_x, origin_y, origin_z = origin/BOHR
        self.f.write('%d %f %f %f\n' % (at.n, origin_x, origin_y, origin_z))

        if extent is None and 'extent' in at.params:
            extent = at.params['extent']
        if extent is None:
            extent = at.lattice
        if extent.shape == (3,):
            extent = np.diag(extent)
        extent = farray(extent)

        for i in (1,2,3):
            n = data.shape[i-1]
            voxel_x, voxel_y, voxel_z = extent[:,i]/BOHR/n
            n += 1 # for padding
            self.f.write('%d %f %f %f\n' % (n, voxel_x, voxel_y, voxel_z))

        for i in frange(at.n):
            self.f.write('%d 0.0 %f %f %f\n' % (at.z[i], at.pos[1,i]/BOHR, at.pos[2,i]/BOHR, at.pos[3,i]/BOHR))

        padded_data = np.zeros([s+1 for s in data.shape])
        padded_data[:-1,:-1,:-1] = data
        padded_data[-1,:,:] = padded_data[0,:,:]
        padded_data[:,-1,:] = padded_data[:,0,:]
        padded_data[:,:,-1] = padded_data[:,:,0]
        for d in padded_data.flat:
            self.f.write('%f\n' % d)

    def close(self):
        if self.opened:  self.f.close()


def CubeReader(f, property_name='charge', discard_repeat=True):

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
    origin = farray([origin_x, origin_y, origin_z])*BOHR

    # Next three lines define number of voxels and shape of each element
    shape = [0,0,0]
    voxel = fzeros((3,3))
    for i in (1,2,3):
        shape[i-1], voxel[1,i], voxel[2,i], voxel[3,i] = convert_line(f.readline(), int, float, float, float)

    at = Atoms(n=n_atom, lattice=voxel*BOHR*shape)
    at.add_property(property_name, 0.0)
    prop_array = getattr(at, property_name)

    # Now there's one line per atom
    for i in frange(at.n):
        at.z[i], prop_array[i], at.pos[1,i], at.pos[2,i], at.pos[3,i] = convert_line(f.readline(), int, float, float, float, float)
        at.pos[:,i] *= BOHR
    at.set_atoms(at.z)

    # Rest of file is volumetric data
    data = np.fromiter((float(x) for x in f.read().split()),float,count=-1)
    if data.size != shape[0]*shape[1]*shape[2]:
        raise IOError("Bad array length - expected shape %r, but got size %d" % (shape, data.size))

    # Save volumetric data in at.data
    data = farray(data.reshape(shape))

    # Discard periodic repeats?
    if discard_repeat:
        at.data = data[:-1,:-1,:-1]
        shape = [s-1 for s in shape]
        at.set_lattice(voxel*BOHR*shape, False)

    at.params['comment1'] = comment1
    at.params['comment2'] = comment2
    at.params['origin'] = origin
    at.params['shape'] = shape

    # save grids in at.grid_x, at.grid_y, at.grid_z
    if at.is_orthorhombic:
        at.grid_x, at.grid_y, at.grid_z = np.mgrid[origin[1]:origin[1]+at.lattice[1,1]:shape[0]*1j,
                                                   origin[2]:origin[2]+at.lattice[2,2]:shape[1]*1j,
                                                   origin[3]:origin[3]+at.lattice[3,3]:shape[2]*1j]

    if opened:
        f.close()

    yield at



AtomsReaders['cube'] = CubeReader
AtomsWriters['cube'] = CubeWriter
