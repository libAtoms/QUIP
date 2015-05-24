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
from farray import *
import sys
import numpy as np

__all__ = ['IMDReader']

class IMDReader:

    def __init__(self, filename, z=14, fix_tags=None, vacuum=None):
        self.filename = filename
        self.z = z
        self.fix_tags = fix_tags
        self.vacuum = vacuum

    def __len__(self):
        return 1

    def __getitem__(self, index):
        if index == 0:
            return iter(self).next()
        else:
            raise IndexError

    def __iter__(self):
        if type(self.filename) == type(''):
            f = open(self.filename)
            opened = True

        # First eight lines are header
        header = [f.readline() for i in range(8)]

        # Rest of file is data array
        data = np.loadtxt(f)

        lattice = fzeros((3,3))
        for i in [1,2,3]:
            lattice[:,i] = [float(x) for x in header[i+1].split()[1:]]
        if self.vacuum is not None:
            lattice += np.diag(self.vacuum)
        at = Atoms(n=len(data), lattice=lattice)
        at.pos[...] = data[:,3:6].T
        at.set_atoms(self.z)
        at.add_property('tag', data[:,1].astype(int))
        at.add_property('mass', data[:,2])
        at.add_property('velo', data[:,6:9].T)
        at.add_property('epot', data[:,9])
        if self.fix_tags is not None:
            at.add_property('move_mask', 1)
            for tag in self.fix_tags:
                at.move_mask[at.tag == tag] = 0
        if opened:
            f.close()

        yield at



AtomsReaders['chkpt'] = IMDReader
