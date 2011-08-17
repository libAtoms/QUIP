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
import sys, numpy

class IMDReader:

    def __init__(self, filename, z=14):
        self.filename = filename
        self.z = z

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
        data = numpy.loadtxt(f)

        lattice = fzeros((3,3))
        for i in [1,2,3]:
            lattice[:,i] = [float(x) for x in header[i+1].split()[1:]]
        at = Atoms(n=len(data), lattice=lattice)
        at.pos[...] = data[:,3:6].T
        at.set_atoms(self.z)
        at.add_property('tag', data[:,1])
        at.add_property('mass', data[:,2])
        at.add_property('velo', data[:,6:9].T)
        at.add_property('epot', data[:,9])
        if opened:
            f.close()

        yield at



AtomsReaders['chkpt'] = IMDReader
