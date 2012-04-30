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

from quippy import _cinoutput
from quippy._cinoutput import *

from quippy import netcdf_file
from quippy.system import INPUT, OUTPUT, INOUT
from quippy.atoms import Atoms, AtomsReaders, AtomsWriters, atoms_reader
from quippy.farray import farray, padded_str_array
from quippy.extendable_str import Extendable_str
import numpy as np

__all__ = _cinoutput.__all__

class CInOutput(_cinoutput.CInOutput):

    __doc__ = _cinoutput.CInOutput.__doc__

    def __init__(self, filename=None, action=INPUT, append=False, netcdf4=True, no_compute_index=None,
                 frame=None, one_frame_per_file=None, mpi=None, zero=False, range=None, indices=None,
                 fpointer=None, finalise=True):
        _cinoutput.CInOutput.__init__(self, filename, action, append, netcdf4, no_compute_index,
                                      frame, one_frame_per_file, mpi, fpointer=fpointer, finalise=finalise)
        self.zero = zero
        self.range = range
        self.indices = indices

    __init__.__doc__ = _cinoutput.CInOutput.__init__.__doc__

    def __len__(self):
        return int(self.n_frame)

    def __getitem__(self, index):
        return self.read(frame=index, zero=self.zero, range=self.range, indices=self.indices)

    def __setitem__(self, index, value):
        self.write(value, frame=index)

    def read(self, properties=None, properties_array=None, frame=None,
             zero=None, range=None, str=None, estr=None, indices=None):
        at = Atoms()
        if range == 0 or range == 'empty':
            range = [-1,-1]

        _cinoutput.CInOutput.read(self, at, properties=properties,
                                  properties_array=properties_array, frame=frame,
                                  zero=zero, range=range, str=str, estr=estr,
                                  indices=indices)
        return at

    def write(self, at, properties=None, prefix=None, int_format=None, real_format=None, frame=None,
              shuffle=None, deflate=None, deflate_level=None, estr=None, update_index=None):

        if properties is not None and (not hasattr(properties, 'dtype') or properties.dtype != dtype('S1')):
            properties = padded_str_array(properties, max([len(x) for x in properties])).T

        _cinoutput.CInOutput.write(self, at, properties_array=properties, prefix=prefix, int_format=int_format,
                                   real_format=real_format, frame=frame, shuffle=shuffle, deflate=deflate,
                                   deflate_level=deflate_level, estr=estr, update_index=update_index)

from quippy import FortranDerivedTypes
FortranDerivedTypes['type(cinoutput)'] = CInOutput

class CInOutputReader(object):
    """Class to read atoms from a CInOutput. Supports generator and random access via indexing."""

    def __init__(self, source, frame=None, range=None, start=0, stop=None, step=1, no_compute_index=False,
                 zero=False, one_frame_per_file=False, indices=None):
        if isinstance(source, basestring):
            self.opened = True
            self.source = CInOutput(source, action=INPUT, append=False, zero=zero, range=range,
                                    no_compute_index=no_compute_index, one_frame_per_file=one_frame_per_file,
                                    indices=indices)
            try:
                self.netcdf_file = netcdf_file(source)
            except (RuntimeError, AssertionError, IOError):
                self.netcdf_file = None
        else:
            self.opened = False
            self.source = source
            self.netcdf_file = None

        if frame is not None:
            self.start = frame
            self.stop = frame+1
            self.step = 1
        else:
            self.start = start
            if stop is None:
                stop = len(self.source)
            self.stop = stop
            self.step = step

    def __iter__(self):
        for frame in range(self.start,self.stop,self.step):
            yield self.source[frame]

    def close(self):
        if self.opened: self.source.close()

    def __len__(self):
        return len(self.source)

    def __getitem__(self, idx):
        return self.source[idx]

    def __getattr__(self, name):
        if self.netcdf_file is not None:
            try:
                return self.netcdf_file.__getattr__(name)
            except AttributeError:
                try:
                    return farray(self.netcdf_file.variables[name][:])
                except KeyError:
                    raise AttributeError('Attribute %s not found' % name)
        else:
            raise AttributeError('Attribute %s not found' % name)


AtomsReaders['xyz'] = AtomsReaders['nc'] = AtomsReaders[CInOutput] = CInOutputReader

@atoms_reader('stdin')
def CInOutputStdinReader(source='stdin'):
    assert source == 'stdin'
    source = CInOutput(source, action=INPUT)
    while True:
        try:
            yield source.read()
        except RuntimeError:
            break

class CInOutputWriter(object):
    """Class to write atoms sequentially to a CInOutput stream"""

    def __init__(self, dest, append=False, netcdf4=True, one_frame_per_file=False, **write_kwargs):
        self.opened = False
        self.write_kwargs = {}
        self.write_kwargs.update(write_kwargs)
        if isinstance(dest, basestring):
            self.opened = True
            self.dest = CInOutput(dest, action=OUTPUT, append=append, netcdf4=netcdf4,
                                  one_frame_per_file=one_frame_per_file)
        else:
            self.dest = dest

    def write(self, at, **kwargs):
        kwargs.update(self.write_kwargs)
        self.dest.write(at, **kwargs)

    def close(self):
        self.dest.close()

AtomsWriters['xyz'] = AtomsWriters['nc'] = AtomsWriters[CInOutput] = CInOutputWriter

    
