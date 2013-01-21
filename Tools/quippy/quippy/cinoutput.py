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
from quippy.atoms import Atoms
from quippy.io import AtomsReaders, AtomsWriters, atoms_reader
from quippy.farray import farray, padded_str_array
from quippy.extendable_str import Extendable_str
import numpy as np
import __builtin__
import itertools

__all__ = _cinoutput.__all__

class CInOutput(_cinoutput.CInOutput):

    __doc__ = _cinoutput.CInOutput.__doc__

    def __init__(self, filename=None, action=INPUT, append=False, netcdf4=True, no_compute_index=None,
                 frame=None, one_frame_per_file=None, mpi=None, zero=False, range=None, indices=None,
                 fpointer=None, finalise=True, string=False):

        self.string = None
        if string:
            self.string = filename
            filename = None
        _cinoutput.CInOutput.__init__(self, filename, action, append, netcdf4, no_compute_index,
                                      frame, one_frame_per_file, mpi, fpointer=fpointer, finalise=finalise)
        self.zero = zero
        self.range = range
        self.indices = indices

    __init__.__doc__ = _cinoutput.CInOutput.__init__.__doc__

    def __len__(self):
        if self.got_index:
            return int(self.n_frame)
        else:
            raise AttributeError('This CInOutput does not have an index')

    def __getitem__(self, index):
        return self.read(frame=index, zero=self.zero, range=self.range, indices=self.indices)

    def __setitem__(self, index, value):
        self.write(value, frame=index)

    def read(self, properties=None, properties_array=None, frame=None,
             zero=None, range=None, str=None, estr=None, indices=None):
        at = Atoms()
        if range == 0 or range == 'empty':
            range = [-1,-1]

        if self.string is not None and str is None:
            str = self.string

        try:
            _cinoutput.CInOutput.read(self, at, properties=properties,
                                      properties_array=properties_array, frame=frame,
                                      zero=zero, range=range, str=str, estr=estr,
                                      indices=indices)
        except RuntimeError, re:
            if 'kind IO EOF' in __builtin__.str(re):
                raise EOFError(__builtin__.str(re))
            else:
                raise
            
        return at

    def write(self, at, properties=None, prefix=None, int_format=None, real_format=None, frame=None,
              shuffle=None, deflate=None, deflate_level=None, estr=None, update_index=None):

        if properties is not None and (not hasattr(properties, 'dtype') or properties.dtype != dtype('S1')):
            properties = padded_str_array(properties, max([len(x) for x in properties])).T

        if self.string and estr is None:
            estr = Extendable_str()
            
        _cinoutput.CInOutput.write(self, at, properties_array=properties, prefix=prefix, int_format=int_format,
                                   real_format=real_format, frame=frame, shuffle=shuffle, deflate=deflate,
                                   deflate_level=deflate_level, estr=estr, update_index=update_index)
        if estr is not None:
            return str(estr)
    

from quippy import FortranDerivedTypes
FortranDerivedTypes['type(cinoutput)'] = CInOutput

class CInOutputReader(object):
    """Class to read atoms from a CInOutput. Supports generator and random access via indexing."""

    def __init__(self, source, frame=None, range=None, start=0, stop=None, step=1, no_compute_index=False,
                 zero=False, one_frame_per_file=False, indices=None, string=False):
        if isinstance(source, basestring):
            self.opened = True
            self.source = CInOutput(source, action=INPUT, append=False, zero=zero, range=range,
                                    no_compute_index=no_compute_index, one_frame_per_file=one_frame_per_file,
                                    indices=indices, string=string)
            self.netcdf_file = None
            if self.source.string is None:
                try:
                    self.netcdf_file = netcdf_file(source)
                except (RuntimeError, AssertionError, IOError):
                    pass
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
                try:
                    stop = len(self.source)
                except AttributeError:
                    stop = None
            self.stop = stop
            self.step = step

    def __iter__(self):
        frames = itertools.islice(itertools.count(0), self.start,self.stop,self.step)
        for frame in frames:
            try:
                yield self.source[frame]
            except EOFError:
                break

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

    def __init__(self, dest, append=False, netcdf4=True, one_frame_per_file=False, string=False, **write_kwargs):
        self.opened = False
        self.write_kwargs = {}
        self.write_kwargs.update(write_kwargs)
        if isinstance(dest, basestring):
            self.opened = True
            self.dest = CInOutput(dest, action=OUTPUT, append=append, netcdf4=netcdf4,
                                  one_frame_per_file=one_frame_per_file, string=string)
        else:
            self.dest = dest

    def write(self, at, **kwargs):
        kwargs.update(self.write_kwargs)
        return self.dest.write(at, **kwargs)

    def close(self):
        self.dest.close()

AtomsWriters['xyz'] = AtomsWriters['nc'] = AtomsWriters[CInOutput] = CInOutputWriter

    

class CInOutputStringReader(CInOutputReader):
    def __init__(self, source, frame=None, range=None, start=0, stop=None, step=1, no_compute_index=False,
                 zero=False, one_frame_per_file=False, indices=None):
        CInOutputReader.__init__(self, source=source, frame=frame, range=range, start=start, stop=stop, step=step,
                                 no_compute_index=no_compute_index, zero=zero, one_frame_per_file=one_frame_per_file,
                                 indices=indices, string=True)

AtomsReaders['string'] = CInOutputStringReader

class CInOutputStringWriter(CInOutputWriter):
    def __init__(self, dest, append=False, netcdf4=True, one_frame_per_file=False, **write_kwargs):
        assert dest == 'string'
        CInOutputWriter.__init__(self, dest=dest, append=append, netcdf4=netcdf4, one_frame_per_file=one_frame_per_file,
                                 string=True, **write_kwargs)

AtomsWriters['string'] = CInOutputStringWriter


@atoms_reader('out')
def EvalOutputReader(source):
    """
    Reader for QUIP `eval` output files

    Extracts XYZ configurations from lines beginnging with "AT "
    """
    
    lines = open(source, 'r').readlines()
    lines = [line[3:] for line in lines if line.startswith('AT ')]
    reader = CInOutputStringReader(''.join(lines))
    for at in reader:
        yield at
