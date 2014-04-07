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

"""
There are two classes for reading trajectories: :class:`AtomsReader`
and :class:`AtomsList`. Use an :class:`AtomsReader` for quick
read-only access to a trajectory or if you only want to access some of
the frames. If you want to load the entire file into memory and
manipulate it use an :class:`AtomsList`.
"""

import sys, os, fnmatch, re, itertools, glob, operator, warnings, math, logging

import numpy as np

from quippy.system import mem_info
from quippy.util import infer_format, parse_slice, time_ordered_glob
from quippy.mockndarray import mockNDarray
from quippy.atoms import Atoms
from quippy.farray import FortranArray

__all__ = ['AtomsReaders', 'AtomsWriters', 'atoms_reader',
           'AtomsReader', 'AtomsWriter', 'AtomsList',
           'read_dataset', 'time_ordered_series', 'read', 'write', 'dict2atoms']

AtomsReaders = {}
AtomsWriters = {}

major, minor = sys.version_info[0:2]
assert (major, minor) >= (2, 4)
if (major, minor) < (2, 5):
    all = lambda seq: not False in seq
    any = lambda seq: True in seq
    __all__.extend(['all', 'any'])
del major, minor

def atoms_reader(source):
    """Decorator to mark a function as a reader for a particular file extension"""
    def decorate(func):
        global AtomsReaders
        if not source in AtomsReaders:
            AtomsReaders[source] = func
        return func
    return decorate

class AtomsReaderMixin(object):
    def __repr__(self):
        try:
            len_self = len(self)
        except:
            len_self = None
        return '<%s source=%r format=%r start=%r stop=%r step=%r random_access=%r len=%r>' % \
               (self.__class__.__name__, self.source, self.format, self._start, self._stop, self._step,
                self.random_access, len_self)

    def write(self, dest=None, format=None, properties=None, prefix=None,
              progress=False, progress_width=80, update_interval=None,
              show_value=True, **kwargs):
        """
        Write all frames to `dest`. If `format` is not
        given it is inferred from the file extension of `dest` (see
        :ref:`fileformats`). If `properties` is present, it should be a list
        of property names to include in the output file, e.g. `['species', 'pos']`.
      
        `progress`, `progress_width`, `update_interval` and `show_value`
        are used to control a textual progress bar. The extra arguments
        in `*args` and `**kwargs` are passed along to the underlying
        writer routine constructed for writing to `dest`.

        See :ref:`fileformats` for a list of supported file formats.
        """
        opened = False
        if dest is None:
            dest = self.source
        filename, dest, format = infer_format(dest, format, AtomsWriters)

        if progress:
            from progbar import ProgressBar
            pb = ProgressBar(0,len(self),progress_width,showValue=show_value)
            update_interval = update_interval or max(1, len(self)/progress_width)

        if format in AtomsWriters:
            dest = AtomsWriters[format](dest, **kwargs)

        if not hasattr(dest, 'write'):
            raise ValueError("Don't know how to write to destination \"%s\" in format \"%s\"" % (dest, format))

        res = []
        for i, a in enumerate(self):
            write_kwargs = {}
            if properties is not None: write_kwargs['properties'] = properties
            if prefix is not None: write_kwargs['prefix'] = prefix
            try:
                res.append(dest.write(a, **write_kwargs))
            except TypeError:
                raise ValueError('destination does not support specifying arguments %r' % write_kwargs)

            if progress and i % update_interval == 0: pb(i)

        if opened:
            dest.close()

        # Special case for writing to a string
        if format == 'string':
            return ''.join(res)
        else:
            if res is not None and not all(el is None for el in res):
                return res

class AtomsReader(AtomsReaderMixin):
    """
    An :class:`AtomsReader` reads a series of :class:`Atoms` objects
    from the trajectory `source` which should be one of the following:
    
     * a filename - in this case `format` is inferred from the file
       extension -- see :ref:`fileformats`
     * a shell-style glob pattern e.g. `"*.xyz"`
     * a list of filenames or glob patterns e.g. `["foo*.xyz",
       "bar*.xyz"]`
     * an open file or file-like object (e.g. a :class:`CInOutput`
       object)
     * any Python `iterator
       <http://docs.python.org/library/stdtypes.html#iterator-types>`_
       which yields a sequence of :class:`Atoms` objects
    
    `start`, `stop` and `step` can be used to restrict the range of frames
    read from `source`. The first frame in the file has index zero.
    
    `cache_limit` determines how many configurations will be stored in
    memory. If more than `cache_limit` configurations are read in, the
    least recently accessed configurations are thrown away. To store
    everything, use an :class:`AtomsList` instead.
    
    Some `sources` understand additional keyword arguments from
    `**kwargs`. For example the CASTEP file reader can take an
    `atoms_ref` argument which is a reference :class:`Atoms` object
    which is used to fill in information which is missing from the
    input file.
    
    All :class:`AtomsReaders` support iteration, so you can loop over
    the contents using a :keyword:`for` loop::
    
       al = AtomsReader('input-file.xyz')
       for at in al:
          # process Atoms object `at`
          print at.energy
    
    or using list comprehension::
    
       print [at.energy for at in al]
    
    In addition to iteration, some sources allow random access. To find
    out if an :class:`AtomsReader` supports random access, either try
    to get it's length with :func:`len`, or check if the
    :attr:`random_access` property is true. If `cache_limit` is large
    enough to store all the frames in the file, all
    :class:`AtomsReaders` will allow random access once the entire
    trajectory has been loaded.
    
    If :attr:`randomaccess` is true, you can access individual frames
    by indexing and slicing, e.g. ``al[i]`` is the i\ :sup:`th`
    :class:`Atoms` object within ``al`` and ``al[i:j]`` returns objects
    from `i` upto but not including `j`. Like ordinary Python lists,
    indices start from 0 and run up to ``len(al)-1``.
    """

    def __init__(self, source, format=None, start=None, stop=None, step=None,
                 cache_mem_limit=-1, rename=None, **kwargs):

        def file_exists(f):
            return f == "stdin" or os.path.exists(f) or len(glob.glob(f)) > 0

        self.source = source
        self.format = format
        self._start = start
        self._stop = stop
        self._step = step

        self.cache_mem_limit = cache_mem_limit
        logging.debug('AtomsReader memory limit %r' % self.cache_mem_limit)

        self._source_len = None
        self._cache_dict = {}
        self._cache_list  = []
        self._cache_mem_usage = []

        self.opened = False
        self.reader = source

        self.rename = rename

        if isinstance(self.reader, basestring):
            if '@' in self.reader:
                self.reader, frames = self.reader.split('@')
                if self.reader.endswith('.db'):
                    self.reader = self.reader+'@'+frames
                    format = 'db'
                else:
                    frames = parse_slice(frames)
                    if start is not None or stop is not None or step is not None:
                        raise ValueError('Conflicting frame references start=%r stop=%r step=%r and @-syntax %r' %
                                         (start, stop, step, frames))
                    if isinstance(frames, int):
                        if frames >= 0:
                            frames = slice(frames, frames+1,+1)
                        else:
                            frames = slice(frames, frames-1,-1)

                    self._start, self._stop, self._step = frames.start, frames.stop, frames.step
                
            self.filename = self.reader
            self.opened = True
            if self.reader in AtomsReaders:
                if format is None:
                    format = self.reader
            elif format != 'string' and format != 'db':
                self.reader = os.path.expanduser(self.reader)
                glob_list = sorted(glob.glob(self.reader))
                if len(glob_list) == 0:
                    raise IOError("input file '%s' not found" % self.reader)
                if len(glob_list) > 1:
                    self.reader = glob_list
                else:
                    self.reader = glob_list[0]
                    filename, self.reader, new_format = infer_format(self.reader, format, AtomsReaders)
                    
                    if format is None:
                        format = new_format

        # special cases if source is a list or tuple of filenames or Atoms objects
        is_filename_sequence = False
        is_list_of_atoms = False
        if isinstance(self.reader, list) or isinstance(self.reader, tuple):
            is_filename_sequence = True
            is_list_of_atoms = True
            for item in self.reader:
                if '@' in item:
                    item = item[:item.index('@')]
                if not isinstance(item, basestring) or not file_exists(item):
                    is_filename_sequence = False
                if not isinstance(item, Atoms):
                    is_list_of_atoms = False

        if is_filename_sequence:
            self.reader = AtomsSequenceReader(self.reader, format=format, **kwargs)
        elif is_list_of_atoms:
            # dummy reader which copies from an existing list or tuple of Atoms objects
            self.reader = [ at.copy() for at in self.reader ]
        else:
            if format is None:
                format = self.reader.__class__
            if format in AtomsReaders:
                self.reader = AtomsReaders[format](self.reader, **kwargs)

        # check if reader is still a string or list of strings - indicates missing files or unknown format
        if isinstance(self.reader, basestring):
            raise IOError("Cannot read Atoms from file %s" % self.reader)
        elif isinstance(self.reader, list):
            is_list_of_strings = True
            for item in self.reader:
                if not isinstance(item, basestring):
                    is_list_of_strings = False
                    break
            if is_list_of_strings:
                raise IOError("Cannot read Atoms from files %s" % self.reader)

        if isinstance(self.reader, AtomsReader):
            self.reader = AtomsReaderCopier(self.reader)

        if not hasattr(self.reader, '__iter__'):
            # call Atoms constructor - this has beneficial side effect of making a copy
            self.reader = [Atoms(self.reader)]

    def __len__(self):
        if self._source_len is not None:
            return len(range(*slice(self._start, self._stop, self._step).indices(self._source_len)))
        elif hasattr(self.reader, '__len__') and hasattr(self.reader, '__getitem__'):
            try:
                return len(range(*slice(self._start, self._stop, self._step).indices(len(self.reader))))
            except:
                raise AttributeError('This AtomsReader does not support random access')
        else:
            raise AttributeError('This AtomsReader does not support random access')

    @property
    def random_access(self):
        """
        Read only property: True if this source supports random access, False if it does not
        """
        try:
            len(self)
            return True
        except:
            return False

    def close(self):
        """
        Close any open files associated with this :class:`AtomsReader`
        """
        if self.opened and hasattr(self.reader, 'close'):
            self.reader.close()

    def __getslice__(self, first, last):
        return self.__getitem__(slice(first,last,None))

    def _cache_fetch(self, frame):
        # least recently used (LRU) cache
        try:
            self._cache_list.append(self._cache_list.pop(self._cache_list.index(frame)))
            return True   # cache hit
        except ValueError:
            return False  # cache miss

    def _cache_store(self, frame, at):
        self._cache_list.append(frame)
        self._cache_dict[frame] = at

        if self.cache_mem_limit is not None:
            if self.cache_mem_limit == -1:
                self.cache_mem_limit = min(10*at.mem_estimate(), 100*1024**2)
            
            if self.cache_mem_limit == 0:
                while len(self._cache_dict) > 1:
                    logging.debug('Reducing AtomsReader cache size from %d' % len(self._cache_dict))
                    del self._cache_dict[self._cache_list.pop(0)]
            else:
                self._cache_mem_usage.append(at.mem_estimate())
                while len(self._cache_dict) > 1 and sum(self._cache_mem_usage) > self.cache_mem_limit:
                    logging.debug('Reducing AtomsReader cache size from %d' % len(self._cache_dict))
                    self._cache_mem_usage.pop(0)
                    del self._cache_dict[self._cache_list.pop(0)]

    def __getitem__(self, frame):
        if not self.random_access:
            raise IndexError('This AtomsReader does not support random access')

        if isinstance(frame, int) or isinstance(frame, np.integer):
            source_len = self._source_len or len(self.reader)
            if self._start is not None or self._stop is not None or self._step is not None:
                frame = range(*slice(self._start, self._stop, self._step).indices(source_len))[frame]
            if frame < 0: frame = frame + len(self)

            if not self._cache_fetch(frame):
                at = self.reader[frame]
                at = self.filter(at)
                self._cache_store(frame, at)

            at = self._cache_dict[frame]
            if not hasattr(at, 'source'):
                at.source = self.source
            if not hasattr(at, 'frame'):
                at.frame = frame
            return at

        elif isinstance(frame, slice):
            return self.__class__([self[f] for f in range(*frame.indices(len(self))) ])
        else:
            raise TypeError('frame should be either an integer or a slice')

    def __setitem__(self, frame, at):
        self._cache_store(frame, at)

    def iterframes(self, reverse=False):
        """
        Return an interator over all the frames in this trajectory. This
        is the default iterator for an :class:`AtomsReader` instance
        `al`, and can be accessed with ``iter(al)``. 

        If `reverse=True` then the iteration starts with the last frame
        and goes backwards through the file. This is only possible if
        :attr:`random_access` is true.
        """
        if self.random_access:
            # iterate using __getitem__, which automatically goes through LRU cache

            frames = range(len(self))
            if reverse: frames = reversed(frames)
            for f in frames:
                yield self[f]

        else:
            # source does not support random access
            # for each call, we make a new iterator from self.reader

            if reverse:
                raise IndexError('Cannot reverse iterate over an AtomsReader which does not support random access')

            frames = itertools.count()
            atoms = iter(self.reader)

            if self._start is not None or self._stop is not None or self._step is not None:
                if (self._start is not None and self._start < 0) or (self._stop is not None and self._stop < 0):
                    raise IndexError('Cannot use negative start or stop indices for an AtomsReader which does not support random access')
                
                frames = itertools.islice(frames, self._start or 0, self._stop or None, self._step or 1)
                atoms  = itertools.islice(atoms, self._start or 0, self._stop or None, self._step or 1)

            n_frames = 0
            last_frame = 0
            for (frame,at) in itertools.izip(frames, atoms):
                self._cache_store(frame, at)
                n_frames += 1
                last_frame = frame
                if not hasattr(at, 'source'):
                    at.source = self.source
                if not hasattr(at, 'frame'):
                    at.frame = frame
                yield at

            # once iteration is finished, random access will be possible if all frames fitted inside cache
            if len(self._cache_dict) == n_frames:
                self.reader = self._cache_dict
                self._source_len = last_frame+1

    def __iter__(self):
        return self.iterframes()

    def __reversed__(self):
        return self.iterframes(reverse=True)

    def filter(self, at):
        """
        Apply read-time filters to `at`
        """
        if self.rename is not None:
            for (old, new) in self.rename:
                if old in at.properties:
                    at.properties[new] = at.properties[old]
                    del at.properties[old]
                elif old in at.params:
                    at.params[new] = at.params[old]
                    del at.params[old]
                else:
                    raise AttributeError('Cannont rename: no property or parameter named "%s" exists' % old)
        return at


class AtomsList(AtomsReaderMixin, list):
    """
    An :class:`AtomsList` is just like an :class:`AtomsReader` except
    that all frames are read in on initialiased and then stored in
    memory. This is equivalent to an :class:`AtomsReader` with a
    `cache_limit` of `None` so an :class:`AtomsList` always
    supports random access.  
    
    The :class:`AtomsList` allows configurations to be added, removed
    or reordered using the standard Python methods for `mutable
    sequence types
    <http://docs.python.org/library/stdtypes.html#mutable-sequence-types>`_
    (e.g. :meth:`append`, :meth:`extend`, :meth:`index`, etc).
    
    The attributes of the component :class:`Atoms` can be accessed as a
    single array, using the frame number as the first array index. Note
    that the first index runs from 0 to `len(al)-1`, unlike the other
    indices which are one-based since the :class:`Atoms` attributes are
    stored in a :class:`FortranArray`.
    
    For example the following statements are all true::
    
       al.energy      ==  [at.energy for at in al] # energies of all atoms
       al.energy[0]   ==  al[0].energy             # energy of first frame
       all(al.velo[0] ==  al[0].velo)              # velocities of all atoms in first frame
       al.velo[0,-1]  ==  al[0].velo[-1]           # velocity of last atom in first frame
    
    In addition to the standard Python list methods and those of
    :class:`AtomsReader`, :class:`AtomsList` defined a couple of extras
    methods.
    """
    
    def __init__(self, source=[], format=None, start=None, stop=None, step=None,
                 rename=None, **kwargs):
        self.source = source
        self.format = format
        self._start  = start
        self._stop   = stop
        self._step   = step
        tmp_ar = AtomsReader(source, format, start, stop, step,
                             rename=rename,
                             **kwargs)
        list.__init__(self, list(iter(tmp_ar)))
        tmp_ar.close()

    def __getattr__(self, name):
        if name.startswith('__'):
            # don't override any special attributes
            raise AttributeError

        try:
            return self.source.__getattr__(name)
        except AttributeError:
            try:
                seq = [getattr(at, name) for at in iter(self)]
            except AttributeError:
                raise
            if seq == []:
                return None
            elif type(seq[0]) in (FortranArray, np.ndarray):
                return mockNDarray(*seq)
            else:
                return seq

    def __getslice__(self, first, last):
        return self.__getitem__(slice(first,last,None))

    def __getitem__(self, idx):
        if isinstance(idx, list) or isinstance(idx, np.ndarray):
            idx = np.array(idx)
            if idx.dtype.kind not in ('b', 'i'):
                raise IndexError("Array used for fancy indexing must be of type integer or bool")
            if idx.dtype.kind == 'b':
                idx = idx.nonzero()[0]
            res = []
            for i in idx:
                at = list.__getitem__(self,i)
                res.append(at)
        else:
            res = list.__getitem__(self, idx)
        if isinstance(res, list):
            res = AtomsList(res)
        return res

    def iterframes(self, reverse=False):
        if reverse:
            return reversed(self)
        else:
            return iter(self)

    @property
    def random_access(self):
        return True

    def sort(self, cmp=None, key=None, reverse=False, attr=None):
        """
        Sort the AtomsList in place. This is the same as the standard
        :meth:`list.sort` method, except for the additional `attr`
        argument. If this is present then the sorted list will be
        ordered by the :class:`Atoms` attribute `attr`, e.g.::
        
           al.sort(attr='energy')
        
        will order the configurations by their `energy` (assuming that
        :attr:`Atoms.params` contains an entry named `energy` for each
        configuration; otherwise an :exc:`AttributError` will be raised).
        """
        if attr is None:
            list.sort(self, cmp, key, reverse)
        else:
            if cmp is not None or key is not None:
                raise ValueError('If attr is present, cmp and key must not be present')
            list.sort(self, key=operator.attrgetter(attr), reverse=reverse)

    def apply(self, func):
        return np.array([func(at) for at in self])

def AtomsWriter(dest, format=None, **kwargs):
    """
    Returns a file-like object for writing Atoms to `dest` which
    should be either a filename or an initiliased output object.  If
    `format` is not given it is inferred from the file extension of
    `dest`. Example usage::

       out = AtomsWriter('out_file.xyz')
       for at in seq:
          out.write(at)
       out.close()
    """
    filename, dest, format = infer_format(dest, format, AtomsWriters)
    if format in AtomsWriters:
        writer = AtomsWriters[format](dest, **kwargs)
        return writer
    else:
        raise ValueError("Don't know how to write Atoms to format %r" % format)


class AtomsSequenceReader(object):
    """Read Atoms from a list of sources"""

    def __init__(self, sources, **kwargs):
        self.sources = sources
        self.readers = []
        self.lengths = []
        for source in sources:
            reader = AtomsReader(source, **kwargs)
            self.readers.append(reader)
            try:
                self.lengths.append(len(reader))
            except AttributeError:
                self.lengths.append(None)

    def __len__(self):
        if None in self.lengths:
            raise IndexError('One or more sources in %r do not support random access' % self.sources)
        return sum(self.lengths)

    def __getitem__(self, index):
        if None in self.lengths:
            raise IndexError('One or more sources in %r do not support random access' % self.sources)

        if isinstance(index,int) or isinstance(index, np.integer):
            if index < 0: index = index + sum(self.lengths)

            idx = 0
            for len, reader, source in zip(self.lengths, self.readers, self.sources):
                if index >= idx and index < idx+len:
                    at = reader[index-idx]
                    at.source = source
                    return at
                idx = idx+len

            raise IndexError('index %d out of range 0..%d' % (index, sum(self.lengths)))

        elif isinstance(index,slice):
            return [self[f] for f in range(*index.indices(sum(self.lengths)))]
        else:
            raise IndexError('indexing object should be either int or slice')


    def __iter__(self):
        for source, reader in zip(self.sources, self.readers):
            for at in reader:
                at.source = source
                yield at


class AtomsReaderCopier(object):
    def __init__(self, source):
        self.source = source

    def __len__(self):
        return len(self.source)

    def __iter__(self):
        for item in self.source:
            yield item.copy()

    def __getitem__(self, index):
        if isinstance(index, int) or isinstance(index, np.integer):
            return self.source[index].copy()
        elif isinstance(index, slice):
            return [at.copy() for at in self.source[index]]
        else:
            raise IndexError('indexing object should be either int or slice')            

    

def read(filename, **readargs):
    """
    Read Atoms from file `filename`

    File format is inferred from file extension, see :ref:`fileformats`.
    """
    return iter(AtomsReader(filename, **readargs)).next()


def write(filename, atoms, **writeargs):
    """
    Write `atoms` to the file `filename`

    File format is inferred from file extension, see :ref:`fileformats`.
    """
    if not isinstance(atoms, Atoms):
        atoms = Atoms(atoms)
    AtomsWriter(filename, **writeargs).write(atoms)


def read_dataset(dirs, pattern, **kwargs):
    """
    Read atomic configurations matching glob `pattern` from each of
    the directories in `dir` in turn. All kwargs are passed along
    to AtomsList constructor.

    Returns an dictionary mapping directories to AtomsList instances.
    """
    dataset = {}
    for dir in dirs:
        dataset[dir] = AtomsList(os.path.join(dir, pattern), **kwargs)
    return dataset


def time_ordered_series(source, dt=None):
    """
    Given a source of Atoms configurations, return a time ordered list of filename and frame references
    """

    if not isinstance(source, AtomsReader):
        if isinstance(source, basestring):
            source = time_ordered_glob(source)
        source = AtomsReader(source, range='empty')

    # traverse backwards over source, choosing most recently modified version of each frame
    revsource = reversed(source)
    last = revsource.next()
    current_time = last.time
    revorder = [(last.source, last.frame)]
    for at in revsource:
        try:
            if (at.time >= current_time) or (dt is not None and current_time - at.time < dt):
                continue
        except AttributeError:
            continue
        current_time = at.time
        revorder.append((at.source, at.frame))

    filenames = []
    # group first by filename
    for key, group in itertools.groupby(reversed(revorder), lambda x: x[0]):
        # annotate group with stride between frames
        group_with_stride = []
        group = list(group)
        for (s1, f1), (s2, f2) in zip(group, group[1:]):
            stride = f2 - f1
            group_with_stride.append((s1, f1, stride))
        group_with_stride.append(group[-1] + (stride,))

        # now group again by stride
        for key, group in itertools.groupby(group_with_stride, lambda x: x[2]):
            group = list(group)
            filenames.append('%s@%d:%d:%d' % (group[0][0], group[0][1], group[-1][1], group[0][2]))

    return filenames


@atoms_reader('traj')
@atoms_reader('cfg')
def ASEReader(source):
    """
    Helper routine to load from ASE trajectories
    """
    from ase.io import read
    from ase.atoms import Atoms as AseAtoms
    from quippy.elasticity import stress_matrix

    images = read(source, index=slice(None,None,None))
    if isinstance(images, AseAtoms):
        images = [images]
    
    for at in images:
        f = None
        try:
            f = at.get_forces()
        except RuntimeError:
            pass

        e = None
        try:
            e = at.get_potential_energy()
        except RuntimeError:
            pass

        v = None
        try:
            v = -stress_matrix(at.get_stress())*at.get_volume()
        except RuntimeError:
            pass

        at = Atoms(at) # convert to quippy Atoms
        if f is not None:
            at.add_property('force', f.T)
        if e is not None:
            at.params['energy'] = e
        if v is not None:
            at.params['virial'] = v
            
        yield at

class ASEWriter(object):

    def __init__(self, filename, translate=None):
        self.filename = filename
        self.translate = translate

    def write(self, at, **writeargs):
        from ase.io import write
        from ase.atoms import Atoms as AseAtoms
        ase_at = AseAtoms(at)

        if self.translate is not None:
            disp = np.dot(ase_at.cell, self.translate)
            ase_at.positions += disp

        write(self.filename, ase_at, **writeargs)

    def close(self):
        pass

AtomsWriters['traj'] = ASEWriter
AtomsWriters['cfg'] = ASEWriter

def dict2atoms(dct):
    from ase.db.core import dict2atoms
    from quippy.elasticity import stress_matrix
    
    at = dict2atoms(dct)

    f = None
    try:
        f = at.get_forces()
    except RuntimeError:
        pass

    e = None
    try:
        e = at.get_potential_energy()
    except RuntimeError:
        pass

    v = None
    try:
        v = -stress_matrix(at.get_stress())*at.get_volume()
    except RuntimeError:
        pass        

    at = Atoms(at) # convert to quippy Atoms
    if f is not None:
        at.add_property('force', f.T)
    if e is not None:
        at.params['energy'] = e
    if v is not None:
        at.params['virial'] = v

    # extract additional info from database
    at.params['id'] = dct['id']
    at.params['unique_id'] = dct['unique_id']
    if 'keywords' in dct:
        for key in dct['keywords']:
            at.params[key] = True
    if 'key_value_pairs' in dct:
        at.params.update(dct['key_value_pairs'])
    if 'data' in dct:
        for (key, value) in dct['data'].items():
            key = str(key) # avoid unicode strings
            value = np.array(value)
            if value.dtype.kind == 'U':
                value = value.astype(str)
            if value.dtype.kind != 'S':
                value = value.T
            try:
                at.add_property(key, value)
            except (TypeError, RuntimeError):
                at.params[key] = value

    return at

@atoms_reader('db')
@atoms_reader('json')
def ASEDatabaseReader(filename):
    from ase.db.core import connect

    index = None
    if isinstance(filename, basestring) and ('.json@' in filename or
                                             '.db@' in filename):
        filename, index = filename.rsplit('@', 1)
        if index.isdigit():
            index = int(index)

    conn = connect(filename)
            
    for dct in conn.select(index):
        at = dict2atoms(at)
        yield at
    
    
class ASEDatabaseWriter(object):
    """
    Write Atoms to :mod:`ase.db` database.
    """
    
    def __init__(self, filename, **kwargs):
        self.dbfile = filename
        self.kwargs = kwargs

    def write(self, at, **kwargs):
        import ase.db
        from ase.calculators.singlepoint import SinglePointCalculator

        all_kwargs = self.kwargs.copy()
        all_kwargs.update(kwargs)
        all_kwargs.update(at.params)
        
        energy = at.params.get('energy', None)
        forces = getattr(at, 'force', None)
        if forces is not None:
            forces = forces.T
        stress = at.params.get('virial', None)
        if stress is not None:
            stress = -stress.view(np.ndarray)/at.get_volume()
        
        keywords = []
        orig_calc = at.get_calculator()

        params = {}
        data = {}
        skip_params = ['energy', 'virial', 'calculator'] # filter out duplicate data
        for (key, value) in all_kwargs.items():
            key = key.lower()
            if key in skip_params:
                continue
            if value is None:
                keywords.append(key)
            if key == 'config_type':
                keywords.append(value)
                
            if (isinstance(value, int)   or isinstance(value, basestring) or
                isinstance(value, float) or isinstance(value, bool)):
                # scalar key/value pairs
                params[key] = value
            else:
                # more complicated data structures
                data[key] = value

        skip_arrays = ['numbers', 'positions', 'species']
        for (key, value) in at.arrays.items():
            if key in skip_arrays:
                continue
            key = key.lower()
            data[key] = value

        try:
            calc = SinglePointCalculator(atoms=at, energy=energy, forces=forces, stress=stress)
            if orig_calc is not None:
                calc.name = orig_calc.name
            else:
                calc.name = all_kwargs.get('calculator', '(unknown)')
            at.set_calculator(calc)

            database = ase.db.connect(self.dbfile)
            database.write(at, keywords, data=data, **params)
        finally:
            at.set_calculator(orig_calc)

AtomsWriters['db'] = ASEDatabaseWriter
AtomsWriters['json'] = ASEDatabaseWriter
