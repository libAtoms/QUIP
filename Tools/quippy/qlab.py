#!/usr/bin/env python

"""
pylab-style interface to quippy and AtomEye, designed for interactive use
"""

from numpy import *
from quippy import *
from atomeye import *

from quippy.cinoutput import CInOutputReader
import os
import sys
import inspect

_viewers = {}
_current_viewer = None

class QuippyViewer(AtomEyeViewer):
    """quippy-specific extensions to AtomEyeViewer"""

    _fields = ['verbose', 'echo', 'fortran_indexing',
               'filename', 'source', 'frame', 'name',
               'cache_mem_limit', 'block']
    
    def __init__(self, name):
        global _viewers, _current_viewer

        AtomEyeViewer.__init__(self, self, fortran_indexing=True)
        self.selection_mark = None
        self.name = name
        _current_viewer = self
        _viewers[name] = _current_viewer
                        
        return self

    def _property_hook(self, at, auxprop):
        if not hasattr(at, 'properties'):
            return

        if auxprop is not None and not isinstance(auxprop,str):
            if isinstance(auxprop,int):
                _show = [i == auxprop for i in at.indices]
            elif isinstance(auxprop, list) or isinstance(auxprop, tuple) or isinstance(auxprop, set):
                _show = [i in auxprop for i in at.indices]
            else:
                _show = auxprop

            if at.has_property('_show'):
                del at.properties['_show']
            at.add_property('_show', _show, overwrite=True)
            auxprop = '_show'

        # Ensure auxprop is one of the first AUX_PROPERTY_COLORING
        # columns by swapping it with an inessential property if necessary.
        prop_idx = 0
        for key, value in at.properties.iteritems():
            ncols = len(value.shape) == 2 and value.shape[0] or 1
            prop_idx += ncols
            if key.lower() == auxprop.lower():
                break
        else:
            raise ValueError('Unknown Atoms property %s' % auxprop)

        if prop_idx >= self.CONFIG_MAX_AUXILIARY:
            for swapprop in at.properties:
                if swapprop.lower() not in ['pos', 'z', 'species']:
                    break
            at.properties.swap(auxprop, swapprop)

        return auxprop

    # def _enter_hook(self, atoms):
    #     def hook(at):
    #         if self.is_alive:
    #             self.redraw()
    #     atoms.update_redraw_hook = hook
    #     atoms.add_hook(atoms.update_redraw_hook)

    # def _exit_hook(self, atoms):
    #     atoms.remove_hook(atoms.update_redraw_hook)

    def _close_hook(self):
        global _viewers, _current_viewer
        if self.name in _viewers:
            del _viewers[self.name]
        if _current_viewer is self:
            _current_viewer = None

    def _click_hook(self, atoms, idx):
        if self.selection_mark is not None:
            self.selection_mark[idx] = self.selection_value
            print idx,
            sys.stdout.flush()
            self.redraw()
        elif self.verbose:
            print
            atoms.print_atom(idx)
            sys.stdout.flush()

    def _redraw_hook(self, atoms):
        print 'Properties:'
        for key, value  in atoms.properties.iteritems():
            print '%-10s shape %r' % (key, value.shape)
        print '\nParams:'
        for (key, value) in atoms.params.iteritems():
            print '%-20s = %r' % (key, value)

    def show(self, property=None, frame=None, arrows=None):
        """
        Update what is shown in this AtomEye viewer window.

        `property` should be the name of the auxiliary property used to colour the atoms (e.g. "charge")
        `frame` is the (zero-based) index of the frame to show.
        `arrows` is the name of a vector property to use to draw arrows on the atoms (e.g. "force")

        When called with no arguments, show() is equivalent to redraw().
        """
        global _current_viewer
        _current_viewer = self
        AtomEyeViewer.show(self, self, property, frame, arrows)
    
    def clip_visible(self, orig_index=True):
        """
        Remove atoms outside the visible window from the Atoms object. Also sets indices for frames not yet loaded from disk.
        """
        indices = self.get_visible()
        print 'Clipping view to include %d atoms' % len(indices)

        at = self.gcat()
        mask = fzeros(len(at), dtype=bool)
        mask[:] = True
        mask[indices] = False
        at.remove_atoms(mask=mask)

        if orig_index:
            at.add_property('orig_index',
                            logical_not(mask).nonzero()[0],
                            overwrite=True)
        self.redraw()

        if hasattr(self, 'reader') and isinstance(self.reader, CInOutputReader):
            self.reader.source.indices = indices
        elif hasattr(self, 'reader') and hasattr(self.reader, '__iter__'):
            for r in self.reader.readers:
                if hasattr(r, 'reader') and isinstance(r.reader, CInOutputReader):
                    r.reader.source.indices = indices

    def select_atoms(self, reset=True, markname='selection_mark', value=True):
        """
        Select atoms by clicking on them. Returns a list of atom indices.

        Specify reset=False to modify an existing property. The name of the
        property is `markname` (default "selection_mark") and the value of
        clicked atoms is given by the `value` argument (default True).
        """
        at = self.gcat()
        if reset:
            at.add_property(markname, False, overwrite=True)
        self.selection_mark = getattr(at, markname)
        self.selection_value = value
        self.show(markname)
        saved_verbose = self.verbose
        self.verbose = False
        print 'Click to select atoms. Press ENTER to finish.'
        print 'indices = [',
        raw_input()
        print ']'
        indices = list(self.selection_mark.nonzero()[0])
        self.selection_mark = None
        self.verbose = saved_verbose
        return indices

    def render_movie(self, moviefile, start=None, stop=None, step=None, hook=None,
                     offset=0, encoder='ffmpeg -i %s -r 25 -b 30M %s'):
        """
        Render a movie for the trajectory.
        """

        if start is not None or stop is not None or step is not None:
            frames = range(*slice(start, stop, step).indices(len(self)))
        else:
            frames = range(len(self))

        basename, ext = os.path.splitext(moviefile)
        out_fmt = '%s%%05d.jpg' % basename
        
        for frame in frames:
            self.show(frame=frame)
            if hook is not None:
                self.wait()
                hook(self.gcat())
                self.redraw()
            self.capture(out_fmt % (frame + offset))
            self.wait()

        print 'Encoding movie...'
        os.system(encoder % (out_fmt, moviefile))

    def copy(self):
        return atoms(self, recycle=False, inject=False)

    def __getinitargs__(self):
        if hasattr(self.atoms, 'filename'):
            return (self.atoms.filename,)
        elif hasattr(self.atoms, 'source'):
            return (self.atoms.source)
        else:
            raise ValueError("can't work out how to pickle Atoms")

    def __getstate__(self):
        state = AtomEyeViewer.__getstate__(self)
        return state

    def __setstate__(self, state):
        AtomEyeViewer.__setstate__(self, state) 
    

class AtomsViewer(Atoms, QuippyViewer):
    """
    Subclass of Atoms and AtomEyeViewer
    """
    def __init__(self, source=None, name=None, **kwargs):
        Atoms.__init__(self, source, **kwargs)
        QuippyViewer.__init__(self, name)

    def gcat(self, update=False):
        return self

    def scat(self, atoms, frame=None):
        pass

    def update_source(self, source, **kwargs):
        self.read_from(source, **kwargs)
        self.redraw()

    def copy(self):
        return Atoms.copy(self)

class AtomsReaderViewer(AtomsReader, QuippyViewer):
    """
    Subclass of AtomsReader and AtomEyeViewer
    """
    def __init__(self, source=None, name=None, cache=True, **kwargs):
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        QuippyViewer.__init__(self, name)

    def update_source(self, source, cache=True, **kwargs):
        self.close()
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        self.redraw()

class AtomsListViewer(AtomsList, QuippyViewer):
    """
    Subclass of AtomsList and AtomEyeViewer
    """
    def __init__(self, source=None, name=None, **kwargs):
        AtomsList.__init__(self, source, **kwargs)
        QuippyViewer.__init__(self, name)

    def update_source(self, source, **kwargs):
        del self[:]
        AtomsList.__init__(self, source, **kwargs)
        self.redraw()

def find_viewer(source, name=None, recycle=True):
    global _viewers, _current_viewer

    if name is None:
        if hasattr(source, 'name'):
            name = source.name
        elif isinstance(source, basestring):
            name = os.path.splitext(os.path.basename(source))[0]
            name = name.replace('-','_').replace('.','_').replace('*','').replace('?','')
        else:
            name = 'at'

    if name in _viewers and not recycle:
        # find a unique name
        n = 1
        new_name = name
        while new_name in _viewers:
            n += 1
            new_name = '%s_%d' % (name, n)
        name = new_name

    if name in _viewers:
        print 'Reusing viewer named %s' % name
        scv(_viewers[name])
        return (name, _viewers[name])
    else:
        print 'Creating viewer named %s' % name
        return (name, None)
        
def read(source, name=None, recycle=True, loadall=False, inject=True, **kwargs):
    """
    Read atoms from `source` and open in an AtomEye viewer window.

    If not present, `name` is derived from the filename of `source`.

    If recycle is true (default), try to reuse an exising viewer window
    with the same name. Otherwise the name is made unique if necesary
    by appending a number.

    If loadall is false (default) we use an AtomsReader to load the
    frames from the trajectory lazily (i.e., as required). Otherwise
    the entire file is read into an AtomsList.

    If inject is true (default), a new variable called `name` is injected
    into the parent stack frame.
    """

    name, viewer = find_viewer(source, name, recycle)

    if viewer is None:
        tmp_reader = AtomsReader(source, **kwargs)

        if isinstance(source, Atoms) or (tmp_reader.random_access and len(tmp_reader) == 1):
            viewer = AtomsViewer(source, name, **kwargs)
        else:
            if loadall:
                viewer = AtomsListViewer(source, name=name, **kwargs)
            else:
                viewer = AtomsReaderViewer(source, name=name, **kwargs)
    else:
        viewer.update_source(source, **kwargs)

    if inject:
        parent_frame = inspect.currentframe().f_back
        parent_frame.f_globals[viewer.name] = viewer
    return viewer

# make 'view' a synonym for 'load'
view = read

def gcv():
    """
    Return the current (most recently created or used) AtomEye viewer instance
    """
    global _current_viewer
    if _current_viewer is None:
        raise ValueError('No viewers are currently open!')
    return _current_viewer

def scv(viewer):
    """
    Set the current AtomEye viewer to `viewer`.
    """
    global _current_viewer
    _current_viewer = viewer

def iterviewers():
    """
    Generator which yields (name, viewer) pairs.
    """
    global _viewers
    for (name, viewer) in _viewers.iteritems():
        yield (name, viewer)

def current_viewer_method_wrapper(method, *args, **kwargs):
    def call(*args, **kwargs):
        return getattr(gcv(), method.func_name)(*args, **kwargs)
    try:
        import functools
        return functools.update_wrapper(call, method)
    except ImportError:
        return call

# Bring all QuippyViewer methods into the top-level name space
for name, method in inspect.getmembers(QuippyViewer, inspect.ismethod):
    if name.startswith('_'):
        continue
    setattr(sys.modules[__name__], name, current_viewer_method_wrapper(method))

del name, method, current_viewer_method_wrapper


def highlight_qm_region(at=None, run_suffix=''):
    """
    Highlight QM region by replacing Si atoms with Al,
    and O atoms with N, and changing colour of QM atoms to dark blue. Can be used as
    a hook function to render_movie().

    If at is None, uses Atoms associated with current viewer
    (i.e., at = gcat()).
    """
    if at is None:
        at = gcat()
    hybrid_mark = getattr(at, 'hybrid_mark'+run_suffix)
    at.z[(at.z == 14) & (hybrid_mark == 1)] = 13
    at.z[(at.z == 8) & (hybrid_mark == 1)] = 7
    at.set_atoms(at.z)
    if highlight_qm_region.first_time:
        redraw()
        wait()
        rcut_patch('Si', 'Si', +0.3)
        rcut_patch('Al', 'Al', -0.55)
        run_command('change_normal_color 13 0.0 0.0 0.7 1.2')
        run_command('change_normal_color 5 0.9 0.4 0 1.5')
        run_command('change_normal_color 7 0.0 0.7 0.7 0.7')
        highlight_qm_region.first_time = False
    redraw()


highlight_qm_region.first_time = True
