"""
pylab-style interface to quippy and AtomEye, designed for interactive use
"""

from numpy import *
from pylab import *
from quippy import *
from atomeye import *

from quippy.cinoutput import CInOutputReader
import os
import functools
import inspect

_viewers = {}
_current_viewer = None

class AtomEyeViewerMixin(AtomEyeViewer):
    """quippy-specific extensions to AtomEyeViewer"""
    
    def __init__(self, source):
        global _viewers, _current_viewer

        AtomEyeViewer.__init__(self, self)
        self.selection_mark = None
        name = os.path.splitext(os.path.basename(source))[0]
        name = name.replace('-','_').replace('.','_').replace('*','').replace('?','')
        print 'Creating %s named %s for file %s' % (self.__class__.__name__, name, source)
        _current_viewer = self
        _viewers[name] = _current_viewer
        self.name = name

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
                at.remove_auxprop('_show')
            at.add_property('_show', _show)
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

    def _enter_hook(self, atoms):
        def hook(at):
            if self.is_alive:
                self.redraw()
        atoms.update_redraw_hook = hook
        atoms.add_hook(atoms.update_redraw_hook)

    def _exit_hook(self, atoms):
        atoms.remove_hook(atoms.update_redraw_hook)

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

        at = self.gca()
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

    def select_atoms(self, reset=True, markname='selection_mark', value=True):
        """
        Select atoms by clicking on them. Returns a list of atom indices.

        Specify reset=False to modify an existing property. The name of the
        property is `markname` (default "selection_mark") and the value of
        clicked atoms is given by the `value` argument (default True).
        """
        at = self.gca()
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

    def render_movie(self, moviefilename):
        """
        Render a movie for the trajectory. Not yet implemented.
        """
        raise NotImplementedError
        

class AtomsViewer(Atoms, AtomEyeViewerMixin):
    """
    Subclass of Atoms and AtomEyeViewer
    """
    def __init__(self, source, **kwargs):
        Atoms.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self, source)

class AtomsReaderViewer(AtomsReader, AtomEyeViewerMixin):
    """
    Subclass of AtomsReader and AtomEyeViewer
    """
    def __init__(self, source, cache=True, **kwargs):
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self, source)

class AtomsListViewer(AtomsList, AtomEyeViewerMixin):
    """
    Subclass of AtomsList and AtomEyeViewer
    """
    def __init__(self, source, **kwargs):
        AtomsList.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self, source)
        
def atoms(filename, single=False, loadall=False, **kwargs):
    """
    Read atoms from `filename` and create a viewer for it.

    If single=False (default) file is treated as a trajectory. Otherwise
    only a single frame is loaded.

    If loadall=False (the default) we use an AtomsReader to load the
    frames from the trajectory lazily (i.e., as required). Otherwise
    the entire file is read into an AtomsList.

    A new variable is inserted into the parent frame.
    """
    if single:
        a = AtomsViewer(filename, **kwargs)
    else:
        if loadall:
            a = AtomsListViewer(filename, **kwargs)
        else:
            a = AtomsReaderViewer(filename, **kwargs)
    parent_frame = inspect.currentframe().f_back
    parent_frame.f_globals[a.name] = a
    return a

def reset():
    """
    Close all currently open AtomEye viewers
    """
    global _viewers, _current_viewer
    for name, traj in _viewers.iteritems():
        print 'Closing trajectory viewer %s' % name
        traj.close()
    _viewers = {}
    _current_viewer = None

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
    return functools.update_wrapper(call, method)

# Bring all AtomEyeViewerMixin methods into the top-level name space
for name, method in inspect.getmembers(AtomEyeViewerMixin, inspect.ismethod):
    if name.startswith('_'):
        continue
    setattr(sys.modules[__name__], name, current_viewer_method_wrapper(method))

del name, method, current_viewer_method_wrapper
