#!/usr/bin/env python

"""
pylab-style interface to quippy and AtomEye, designed for interactive use
"""

from numpy import *
from pylab import *
from quippy import *
from atomeye import *

from quippy.cinoutput import CInOutputReader
import os
import inspect

_viewers = {}
_current_viewer = None

class AtomEyeViewerMixin(AtomEyeViewer):
    """quippy-specific extensions to AtomEyeViewer"""
    
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
    def __init__(self, source, name=None, **kwargs):
        Atoms.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self, name)

    def update_source(self, source, **kwargs):
        self.read_from(source, **kwargs)
        self.redraw()

class AtomsReaderViewer(AtomsReader, AtomEyeViewerMixin):
    """
    Subclass of AtomsReader and AtomEyeViewer
    """
    def __init__(self, source, name=None, cache=True, **kwargs):
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self, name)

    def update_source(self, source, cache=True, **kwargs):
        self.close()
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        self.redraw()

class AtomsListViewer(AtomsList, AtomEyeViewerMixin):
    """
    Subclass of AtomsList and AtomEyeViewer
    """
    def __init__(self, source, name=None, **kwargs):
        AtomsList.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self, name)

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
        print 'Reusing viewer named %s for file %s' % (name, source)
        return (name, _viewers[name])
    else:
        print 'Creating viewer named %s for file %s' % (name, source)
        return (name, None)
        
def atoms(source, name=None, recycle=True, loadall=False, **kwargs):
    """
    Read atoms from `source` and open in an AtomEye viewer window.

    If not present, `name` is derived from the filename of `source`.

    If recycle=True (default), try to reuse an exising viewer window
    with the same name. Otherwise the name is made unique if necesary
    by appending a number.

    If loadall=False (the default) we use an AtomsReader to load the
    frames from the trajectory lazily (i.e., as required). Otherwise
    the entire file is read into an AtomsList.

    A new variable called `name` inserted into the parent stack frame.
    """

    name, viewer = find_viewer(source, name, recycle)

    if viewer is None:
        tmp_reader = AtomsReader(source, **kwargs)

        if tmp_reader.random_access and len(tmp_reader) == 1:
            viewer = AtomsViewer(source, name=name, **kwargs)
        else:
            if loadall:
                viewer = AtomsListViewer(source, name=name, **kwargs)
            else:
                viewer = AtomsReaderViewer(source, name=name, **kwargs)
    else:
        viewer.update_source(source, **kwargs)

    parent_frame = inspect.currentframe().f_back
    parent_frame.f_globals[viewer.name] = viewer
    return viewer

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

# Bring all AtomEyeViewerMixin methods into the top-level name space
for name, method in inspect.getmembers(AtomEyeViewerMixin, inspect.ismethod):
    if name.startswith('_'):
        continue
    setattr(sys.modules[__name__], name, current_viewer_method_wrapper(method))

del name, method, current_viewer_method_wrapper

if __name__ == '__main__':
    import optparse

    p = optparse.OptionParser(usage='%prog [options] <trajectory file>...')

    p.add_option('-f', '--frame', action='append', help='Initial frame to show. Can be given separately for each trajectory file. Default is last frame in file.', type=int)
    p.add_option('-p', '--property', action='append', help='Property to show. Can be given separately for each trajectory file')
    p.add_option('-a', '--arrows', action='append', help='Property to use to draw arrows. Can be given separately for each trajectory file')
    p.add_option('-l', '--load-view', action='store', help='Load view from AtomEye command script')
    p.add_option('-W', '--width', action='store', help="""Width of output movie, in pixels.""", type='int')
    p.add_option('-H', '--height', action='store', help="""Height of output movie, in pixels.""", type='int')
    p.add_option('-A', '--aspect', action='store', help="""Aspect ratio. Used if only one of --width or --height is given. Default 0.75.""", default=0.75, type='float')
    p.add_option('-R', '--rcut', action='append', help="""Following three arguments should be SYM1 SYM2 INCREMENT, to increment cutoff distance for SYM1-SYM2 bonds.""", nargs=3)
    p.add_option('-s', '--single', action='append', help="""Read next argument as a single frame, not a trajectory""", default=[])
    p.add_option('-L', '--loadall', action='store_true', help="""Load all frames from all trajectories into memory on startup.""")

    opt, args = p.parse_args()

    viewers = []
    for filename in opt.single:
        viewers.append(atoms(filename), single=True)
    for filename in args:
        viewers.append(atoms(filename, loadall=opt.loadall))

    if opt.frame is None:
        opt.frame = [-1 for viewer in viewers]

    show_args_list = [{} for viewer in viewers]
    for arg in ['frame', 'property', 'arrows']:
        values = getattr(opt, arg)
        if values is None:
            continue
        if len(values) == 1:
            values = [values[0] for traj in viewers]

        if len(values) != len(viewers):
            p.error('Number of -%s/--%s options does not match number of trajectory files' % (arg[0], arg))

        for show_args, value in zip(show_args_list, values):
            show_args[arg] = value

    for traj, show_args in zip(viewers, show_args_list):
        if opt.load_view is not None:
            print 'Loading view script %s' % opt.load_view
            traj.load_script(opt.load_view)

        if opt.rcut is not None:
            print 'Applying rcut patches %r' % opt.rcut
            for (sym1, sym2, rcut) in opt.rcut:
                traj.rcut_patch(sym1, sym2, float(rcut))

        if opt.width is not None or opt.height is not None:
            if opt.width  is None: opt.width = int(opt.height/opt.aspect)
            if opt.height is None: opt.height = int(opt.width*opt.aspect)
            print 'Setting viewer size to %d x %d pixels' % (opt.width, opt.height)
            traj.resize(opt.width, opt.height)

        if len(show_args) != 0:
            print 'Applying show_args=%r to trajectory %s' % (show_args, traj.name)
            traj.show(**show_args)

    del traj, show_args, filename, viewers, show_args_list, p, opt, args, arg, optparse

    from IPython import embed
    embed()


