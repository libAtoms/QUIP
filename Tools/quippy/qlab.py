#!/usr/bin/env python

"""
The :mod:`qlab` module provides a :mod:`pylab <matplotlib.pyplot>`\ -style
interface to :mod:`quippy` and :mod:`atomeye`, and is designed for interactive use,
especially in conjunction with `ipython <http://www.ipython.org>`_.

Viewing Atoms
-------------

The most important function defined in this module is :func:`view`, which is
used to open a new viewer or to update what is displayed in an existing viewer.
Here's an example session to give some ideas of how :mod:`qlab` can be used::

   from qlab import *    # Import all numpy, quippy and atomeye functions

   d = diamond(5.44, 14) # Make an 8-atom silicon unit cell
   view(d)               # Pops up an AtomEye window to display our cell
   capture('si8.png')    # Capture a screen shot to file "si8.png"

.. image:: si8.png
   :align: center

Once a viewer window is open, you can run all the methods of the
:class:`atomeye.AtomEyeViewer` class as top-level functions in the
which affect the current viewer. This means that you can type, e.g.::

    aux_property_coloring("z")

to colour the Atoms being displayed in the current viewer by their
atomic number :attr:`property <~.Atoms.atoms.properties>`, or::

    change_bgcolor((0, 0, 0))
    resize(640,480)

to change the background colour and size of the window. To redraw the
display, just repeat the :func:`view` command.

.. _qlab_atom_coloring:

Custom Atom Colouring
---------------------

The :meth:`~atomeye.AtomEyeViewer.aux_property_coloring` method is
overloaded to allow custom arrays to be visualised, e.g. to highlight
the atoms where :math:`x \ge a/2`::

   aux_property_coloring(d.pos[1,:] >= d.lattice[1,1]/2.)

.. image:: qlab-aux-property-colouring.png
   :align: center
   :width: 400

You can quickly find individual atoms or groups of atoms by passing an
integer or list of integers to :func:`aux_property_coloring`::

   aux_property_coloring(1) # highlight first atom (counting from 1 by default)

Note that these indices respect the :attr:`~atoms.Atoms.fortran_indexing`
attribute of the :class:`Atoms` object being viewed, i.e. they are
zero-based if ``gca().fortran_indexing`` is True and one-based if it's False.
This can be overriden by passing a ``fortran_indexing`` argument to :func:`view`::

   view(d, fortran_indexing=False)
   aux_property_coloring(0)  # highlight first atom (counting from zero)

To read off the properties of an individual atom, `right click` on it
in the AtomEye window. Again, the starting point for the indices printed
depends on :attr:`~atoms.Atoms.fortran_indexing`. The values used for
colouring the atom are in the ``_show`` property.

Here is a more advanced example showing how to draw arrows to visualise a vector
property (the forces on the atoms), and to colour the atoms by a component of the
stress tensor::

    d = diamond(5.44, 14)      # Make an 8-atom silicon unit cell
    s = supercell(d, 5, 5, 5)  # Make a supercell of the Si bulk
    view(s)                    # Visualise the system
    s.rattle(0.01)             # Randomise the atomic positions a little
    
    p = Potential('IP SW')     # Make a Stillinger-Weber potential calculator
    s.set_calculator(p)        # Associate it with our supercell

    f = s.get_forces()
    draw_arrows(f)  # draw vectors on the atoms to represent forces

    sigma = s.get_stresses()
    aux_property_coloring(sigma[:, 1, 1])  # colour atoms by sigma_yy

.. image:: qlab-force-stress.png
   :align: center
   :width: 600


Viewing Trajectories
--------------------

You can also pass a list of Atoms objects or a trajectory filename to
:func:`view` to visualise a sequence of frames. For example, here's a
sequence of silicon unit cells with increasing lattice constant::

   ds = [diamond(5.44+0.005*x, 14) for x in range(100)]
   view(ds)

Use `Insert` and `Delete` to move through the frames, or `Ctrl+Insert`
and `Ctrl+Delete` to jump to the first/last frame. Note the frame
number is printed in the window title. There are also :func:`first`,
:func:`last`, :func:`forward`, and :func:`backward` functions.

There is also a command line script `quippy` which starts an `ipython` shell
and imports everything from :mod:`qlab` automatically, and opens
viewers for any file given on the command line, e.g. from the shell ::

   $ quippy traj.nc

will fire up `ipython`, load everything from
:mod:`qlab` and then open a viewer for the file ``traj.nc``.

When working with large trajectory files, the :func:`clip_visible`
function is useful to restrict the number of atoms loaded from disk
and displayed, which can make visualising big systems much more
managable.

There is a :func:`gcv` function to get a reference to the current
viewer (short for"get current viewer", in analogy with
:func:`~matplotlib.pyplot.gcf` and :func:`~matplotlib.pyplot.gca`
in :mod:`pylab <matplotlib.pyplot>`), and a similarily named :func:`gcat`
("get current atoms") to get a reference to the current Atoms being viewed.
For a trajectory with multiple frames, this corresponds to the current frame.

You can set the current frame to 10 directly with::

   set_frame(10)

Or specify that the frame should advance by 5 steps each time `Insert`
is pressed with::

   set_delta(5)

If you would like to make a movie of your simulation, you can use
the :func:`render_movie` function::

   render_movie('movie.mp4')

This function renders each frame to a ``.jpg`` file, before combining the
snapshots with the `ffmpeg <http://www.ffmpeg.org/>`_ tool (which needs to
be installed for this to work). There is a `hook` function which is called at
each frame to allow change to be made to the Atoms object. For example, to run
a 1 ps MD and render the movie of every 10th frame::

   d = diamond(5.44, 14)           # Usual 8-atom cell
   s = supercell(d, 5, 5, 5)       # Make a supercell
   p = Potential('IP SW')          # Make a Stillinger-Weber potential
   s.set_cutoff(p.cutoff()+1.)     # Neighbour cutoff for atoms is pot cutoff + 1 A

   ds = DynamicalSystem(s)         # Construct a DynamicalSystem from our Atoms
   ds.rescale_velo(1000.)          # Set the initial temperature to 1000 K
   traj = ds.run(p, 1.0, 1000, save_interval=10) # Run 1 ps of dynamics

   view(traj)                      # Visualise the trajectory: 100 frames long
   aux_property_coloring('avg_ke') # Colour atoms by their time-averaged kinetic energy
   toggle_aux_property_thresholds_ridid() # Fix the colour scale during movie
   toggle_aux_property_thresholds_saturation() # Don't make any atoms invisible
   render_movie('traj.mp4')        # Make the movie

The movie generated by this script looks like this:

.. video:: traj 640 320

.. _qlab_select_atoms:

Selecting Atoms
---------------

The :func:`select_atoms` function is useful for graphically selecting a
group of atoms, e.g. to manually define an initial QM region for a
:ref:`QM/MM calculation <qmmm_tutorial>`. After invoking the function, right
clicking on atoms you want to select. They will be highlighted in red. When you
are done, press `Enter`. Here's an example `ipython` session::

    In [9]: view(s)
    Reusing viewer named at
    Out[9]: <AtomsViewer object at 0x10b9f19d0 fpointer=(-33474560, 32762)>
    In [10]: atoms_list = select_atoms()
    Right click to select atoms. Press ENTER to finish.
    indices = [ 224 265 261
     ]

If you want to add to the selection, run :func:`select_atoms` again with
the `reset` argument set to False. To customise the value applied to atoms
which are clicked on, use the `value` argument: e.g. to de-select atoms
clicked on by mistake, you could do::

    select_atoms(reset=False, value=False)

By default the selection is saved as a logical property named ``selection_mask``.
This name can be changed with the `markname` argument, e.g. to set the `hybrid`
property, you could could use::

    qm_list = select_atoms(markname='hybrid', value=HYBRID_ACTIVE_MARK, reset=False)


Multiple Viewers
----------------

By default a single viewer window is reused for individual Atoms
objects, while a new window is opened for each different trajectory
filename. You can override this with the `recycle` argument to
:func:`view`, e.g. to open a second window for a copy of `d`::

    d2 = d.copy()
    view(d2, recycle=False)

You can close a viewer by pressing `q` when it has the mouse focus.
Each viewer has a :attr:`~QuippyViewer.name` attribute which can be
use to get a reference to it with the :func:`get_viewer` function.

For trajectories read from files, the names are derived from the
filename. The default viewer for :class:`Atoms` objects is named `at`
and for the default viewers for :class:`AtomsList` and :class:`AtomsReader`
objects is called `al`. :func:`get_viewer_names` returns a list of the names af
all the viewers currently open.

Usually referring to the current viewer with :func:`gcv` is sufficient,
but you may want to change the default focus for console commands with
:func:`scv`, e.g. to set the current viewer to the one visualising the
file ``traj.nc``::

   scv(get_viewer('traj'))


.. note::

    :mod:`qlab` combines all the :mod:`numpy`, :mod:`quippy` and :mod:`atomeye` functions
    into a single namespace, so ``from qlab import *`` is roughly equivalent to::

        import numpy as np
        from numpy import *
        import quippy
        from quippy import *
        from atomeye import *
"""

import os
import sys
import inspect
import itertools

import numpy as np
from numpy import *

import quippy
from quippy import *
from atomeye import *

from quippy.cinoutput import CInOutputReader

__alldoc__ = ['QuippyViewer', 'AtomsViewer', 'AtomsListViewer', 'AtomsReaderViewer',
              'view', 'gcv', 'gcat', 'scv', 'get_viewer_names', 'get_viewers',
              'highlight_qm_region', 'redraw', 'run_command', 'run_script', 'close',
              'setp', 'save_script', 'toggle_coordination_coloring', 'translate',
              'shift_xtal', 'rotate', 'advance', 'shift_cutting_plane', 'change_bgcolor',
              'change_atom_r_ratio', 'change_bond_radius', 'change_view_angle_amplification',
              'toggle_parallel_projection', 'toggle_bond_mode', 'toggle_small_cell_mode',
              'normal_coloring', 'aux_property_coloring', 'central_symmetry_coloring',
              'change_aux_property_threshold', 'reset_aux_property_thresholds',
              'toggle_aux_property_thresholds_saturation', 'toggle_aux_property_thresholds_rigid',
              'rcut_patch', 'select_gear', 'cutting_plane', 'shift_cutting_plane_to_anchor',
              'delete_cutting_plane', 'flip_cutting_plane', 'capture', 'change_wireframe_mode',
              'change_cutting_plane_wireframe_mode', 'get_frame', 'set_frame', 'get_delta',
              'set_delta', 'first', 'last', 'forward', 'backward', 'load_atom_color',
              'load_aux', 'look_at_the_anchor', 'observer_goto', 'xtal_origin_goto',
              'find_atom', 'resize', 'change_aux_colormap', 'draw_arrows',
              'wait', 'get_visible', 'clip_visible', 'select_atoms', 'render_movie',
              'set_cutoffs']

_viewers = {}
_current_viewer = None

class QuippyViewer(AtomEyeViewer):
    """quippy-specific extensions to AtomEyeViewer"""

    _fields = ['verbose', 'echo', 'fortran_indexing',
               'filename', 'source', 'frame', 'name',
               'cache_mem_limit', 'block']
    
    def __init__(self, name, verbose=True, fortran_indexing=None):
        global _viewers, _current_viewer

        AtomEyeViewer.__init__(self, self, verbose=verbose, fortran_indexing=fortran_indexing)
        self._selection_markname = None
        self._selection_value = True
        self.name = name
        self._prev_viewer = _current_viewer
        _current_viewer = self
        _viewers[name] = _current_viewer
        self.set_cutoffs() # syncronise QUIP cutoffs with AtomEye
                        
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

            # allow ASE-style vector arrays
            if isinstance(_show, np.ndarray) and _show.shape == (at.n, 3):
               _show = np.transpose(_show)

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
            _current_viewer = self._prev_viewer

    def _click_hook(self, atoms, idx):
        if self._selection_markname is not None:
            atoms.properties[self._selection_markname][idx] = self._selection_value
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
        if reset or not at.has_property(markname):
            at.add_property(markname, False, overwrite=True)
        self._selection_markname = markname
        self._selection_value = value
        self.show(markname)
        saved_verbose = self.verbose
        self.verbose = False
        print 'Right click to select atoms. Press ENTER to finish.'
        print 'indices = [',
        raw_input()
        print ']'
        indices = list(at.properties[self._selection_markname].nonzero()[0])
        self._selection_markname = None
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

    def set_cutoffs(self, nneighb_only=True):
        """
        Set cutoffs for AtomEye bonds and coordination colouring

        Cutoff lengths are set to match quippy :attr:`~atoms.Atoms.nneightol`
        (if `nneighb_only` is True, the default) or :attr:`~atoms.Atoms.cutoff`
        (otherwise).
        """
        at = self.gcat()
        seen = []
        for Z1 in set(at.z):
            for Z2 in set(at.z):
                if (min(Z1,Z2), max(Z1, Z2)) in seen:
                    continue
                seen.append((min(Z1, Z2), max(Z1, Z2)))
                sym1, sym2 = ElementName[Z1], ElementName[Z2]
                print sym1, sym2, 
                if nneighb_only:
                    cutoff = at.nneightol*bond_length(Z1, Z2)
                    print 'nneigbb', cutoff
                elif at.use_uniform_cutoff:
                    cutoff = at.cutoff
                    print 'uniform', cutoff
                else:
                    cutoff = at.cutoff*bond_length(Z1, Z2)
                    print 'relative', cutoff
                self.rcut_patch(sym1, sym2, cutoff, absolute=True)
    

class AtomsViewer(Atoms, QuippyViewer):
    """
    Subclass of Atoms and AtomEyeViewer
    """
    def __init__(self, source=None, name=None, verbose=True, fortran_indexing=None, **kwargs):
        if fortran_indexing is None:
            fortran_indexing = False
            if hasattr(source, 'fortran_indexing'):
                fortran_indexing = source.fortran_indexing
        Atoms.__init__(self, fortran_indexing=fortran_indexing)
        if isinstance(source, Atoms):
            self.shallow_copy_from(source)
        else:
            self.read_from(source, **kwargs)
        QuippyViewer.__init__(self, name, verbose=verbose,
                              fortran_indexing=fortran_indexing,
                              **kwargs)

    def gcat(self, update=False):
        return self

    def scat(self, atoms, frame=None):
        pass

    def update_source(self, source, **kwargs):
        if not isinstance(source, Atoms):
            source = Atoms(source)
        self.shallow_copy_from(source)
        self.redraw()

    def reload(self):
        self.update_source(self.source)

    def copy(self):
        return Atoms.copy(self)


class AtomsReaderViewer(AtomsReader, QuippyViewer):
    """
    Subclass of AtomsReader and AtomEyeViewer
    """
    def __init__(self, source=None, name=None, cache=True, verbose=True,
                 fortran_indexing=None, **kwargs):
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, fortran_indexing=fortran_indexing, **kwargs)
        QuippyViewer.__init__(self, name, verbose=verbose, fortran_indexing=fortran_indexing)

    def update_source(self, source, cache=True, **kwargs):
        self.close()
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, fortran_indexing=self.fortran_indexing,
                             **kwargs)
        self.redraw()

    def reload(self):
        self.update_source(self.source)
        self.last()


class AtomsListViewer(AtomsList, QuippyViewer):
    """
    Subclass of AtomsList and AtomEyeViewer
    """
    def __init__(self, source=None, name=None, fortran_indexing=None, **kwargs):
        AtomsList.__init__(self, source, fortran_indexing=fortran_indexing, **kwargs)
        QuippyViewer.__init__(self, name, fortran_indexing=fortran_indexing)

    def update_source(self, source, **kwargs):
        del self[:]
        AtomsList.__init__(self, source, fortran_indexing=self.fortran_indexing, **kwargs)
        self.redraw()

    def reload(self):
        self.update_source(self.source)
        self.last()


def find_viewer(source, name=None, recycle=True):
    global _viewers, _current_viewer

    if name is None:
        if hasattr(source, 'name'):
            name = source.name
        elif isinstance(source, basestring):
            name = os.path.splitext(os.path.basename(source))[0]
            name = name.replace('-','_').replace('.','_').replace('*','').replace('?','')
        elif hasattr(source, '__iter__'):
            name = 'al'
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
        

def view(source, name=None, recycle=True, loadall=False, inject=True,
         fortran_indexing=None, **kwargs):
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
        tmp_reader = AtomsReader(source, fortran_indexing=fortran_indexing, **kwargs)

        if isinstance(source, Atoms) or (tmp_reader.random_access and len(tmp_reader) == 1):
            viewer = AtomsViewer(source, name, fortran_indexing=fortran_indexing, **kwargs)
        else:
            if loadall or not tmp_reader.random_access:
                viewer = AtomsListViewer(source, name=name,
                                         fortran_indexing=fortran_indexing,
                                         **kwargs)
            else:
                viewer = AtomsReaderViewer(source, name=name,
                                           fortran_indexing=fortran_indexing,
                                           **kwargs)
    else:
        viewer.update_source(source, **kwargs)

    if inject:
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

def gcat():
    """
    Return the current Atoms object being visualised by the current viewer
    """
    return gcv().gcat()

def scv(viewer):
    """
    Set the current AtomEye viewer to `viewer`.
    """
    global _current_viewer
    _current_viewer = viewer

def get_viewer(name):
    """
    Return the viewer identified by `name`
    """
    global _viewers
    return _viewers[name]

def get_viewer_names():
    """
    Return the current list of viewer names
    """
    global _viewers
    return _viewers.keys()

def get_viewers():
    """
    Return the current list of viewers
    """
    global _viewers
    return _viewers.values()

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
        save_block = gcv().block
        gcv().block = True
        redraw()
        rcut_patch('Si', 'Si', +0.3)
        rcut_patch('Al', 'Al', -0.55)
        run_command('change_normal_color 13 0.0 0.0 0.7 1.2')
        run_command('change_normal_color 5 0.9 0.4 0 1.5')
        run_command('change_normal_color 7 0.0 0.7 0.7 0.7')
        highlight_qm_region.first_time = False
        gcv().block = save_block
    redraw()

highlight_qm_region.first_time = True

def redraw():
    """
    Redraw current AtomEye window, keeping Atoms and settings the same.
    """
    gcv().redraw()

def run_command(command):
    """
    Run a command in current AtomEye thread.

    The command is queued for later execution, unless :attr:`block` is True.

    Parameters
    ----------

    command : string
       The command to pass to AtomEye
    """
    gcv().run_command(command)

def run_script(script):
    """
    Run commands from the file script, in a blocking fashion.
    """
    gcv().run_script(script)

def close():
    """
    Close the current viewer window.
    """
    gcv().close()

def setp(self, key, value):
    """
    Run the AtomEye command "set key value"
    """
    gcv().setp(key, value)

def save_script(filename):
    """
    Save AtomEye viewer settings to a file.
    """
    gcv().save(filename)

def toggle_coordination_coloring():
    """
    Turn on or off colouring by coordination number (key "k")
    """
    gcv().toggle_coordination_coloring()

def translate(axis, delta):
    """
    Translate system along `axis` by an amount `delta` (key "Ctrl+left/right/up/down")
    """
    gcv().translate(axis, delta)

def shift_xtal(axis, delta):
    """
    Shift crystal within periodic boundaries along `axis` by `delta` (key "Shift+left/right/up/down").
    """
    gcv().shift_xtal(axis, delta)

def rotate(axis, theta):
    """
    Rotate around `axis` by angle `theta`.
    """
    gcv().rotate(axis, theta)

def advance(delta):
    """
    Move the camera forward by `delta`.
    """
    gcv().advance(delta)

def shift_cutting_plane(delta):
    """
    Shift the current cutting plane by an amount `delta`.
    """
    gcv().shift_cutting_plane(delta)

def change_bgcolor(color):
    """
    Change the viewer background colour to `color`, which should be a RGB tuple with three floats in range 0..1.
    """
    gcv().change_bgcolor(color)

def change_atom_r_ratio(delta):
    """
    Change the size of the balls used to draw the atoms by `delta`.
    """
    gcv().change_atom_r_ratio(delta)

def change_bond_radius(delta):
    """
    Change the radius of the cylinders used the draw bonds by `delta`.
    """
    gcv().change_bond_radius(delta)

def change_view_angle_amplification(delta):
    """
    Change the amplification of the view angle by `delta`.
    """
    gcv().change_view_angle_amplification(delta)

def toggle_parallel_projection():
    """
    Toggle between parallel and perspective projections.
    """
    gcv().toggle_parallel_projection()

def toggle_bond_mode():
    """
    Turn on or off bonds.
    """
    gcv().toggle_bond_mode()

def toggle_small_cell_mode():
    """
    Toggle between two different behaviours for when cell is smaller than r_cut/2:
     1. clip cell - some neigbours may be lost (default)
     2. replicate cell along narrow directions
    """
    gcv().toggle_small_cell_mode()

def normal_coloring():
    """
    Return to normal colouring of the atoms (key "o").
    """
    gcv().normal_coloring()

def aux_property_coloring(auxprop):
    """
    Colour the currently viewed atoms according to `auxprop`.

    Overloaded to allow 
    See :ref:`qlab_atom_coloring` for more details and examples.

    Parameters
    ----------
    auxprop : str, array_like, int or list
       Values to use to colour the atoms. Should be either the
       name of a scalar field entry in :attr:`~.Atoms.properties`
       (or equivalently, :attr:`~Atoms.arrays`) such as ``"charge"``,
       a float, int or bool array of shape ``(len(gcat()),)``, or an
       atom index or list of atom indices to highlight particular atoms.
    """
    gcv().aux_property_coloring(auxprop)

def central_symmetry_coloring():
    """
    Colour atoms by centro-symmetry parameter.
    """
    gcv().central_symmetry_coloring()

def change_aux_property_threshold(lower, upper):
    """
    Change the lower and upper aux property thresholds.
    """
    gcv().change_aux_property_threshold(lower, upper)

def reset_aux_property_thresholds():
    """
    Reset aux property thresholds to automatic values.
    """
    gcv().reset_aux_property_thresholds()

def toggle_aux_property_thresholds_saturation():
    """
    Toggle between saturated colouring and invisibility for values outside aux prop thresholds.
    """
    gcv().toggle_aux_property_thresholds_saturation()

def toggle_aux_property_thresholds_rigid():
    """
    Toggle between floating and rigid aux property thresholds when moving between frames
    """
    gcv().toggle_aux_property_thresholds_rigid()

def rcut_patch(sym1, sym2, value, absolute=False):
    """
    Change the cutoff distance for `sym1`--`sym2` bonds by `delta`.

    e.g. to increase cutoff for Si-Si bonds by 0.5 A use::

         viewer.rcut_patch('Si', 'Si', 0.5)

    With `absolute` set to True, `value` is used to set the
    absolute cutoff distance for `sym1`--`sym2` bonds, e.g.::

         viewer.rcut_patch('Si', 'Si', 2.50, True)
    """
    gcv().rcut_patch(sym1, sym2, value, absolute)

def select_gear(gear):
    """
    Change the AtomEye gear to `gear`

    Equivalent to pressing the one of the numeric keys 0..9
    """
    gcv().select_gear(gear)

def cutting_plane(n, d, s):
    """
    Create a new cutting plane with index `n`, normal `d`, and fractional displacement `s`.
    """
    gcv().cutting_plane(n, d, s)

def shift_cutting_plane_to_anchor(n):
    """
    Move the cutting plane with index `n` to the anchor
    """
    gcv().shift_cutting_plane_to_anchor(n)

def delete_cutting_plane(n):
    """
    Delete the cutting plane with index `n`
    """
    gcv().delete_cutting_plane(n)

def flip_cutting_plane(n):
    """
    Flip the cutting plane with index `n`
    """
    gcv().flip_cutting_plane(n)

def capture(filename, resolution=None):
    """
    Render the current view to image `filename`

    Format is determined from file extension: .png, .jpeg, or .eps.
    """
    gcv().capture(filename, resolution)

def change_wireframe_mode():
    """
    Change the display mode for the unit cell box.

    Equivalent to pressing the `i` key.
    """
    gcv().change_wireframe_mode()

def change_cutting_plane_wireframe_mode():
    """
    Change the display mode for cutting planes
    """
    gcv().change_cutting_plane_wireframe_mode()

def get_frame():
    """
    Get index of frame currently being viewed
    """
    return gcv().frame

def set_frame(frame):
    """
    Set current frame index to `frame`
    """
    gcv().frame = frame

def get_delta():
    """
    Get frame increment rate
    """
    return gcv().delta

def set_delta(delta):
    """
    Set frame increment rate
    """
    gcv().delta = delta

def first():
    """
    Show the first frame (frame 0).
    """
    gcv().first()

def last():
    """
    Show the last frame, i.e. len(gcv())-1
    """
    gcv().first()

def forward(delta=None):
    """
    Move forward by `delta` frames (default value is gcv().delta).
    """
    gcv().forward(delta)

def backward(delta=None):
    """
    Move backward by `delta` frames (default values is gcv().delta).
    """
    gcv().backward(delta)

def load_atom_color(filename):
    """
    Load atom colours from a .clr file.
    """
    gcv().load_atom_color(filename)

def load_aux(filename):
    """
    Load aux property values from a .aux file.
    """
    gcv().load_aux(filename)

def look_at_the_anchor():
    """
    Equivalent to pressing the `a` key
    """
    gcv().look_at_the_anchor()

def observer_goto():
    """
    Prompt for fractional position and move the observer there

    Equivalent to pressing the `g` key.
    """
    gcv().observer_goto()

def xtal_origin_goto(s):
    """
    Move the crystal origin to fractional coordinates `s`

    For example, use ``s=[0.5, 0.5, 0.5]`` to shift by half the cell along
    the :math:`\mathbf{a}`, :math:`\mathbf{b}` and :math:`\mathbf{c}`
    lattice vectors.
    """
    gcv().xtal_origin_goto(s)

def find_atom(i):
    """
    Set the anchor to the atom with index `i`.
    """
    gcv().find_atom(i)

def resize(width, height):
    """
    Resize the current window to `width` x `height` pixels.
    """
    gcv().resize(width, height)

def change_aux_colormap(n):
    """
    Select the `n`\ -th auxiliary property colourmap. 
    """
    gcv().change_aux_colormap(n)

def draw_arrows(property, scale_factor=0.0, head_height=0.1,
                head_width=0.05, up=(0.0,1.0,0.0)):
    """
    Draw arrows on each atom, based on a vector property

    Parameters
    ----------
    property : string
       Name of the array to use for arrow vectors.
       Use ``None`` to turn off previous arrows.
    scale_factor : float
       Override length of arrows. 1 unit = 1 Angstrom; default
       value of 0.0 means autoscale.
    head_height : float
       Specify height of arrow heads in Angstrom. 
    head_width : float
    up : 3-vector (tuple, list or array)
       Specify the plane in which the arrow heads are
       drawn. Arrows are drawn in the plane which is common
       to their direction and this vector.
       Default is ``[0.,1.,0.]``.
    """
    gcv().draw_arrows(property, scale_factor, head_height,
                      head_width, up)

def wait():
    """Sleep until current AtomEye viewer has finished processing all queued events."""
    gcv().wait()

def get_visible():
    """Return list of indices of atoms currently visible in the current viewer."""
    return gcv().get_visible()


def clip_visible(orig_index=True):
    """
    Remove atoms outside the visible window from the Atoms object. Also sets indices for frames not yet loaded from disk.
    """
    gcv().clip_visible(orig_index)

def select_atoms(reset=True, markname='selection_mark', value=True):
    """
    Select atoms by clicking on them. Returns a list of atom indices.

    Specify reset=False to modify an existing property. The name of the
    property is `markname` (default "selection_mark") and the value of
    clicked atoms is given by the `value` argument (default True).

    See also :ref:`qlab_select_atoms`
    """
    gcv().select_atoms(reset, markname, value)

def render_movie(moviefile, start=None, stop=None, step=None, hook=None,
                 offset=0, encoder='ffmpeg -i %s -r 25 -b 30M %s'):
    """
    Render a movie for the trajectory.
    """
    gcv().render_movie(moviefile, start, stop, step, hook, offset, encoder)

def set_cutoffs(nneighb_only=True):
    """
    Set cutoffs for AtomEye bonds and coordination colouring

    Cutoff lengths are set to match quippy :attr:`~atoms.Atoms.nneightol`
    (if `nneighb_only` is True, the default) or :attr:`~atoms.Atoms.cutoff`
    (otherwise).
    """
    gcv().set_cutoffs(nneighb_only)
