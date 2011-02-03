.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
.. HQ X
.. HQ X   quippy: Python interface to QUIP atomistic simulation library
.. HQ X
.. HQ X   Copyright James Kermode 2010
.. HQ X
.. HQ X   These portions of the source code are released under the GNU General
.. HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
.. HQ X
.. HQ X   If you would like to license the source code under different terms,
.. HQ X   please contact James Kermode, james.kermode@gmail.com
.. HQ X
.. HQ X   When using this software, please cite the following reference:
.. HQ X
.. HQ X   http://www.jrkermode.co.uk/quippy
.. HQ X
.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

:class:`quippy.AtomsList` -- manipulate a list of :class:`Atoms` objects
========================================================================

.. currentmodule:: quippy

.. class:: AtomsList(source[, format, lazy, *args, **kwargs])

   An :class:`AtomsList` manages a list of :class:`Atoms` objects. It
   is constructed from the argument `source`. If format is not given,
   and `source` is a filename, then it is inferred from the file
   extension -- see :ref:`fileformats`. If format is given, then
   `source` should be either a filename or an open file, or an
   appropriate file-like object (for example a :class:`CInOutput` for
   the case of XYZ or NetCDF files). Alternatively, `source` can be an
   iterable type containing a sequence of :class:`Atoms` objects, for
   example a generator, tuple or list. If `lazy` given and set to
   ``False`` then all frames in `source` will be read immediately;
   otherwise the default lazyiness is determined by the file type.

   You can access the component objects by indexing and slicing,
   e.g. ``al[i]`` is the ith Atoms object within ``al`` and
   ``al[i:j]`` returns objects i to j inclusive. Like ordinary Python
   lists, indices start from 0 and run up to ``len(al)-1``.

   AtomsList objects support iteration, so you can loop over the
   contents using a :keyword:`for` loop::

      for at in al:
         # process Atoms object `at`
	 pass

   .. note:: 

      The iteration behaviour differs depending on the source from
      which this AtomsList was constructed. If the source supports
      random access then its length will be known in advance, and
      ``al.iteratoms()`` will return an iterator containing exactly
      ``len(al)`` items. If the source does not support random access
      -- for example if it's a generator -- then :meth:`iteratoms`
      will continue to return configurations until the source
      iterator is exhausted. XYZ and NetCDF files both support
      random access (with the exception of reading XYZ from
      standard input). 

      This also affects indexing and slicing: if you try to index
      beyond the end of a AtomList contrusted from a source
      which doesn't support random access, the AtomsList will first
      try to read more frames before raising an :exc:`IndexError`
      exception if insufficent frames are available.


   List comprehensions can also be use, for example to extract
   information from a subset of the Atoms objects within an
   AtomsList::

      min_energy_threshold = -100.0
      print [at.energy for at in al if at.energy > min_energy_threshold]

   AtomsLists provide access to the attributes of their components as
   a single array, with the last array index (slowest varying)
   corresponding to the frame number (plus one, since this is now a
   fortran array). For example the following
   equalities all hold. ::

     al.energy == farray([at.energy for at in al])
     al.energy[0] == al[0].energy         # energy of first frame
     al.velo[0] == al[0].velo             # velocities of all atoms in first frame
     al.velo[:,-1,1] == al[0].velo[:,-1]  # velocity of last atom in first frame

   Methods and attributes:

   .. method:: iteratoms

      Return an interator over the contents of this AtomsList. This is
      the default iterator for an AtomsList instance `al`, and can be
      accessed with ``iter(al)``. See note above regarding iteration
      over sources with do not support random access.

   .. method:: loadall([progress, progress_width, update_interval, show_value])

      Load all frames, equivalent to ``list(al)``. If `progress` is
      present and set to ``True``, a text progress bar is shown to
      indicate the speed of loading. ``progress_width`` (default 80)
      sets the width of this bar in characters, `update_interval` the
      frequency of updates in frames (default is `len(al)` /
      `progress_width`), and `show_value` (default True) controls
      whether the current value should be overlaid on the progress
      bar.


   .. method:: show([property, frame, arrows])
      
      Show this AtomsList using :func:`atomeye.show`. The optional
      arguments `property`, `frame` and `arrows` control the atom
      colouring, initial frame and vector arrows respectively.


   .. method:: write(dest[, format, properties, progress, progress_width, update_interval, show_value, *args, **kwargs])

      Write all frames in this AtomsList to `dest`. If `format` is not
      given it is inferred from the file extension of `dest` (see
      :ref:`fileformats`). If `properties` is present, it should be a list
      of property names to include in the output file, e.g. `['species', 'pos']`.
      
      `progress`, `progress_width`, `update_interval` and
      `show_value` are used to control a textual progress bar, as described
      in :meth:`loadall`. The extra arguments in `*args` and `**kwargs`
      are passed along to the underlying writer routine constructed
      for writing to `dest`.

      See :ref:`fileformats` for a list of supported file formats.

Utility functions
-----------------

.. function:: find_files(filepat[, top])

   Generator which yields a sequence of filenames which match the
   pattern `filepat`, similar to the unix ``find(1)`` command. If
   `top` is absent it defaults to the current working directory.


.. function:: read_files(filenames[, frame, *args, **kwargs])

   Return a generator which reads :class:`Atoms` objects from the
   sequence of files with names given by `filenames`. If `frame` is
   absent the first frame is read from each file, otherwise it
   specifies the frame to load. The extra arguments in `*args` and
   `**kwargs` are passed along to :class:`Atoms.read`. The Atoms objects
   are labelled with a `filename` parameter to indicate the file
   from which they were read.

   This function can be usefully combined with :func:`find_files`::

      filenames = find_files('cluster*.castep')
      source    = read_files(filenames, abort=False)
      al = AtomsList(source)
      al.loadall()


.. function:: AtomsReader(source[, format, *args, **kwargs])

   Simpler generator version of :class:`AtomsList` which does not
   support indexing or slicing, only iteration. All frames in `source`
   are yielded one-by-one as :class:`Atoms` objects.  If `format` is
   given it overrides the automatic guess of the format from the file
   extension of `source`. `*args` and `**kwargs` are passed on to the
   actual reading routine unmodified.
   
