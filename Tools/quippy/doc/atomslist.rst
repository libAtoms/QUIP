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

:class:`quippy.AtomsReader` and :class:`quippy.AtomsList` -- manipulate trajectories
====================================================================================

.. currentmodule:: quippy

There are two classes for reading trajectories: :class:`AtomsReader`
and :class:`AtomsList`. Use an :class:`AtomsReader` for quick
read-only access to a trajectory or if you only want to access some of
the frames. If you want to load the entire file into memory and
manipulate it use an :class:`AtomsList`.

.. class:: AtomsReader(source[, format, start, stop, step, cache_limit=10, **kwargs])

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

   Methods and attributes:

   .. attribute:: random_access

      True if this source supports random access, False if it does not.

   .. method:: iterframes([reverse])

      Return an interator over all the frames in this trajectory. This
      is the default iterator for an :class:`AtomsReader` instance
      `al`, and can be accessed with ``iter(al)``. 

      If `reverse=True` then the iteration starts with the last frame
      and goes backwards through the file. This is only possible if
      :attr:`random_access` is true.

   .. method:: show([property, frame, arrows])
      
      Visualise using :func:`atomeye.show`. The optional arguments
      `property`, `frame` and `arrows` control the atom colouring,
      initial frame and vector arrows respectively.


   .. method:: write(dest[, format, properties, progress, progress_width, update_interval, show_value, **kwargs])

      Write all frames in this AtomsList to `dest`. If `format` is not
      given it is inferred from the file extension of `dest` (see
      :ref:`fileformats`). If `properties` is present, it should be a list
      of property names to include in the output file, e.g. `['species', 'pos']`.
      
      `progress`, `progress_width`, `update_interval` and `show_value`
      are used to control a textual progress bar. The extra arguments
      in `*args` and `**kwargs` are passed along to the underlying
      writer routine constructed for writing to `dest`.

      See :ref:`fileformats` for a list of supported file formats.


.. class:: AtomsList(source, [format, start, stop, step, *args, **kwargs])

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

   .. method:: sort([cmp, key, reverse, attr])
  
      Sort the AtomsList in place. This is the same as the standard
      :meth:`list.sort` method, except for the additional `attr`
      argument. If this is present then the sorted list will be
      ordered by the :class:`Atoms` attribute `attr`, e.g.::

         al.sort(attr='energy')

      will order the configurations by their `energy` (assuming that
      :attr:`Atoms.params` contains an entry named `energy` for each
      configuration; otherwise an :exc:`AttributError` will be raised).


.. function:: AtomsWriter(dest[, format])

   Returns a file-like object for writing Atoms to `dest` which
   should be either a filename or an initiliased output object.  If
   `format` is not given it is inferred from the file extension of
   `dest`. Example usage::

      out = AtomsWriter('out_file.xyz')
      for at in seq:
         out.write(at)
      out.close()


