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

:mod:`quippy.oo_fortran` --- Fortran 90 derived-type support
============================================================

.. automodule:: quippy.oo_fortran
   :synopsis: Fortran 90 derived-type support


.. autoclass:: FortranDerivedType
   :members:

   Important attributes and methods are as follows:

   .. attribute:: _fpointer

      Pointer to underlying Fortran derived type instrance, or ``None`` if
      not initialised.

   .. attribute:: _routines

      Dictionary mapping routine names to tuples `(fortran_object, specification)`.
      Used by :meth:`_runroutine()`.

   .. attribute:: _interfaces

      Dictionary mapping interface names to a list of 3-tuples `(name, specification, fortran_object)`.
      Used by :meth:`_runinterface()`.

   .. attribute:: _elements

      Dictionary mapping Fortran element names to tuple `(get_routine, set_routine)`. 
      A property is created for each entry in `_elements` which invokes `get_routine`
      and `set_routine` automatically.
      
   .. attribute:: _subobjs

      Dictionary mapping names of component objects to tuples `(get_routine, set_routine)`.
      Used by :meth:`_update()`.

   .. attribute:: _arrays

      Dictionary mapping names of arrays to tuples `(array_func, doc_string)`. Used by
      :meth:`_update()`.


   .. automethod:: _runroutine(name, *args, **kwargs)


   .. automethod:: _runinterface(name, *args, **kwargs)
   
   .. automethod:: _update()

   .. automethod:: _get_array_shape(name)


.. autofunction:: wrap_all(topmod, spec, mods, short_names)


.. attribute:: FortranDerivedTypes

   Dictionary mapping Fortran type names in format ``type(lower_case_name)`` to
   classes derived from :class:`FortranDerivedType`.
