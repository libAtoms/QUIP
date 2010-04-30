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

:mod:`quippy.oo_fortran` -- Fortran 90 derived-type support
===========================================================

.. module:: quippy.oo_fortran
   :platform: Unix
   :synopsis: Fortran 90 derived-type support

.. moduleauthor:: James Kermode <james.kermode@kcl.ac.uk>

This module adds support for derived types to f2py generated modules,
providing the Fortran source code was built using
:mod:`f2py_wrapper_gen`.

.. class:: FortranDerivedType

   This class is an abstract base class for all Fortran
   derived-types. It contains all the magic necessary
   to communicate with a module wrapped using the
   :mod:`f2py_wrapper_gen` module.

   The constructor invokes a Fortran interface or routine to initialise
   the object. When the reference count for instances of this class drops to zero,
   the destructor frees the associated Fortran derived type.

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


   .. method:: _runroutine(name, *args, **kwargs)

      Internal method used to invoke the Fortran routine `name` as a
      method. `name` must be a valid key in
      :attr:`_routines`. Wrapper methods which simply call
      :meth:`_runroutine()` are automatically generated in subclasses
      of :class:`FortranDerivedType` by :func:`wrap_all()`.

      Input arguments which are instances of a subclass of :class:`FortranDerivedType`
      are replaced by their :attr:`_fpointer` integer attribute.

      If there is an keyword argument with the name `args_string`
      then unexpected keyword arguments are permitted. All the
      undefined keyword arguments are collected together to form a
      dictionary which is converted to string form and used as the the
      `arg_string` argument, providing rudimentary support for
      variable numbers and types of arguments. For example::

        p = Potential('IP SW', xml_string)
	p.calc(at, calc_virial=True, calc_energy=True)

      is equivalent to::

        p = Potential('IP SW', xml_string)
	p.calc(at, args_str="calc_virial=T calc_energy=T")
	
      The return value us made up of a tuple of the arguments to the
      Fortran routine which are ``intent(out)``. Pointers to Fortran
      derived-type instances are replaced with new instances of the
      appropriate subclass of :class:`FortranDerivedType`. Arrays are
      converted to use one-based indexing using
      :class:`~quippy.farray.FortranArray`.


   .. method:: _runinterface(name, *args, **kwargs)
   
      Internal method used to invoke the appropriate routine
      within the Fortran interface `name`. If no routine is found 
      matching the names and types of the arguments provided then
      an :exc:`TypeError` exception is raised.

      Arguments and results are handled in the same way as 
      :func:`_runroutine`.


   .. method:: _update()

      Automatically invoked whenever this object needs to be
      updated. This happens when it is first created, when it is
      direcly passed to a Fortran routine with ``intent(in,out)`` or
      when it is a component of another object which is passed with
      ``intent(in,out)``.

   .. method:: _update_hook() 

      Invoked by :meth:`_update()`. Can be overriden in subclasses to
      allow a customised response. For example this mechanism is used
      in :class:`quippy.extras.Atoms` to update Atoms properties.


   .. method:: _get_array_shape(name)

      This method can be used to override Fortran's idea of the shape of 
      arrays within derived types, for example to present only a partial
      view of an array. This is used in :class:`quippy.extras.Table`
      to allow the sizes of arrays within the Table class to correspond to the
      current extent of the Table, rather than the size of the allocated storage
      which will usually be larger.

      If this method returns ``None`` then the full array is presented, otherwise
      the return value should be a tuple `(N_1, N_2, ..., N_d)` where `d`
      is the number of dimensions of the array and `N_1, N_2` etc. are the lengths
      of each dimension.


.. function:: wrap_all(topmod, spec, mods, short_names)

   Returns tuple `(classes, routines, params)` suitable for
   importing into top-level package namespace. `topmod` should be an
   f2py-generated module containing `fortran` objects, and `spec`
   should be the specification dictionary generated by
   :func:`f2py_wrapper_gen.wrap_mod`.  `mods` is a list of the names
   of Fortran modules to wrap, and `short_names` is a dictionary
   mapping shortened Fortran derived-type names to their canonical
   form.

   `classes` and `routines` are lists of `(name, value)` tuples
   where `value` is a newly defined subclass of
   :class:`FortranDerivedType` or newly wrapped routine respectively.
   `params` is a Python dictionary of Fortran parameters (constants).

   Here's how this function is used in quippy's :file:`__init.py__` to
   import the new classes, routines and params into the top-level quippy
   namespace::
   
      classes, routines, params = wrap_all(_quippy, spec, spec['wrap_modules'], spec['short_names'])

      for name, cls in classes:
	 setattr(sys.modules[__name__], name, cls)

      for name, routine in routines:
	 setattr(sys.modules[__name__], name, routine)

      sys.modules[__name__].__dict__.update(params)

.. attribute:: FortranDerivedTypes

   Dictionary mapping Fortran type names in format ``type(lower_case_name)`` to
   classes derived from :class:`FortranDerivedType`.
