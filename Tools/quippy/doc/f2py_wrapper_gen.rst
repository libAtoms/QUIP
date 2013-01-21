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

:mod:`f2py_wrapper_gen` -- Fortran 90 source code wrapping
==========================================================

.. module:: f2py_wrapper_gen
   :platform: Unix
   :synopsis: Generate simplified Fortran 90 source suitable for f2py

.. moduleauthor:: James Kermode <james.kermode@kcl.ac.uk>

This module is used to write intermediary Fortran 90 wrapper files 
to allow a code which makes use of derived-types to be wrapped with f2py.

All routines which accept derived types arguments are wrapped by
equivalent routines which instead accept integer arrays as opaque
handles.  The Fortran `transfer()` intrinsic is used to convert these
handles into pointers to derived types, as described in
[Pletzer2008]_. Using integer arrays of size 12 makes this approach
portable across most currently available Fortran compilers. In this
way we can access the underlying Fortran structures from Python.

:c:func:`initialise()` and :c:func:`finalise()` routines are handled
specially: on initialisation, a derived type pointer is allocated
before the wrapped routine is invoked, and an opaque reference to this
new derived type is returned. On finalisation the underlying
derived type pointer is deallocated after the wrapped routine returns.

Extra routines are generated to access the values of attributes
within derived types. For scalars a pair of get and set routines is
created, whilst for arrays a single routine which returns the shape,
memory location and type of the array is output. These routines are
used by :mod:`quippy.oo_fortran` when constructing the object-oriented
wrapper layer.

There are various other special cases which are handled individually: for 
details see the 
`source code <http://src.tcm.phy.cam.ac.uk//viewvc/jrk33/repo/trunk/QUIP/Tools/quippy/f2py_wrapper_gen.py?view=markup>`_.

This module defines a single function:

.. function:: wrap_mod(mod, type_map, [out, kindlines])

   :param mod: :class:`f90doc.C_module` instance
   :param type_map: dictionary 
   :param out:  file-like object or None
   :param kindlines: list of strings
   :rtype: dictionary
      
   Write a Fortran 90 wrapper file for the module `mod` to the file `out`. 
   Returns a specification dictionary describing the interface which has been
   created.

   `type_map` should be a dictionary mapping Fortran derived-type
   names to the names of the Fortran modules in which they are defined.

   `kind_lines` is a list of strings to be included near the top of
   the generated module to import any Fortran kind definitions that are needed
   from other modules.

   The `spec` dictionary returned by this function is combined with those
   of other modules and saved to disk by quippy's :file:`setup.py` script.
