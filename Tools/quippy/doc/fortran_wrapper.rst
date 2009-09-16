.. _wrapping-fortran-90-code:


Wrapping Fortran 90 Code
************************

This section describes the machinery involved in connecting between
Python and Fortran. You should only read it if you want to extend
or modify quippy, or if you're interested in how it works.


There are four steps in the process of wrapping a QUIP Fortran
routine to allow it to be called from quippy. The first two take place
at compile-time (and so these modules do not form part of the
installed quippy package), and the last two at run-time.

1. :mod:`f90doc` is used to scan the Fortran source file,
   building up a data struture which describes all the types, subroutines
   and functions found in the file.
   
2. :mod:`f2py_wrapper_gen` is used to write a simplified Fortran 90
   prototype for each routine, with derived types arguments replaced
   by pointers. These pointers will be treated as integers by the f2py
   tool when it creates a wrapper, allowing us to pass opaque
   references to the true Fortran derived type data structures back
   and forth between Python and Fortran. The :mod:`patch_f2py` module
   is used to patch the :mod:`numpy.f2py` system at runtime, to
   slightly modify the way it generates C code.
   
3. :mod:`quippy.oo_fortran` acts as a thin object-oriented layer on
   top of the f2py generated wrapper functions which handles
   conversion between Python object instances and Fortran derived-type
   variables, converting arguments back and forth automatically.

4. :mod:`quippy.extras` subclasses some of the automatically generated
   classes to create a more Pythonic look and feel. 

Contents:

.. toctree::

   f90doc.rst
   f2py_wrapper_gen.rst
   patch_f2py.rst
   oo_fortran.rst
   extras.rst
