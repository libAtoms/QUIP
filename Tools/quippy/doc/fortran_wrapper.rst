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

.. _wrapping-fortran-90-code:


Appendix: Wrapping Fortran 90 Code
**********************************

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
   prototype for each routine, with derived type arguments replaced by
   integer arrays containing a representation of a pointer to the
   derived type, in the manner described in [Pletzer2008]_ This allows
   us to pass opaque references to the true Fortran derived type data
   structures back and forth between Python and Fortran. The
   :mod:`patch_f2py` module is used to patch the :mod:`numpy.f2py`
   system at runtime, to slightly modify the way it generates C code.
   
3. :mod:`quippy.oo_fortran` acts as a thin object-oriented layer on
   top of the f2py generated wrapper functions which handles
   conversion between Python object instances and Fortran derived-type
   variables, converting arguments back and forth automatically.

Contents:

.. toctree::

   f90doc.rst
   f2py_wrapper_gen.rst
   patch_f2py.rst
   oo_fortran.rst

References:

.. [Pletzer2008] Pletzer, A et al., Exposing Fortran Derived Types to C and Other Languages,
   *Computing in Science and Engineering*, **10**, 86 (2008).
   http://link.aip.org/link/?CSENFA/10/86/1
