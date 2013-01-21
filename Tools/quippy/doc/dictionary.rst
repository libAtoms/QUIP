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

:class:`quippy.dictionary.Dictionary` -- Fortran dictionary for key/value pairs
===============================================================================

.. currentmodule:: quippy.dictionary

.. class:: Dictionary([D])
   
   The :class:`Dictionary` class is designed to behave as much as
   possible like a true Python dictionary, but since it is implemented
   in Fortran it can only store a restricted range of data types.
   Keys much be strings and values must be one of the following types:

   - Integer
   - Real
   - String
   - Complex
   - Logical
   - 1D integer array 
   - 1D real array 
   - 1D complex array
   - 1D logical array
   - 2D integer array
   - 2D real array

   Trying to store any other type of data will raise a :exc:`ValueError`.

   For :attr:`Atoms.params` entries, there are further restrictions
   imposed by the implementation of the XYZ and NetCDF I/O
   routines. The only types of data that can be stored here are:

   - Integer
   - Real
   - String
   - Integer 3-vector
   - Real 3-vector
   - Integer 3 x 3 matrix
   - Real 3 x 3 matrix
   
   A :class:`Dictionary` can be created from a standard Python dictionary,
   and easily converted back::

      >>> py_dict = {'a':1, 'b':2}
      >>> fortran_dict = Dictionary(py_dict)
      >>> py_dict == dict(fortran_dict)
      True
   
   It also supports all the standard dictionary operations and methods::

      >>> fortran_dict['c'] = 3
      >>> fortran_dict.keys()
      ['a', 'b', 'c']


   An additional feature of the quippy :class:`Dictionary` is that it can 
   read and write itself to a string in the format used within XYZ files::

      >>> str(fortran_dict)
      'a=1 b=2 c=3'
      >>> d2 = Dictionary('a=1 b=2 c=3')
      >>> d2.keys(), d2.values()      
      (['a', 'b', 'c'], [1, 2, 3])

   
