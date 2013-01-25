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

:mod:`quippy.farray` --- Fortran one-based array indexing
=========================================================

.. automodule:: quippy.farray
   :synopsis: Fortran one-based array indexing

.. autoclass:: FortranArray
   :members:

There are also some functions which operate on :class:`FortranArray`
instances:

.. autofunction:: farray

.. autofunction:: fenumerate

.. autofunction:: fidentity

.. autofunction:: frange

.. autofunction:: fvar

.. autofunction:: fzeros

.. autofunction:: padded_str_array

.. autofunction:: convert_farray_to_ndarray

.. autofunction:: convert_ndarray_to_farray

.. autofunction:: n2f

.. autofunction:: f2n
