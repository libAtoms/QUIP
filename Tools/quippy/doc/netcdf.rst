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

NetCDF
======

We use the `NetCDF file format
<http://www.unidata.ucar.edu/software/netcdf/>`_, a flexible binary
file format designed for scientific array data. 

We use a superset of the `AMBER conventions
<http://amber.scripps.edu/netcdf/nctraj.html>`_, so that our
trajectory file can be read directly by VMD. An important distinction
from the Extended XYZ format is the names of some properties:

+---------------------+----------------------+
| Extended XYZ name   | NetCDF name          |
+=====================+======================+
| pos                 | coordinates          |
+---------------------+----------------------+
| velo                | velocities           |
+---------------------+----------------------+

This mapping is handled automatically by, but if you access the data
directly you'll need to be aware of it.

NetCDF versions 3 and 4 are supported. If version 4 is used then it's
possible to use zlib compression, which greatly reduces the file size.

NetCDF Convention
-----------------

All data is either per-atom (referred to as a `property`) or per-frame (refereed to as a `parameter`).

Dimensions
^^^^^^^^^^
 * `frame` - number of frames (this is the unlimited dimension)
 * `spatial` - number of spatial dimensions (i.e. 3)
 * `atom` - number of atoms
 * `cell_spatial` - number of cell lengths (i.e. 3)
 * `cell_angular` - number of cell angles (i.e. 3)
 * `label` - length of string properies (per-atom character data, e.g. species, value is 10)
 * `string` - length of string parameters (per-frame string data, value is 1024)

Variables
^^^^^^^^^

Global variables

  * `spatial (spatial)` - character, equal to `('x','y','z')`
  * `cell_spatial (cell_spatial)` - character, equal to `('a','b','c')`
  * `cell_angular (cell_angular, label)` - character, equal to `('alpha', 'beta', 'gamma')`

Parameters (per-frame variables)
  
 * `cell_lengths (frame, cell_spatial)` - double, cell lengths in Angstrom
 * `cell_angles (frame, cell_angular)` - double, cell angles in degrees

Other parameters can be of type double, integer, logical or
string. Integer, logical and double types can be vectors or scalars
(i.e. dimension `(frame)` or `(frame,spatial)`), but string parameters
must be scalar (i.e. dimension `(frame,string)`. Additionally, real and
integer 3x3 matrices with dimension `(frame,spatial,spatial)` are
supported (e.g. for the virial tensor). In order to distinguish
between integer and logical variables, a special `type` attribute
should be added to the variable, set to one of the following values::

  T_INTEGER = 1
  T_REAL = 2
  T_LOGICAL = 4 
  T_INTEGER_A = 5
  T_REAL_A = 6
  T_LOGICAL_A = 8
  T_CHAR = 9
  T_INTEGER_A2 = 12
  T_REAL_A2 = 13

Properties (per-atom variables)

Properties can be of type integer, real, string or logical. As for parameters,
integer, real and logical properties can be scalar
(dimension `(frame,atom)`) or vector (dimension `(frame,atom,spatial)`), but
string properties must be of dimension `(frame,atom,label)`. Again a `type`
attribute must be added to the NetCDF variable, with one of the following values::

  PROPERTY_INT     = 1
  PROPERTY_REAL    = 2
  PROPERTY_STR     = 3
  PROPERTY_LOGICAL = 4

See the :mod:`xyz_netcdf` module for a reference implementation of the
NetCDF reading and writing routines, in pure Python. The standard
implementation is in C in the :file:`xyz_netcdf.c` file within the
libAtoms directory. It is the version which is wrapped by the
:class:`CInOutput` class.
