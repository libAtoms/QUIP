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

.. _extendedxyz:

Extended XYZ
============

Extended XYZ format is a enhanced version of the `basic XYZ format
<http://en.wikipedia.org/wiki/XYZ_file_format>`_ that allows extra
columns to be present in the file for additonal per atom properties as
well as standardising the format of the comment line to include the
cell lattice and other per frame parameters.

It's easiest to describe the format with an example. Here is a standard XYZ file containing a bulk cubic
8 atom silicon cell ::

  8
  Cubic bulk silicon cell
  Si        0.00000000      0.00000000      0.00000000
  Si        1.36000000      1.36000000      1.36000000
  Si        2.72000000      2.72000000      0.00000000
  Si        4.08000000      4.08000000      1.36000000
  Si        2.72000000      0.00000000      2.72000000
  Si        4.08000000      1.36000000      4.08000000
  Si        0.00000000      2.72000000      2.72000000
  Si        1.36000000      4.08000000      4.08000000

The first line is the number of atoms, followed by a comment and
then one line per atom, giving the element symbol and cartesian
x y, and z coordinates in Angstroms.

Here's the same configuration in extended XYZ format ::

  8
  Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3 Time=0.0
  Si        0.00000000      0.00000000      0.00000000
  Si        1.36000000      1.36000000      1.36000000
  Si        2.72000000      2.72000000      0.00000000
  Si        4.08000000      4.08000000      1.36000000
  Si        2.72000000      0.00000000      2.72000000
  Si        4.08000000      1.36000000      4.08000000
  Si        0.00000000      2.72000000      2.72000000
  Si        1.36000000      4.08000000      4.08000000

In extended XYZ format, the comment line is replaced by a series of
key/value pairs.  The keys should be strings and values can be
integers, reals, logicals (denoted by `T` and `F` for true and false)
or strings. Quotes are required if a value contains any spaces (like
`Lattice` above).  There are two mandatory parameters that any
extended XYZ: `Lattice` and `Properties`. Other parameters --
e.g. `Time` in the example above -- can be added to the parameter line
as needed.

`Lattice` is a Cartesian 3x3 matrix representation of the cell lattice
vectors, in the form ::

  Lattice="R11 R21 R31 R12 R22 R32 R13 R23 R33"

The list of properties in the file is described by the `Properties` parameter, which should take
the form of a series of colon separated triplets giving the name, format (`R` for real, `I` for integer) and number of columns of each property. For example::

  Properties="species:S:1:pos:R:3:vel:R:3:select:I:1"

indicates the first column represents atomic species, the next three
columns represent atomic positions, the next three velcoities, and the
last is an single integer called `select`. With this property
definition, the line ::

  Si        4.08000000      4.08000000      1.36000000   0.00000000      0.00000000      0.00000000       1

would describe a silicon atom at position (4.08,4.08,1.36) with zero velocity and the `select` property set to 1.



