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

AtomsList and AtomsReader objects for I/O
=========================================

.. automodule:: quippy.io
   :synopsis: Read and write atomic configurations
   :members:

.. _fileformats:

Supported File Formats
----------------------

The :attr:`AtomsReaders` and :attr:`AtomsWriters` dictionaries are
used by the :class:`Atoms` constructor, :meth:`Atoms.write`, the
:class:`AtomsList` constructor and :meth:`AtomsList.write` to work out
how to read or write a particular file type, based on the filename
extension.

The quippy native formats are :ref:`extendedxyz` and :ref:`netcdf`.

The standard implementation of both these formats is in C in the files
:file:`xyz.c` and :file:`netcdf.c` in the `libAtoms` package. It is
this version which is wrapped by the
:class:`~quippy.cinoutput.CInOutput` class, which is used by
:class:`~quippy.io.AtomsReader` and :class:`~quippy.io.AtomsList` when
reading from or writing to XYZ or NetCDF files.

.. attribute:: AtomsReaders

   Supported file formats for reading :class:`~quippy.atoms.Atoms` objects from files.

   +-----------------+-----------------------------+
   | File extension  |  Description                |
   +=================+=============================+   
   | ``castep`` or   | :ref:`castep` output        |
   | ``castep_log``  |                             |
   +-----------------+-----------------------------+
   | ``cell``        | :ref:`castep` cell files    |
   +-----------------+-----------------------------+
   | ``chkpt``       | :ref:`imd`                  |
   +-----------------+-----------------------------+
   | ``cp2k_output`` | :ref:`cp2k`                 |
   +-----------------+-----------------------------+
   | ``cube``        | :ref:`cube`                 |
   +-----------------+-----------------------------+
   | ``geom``        | :ref:`castep` geometry      |
   +-----------------+-----------------------------+
   | ``md``          | :ref:`castep` MD file       |
   +-----------------+-----------------------------+
   | ``nc``          | :ref:`netcdf`               |
   +-----------------+-----------------------------+
   | ``pos``         | :ref:`asap`                 |
   +-----------------+-----------------------------+
   | ``POSCAR`` or   | :ref:`vasp` coordinates     |
   | ``CONTCAR``,    |                             |
   +-----------------+-----------------------------+
   | ``OUTCAR``,     | :ref:`vasp` output          |
   +-----------------+-----------------------------+
   | ``stdin``       | Read from stdin in          |
   |                 | :ref:`extendedxyz` format   |
   +-----------------+-----------------------------+
   | ``string``      | Read from string in         |
   |                 | :ref:`extendedxyz` format   |
   +-----------------+-----------------------------+
   | ``xyz``         | :ref:`extendedxyz`          |
   +-----------------+-----------------------------+

.. attribute:: AtomsWriters

   Supported file formats for writing :class:`~quippy.atoms.Atoms` objects to files.

   +----------------+----------------------------+
   | File extension |  Description               |
   +================+============================+
   | ``cell``       | :ref:`castep` cell file    |
   +----------------+----------------------------+
   | ``cube``       | :ref:`cube`                |
   +----------------+----------------------------+
   | ``dan``        | :ref:`dan`                 |
   +----------------+----------------------------+
   | ``eps``,       | Images                     |
   | ``jpg``,       | (via:ref:`atomeyewriter`)  |
   | ``png``,       | 			         |
   +----------------+----------------------------+
   | ``nc``         | :ref:`netcdf`              |
   +----------------+----------------------------+
   | ``pos``        | :ref:`asap`                |
   +----------------+----------------------------+
   | ``pov``        | :ref:`povray`              |
   +----------------+----------------------------+
   | ``POSCAR``     | :ref:`vasp` coordinates    |
   +----------------+----------------------------+
   | ``-``,         | Write to stdout in         |
   | ``stdout``     | :ref:`extendedxyz` format  |
   +----------------+----------------------------+
   | ``string``     | Write to string in         |
   |                | :ref:`extendedxyz` format  |
   +----------------+----------------------------+
   | ``xyz``        | :ref:`extendedxyz`         |
   +----------------+----------------------------+


.. _extendedxyz:

Extended XYZ
------------

.. module:: quippy.xyz
  :synopsis: Python reference implementation of Extended XYZ I/O

Extended XYZ format is a enhanced version of the `basic XYZ format
<http://en.wikipedia.org/wiki/XYZ_file_format>`_ that allows extra
columns to be present in the file for additonal per atom properties as
well as standardising the format of the comment line to include the
cell lattice and other per frame parameters.

It's  easiest to  describe  the  format with  an  example.  Here is  a
standard XYZ file containing a bulk cubic 8 atom silicon cell ::

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
e.g. `Time` in the example above --- can be added to the parameter line
as needed.

`Lattice` is a Cartesian 3x3 matrix representation of the cell lattice
vectors in Fortran column-major order, i.e. in the form ::

  Lattice="R1x R1y R1z R2x R2y R2z R3x R3y R3z"

where `R1x R1y R1z` are the Cartesian x-, y- and z-components of the
first lattice vector (:math:`\mathbf{a}`), `R2x R2y R2z` those of the second
lattice vector (:math:`\mathbf{b}`) and `R3x R3y R3z` those of the
third lattice vector (:math:`\mathbf{c}`).

The list of properties in the file is described by the `Properties`
parameter, which should take the form of a series of colon separated
triplets giving the name, format (`R` for real, `I` for integer) and
number of columns of each property. For example::

  Properties="species:S:1:pos:R:3:vel:R:3:select:I:1"

indicates the first column represents atomic species, the next three
columns represent atomic positions, the next three velcoities, and the
last is an single integer called `select`. With this property
definition, the line ::

  Si        4.08000000      4.08000000      1.36000000   0.00000000      0.00000000      0.00000000       1

would describe a silicon atom at position (4.08,4.08,1.36) with zero
velocity and the `select` property set to 1.

The extended XYZ format is now also supported by the
:func:`ase.io.read` and :func:`ase.io.write` functions in the `Atomic
Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase>`_
toolkit, and by the `Ovito <http://www.ovito.org>`_ visualisation tool
(from `v2.4 beta
<http://www.ovito.org/index.php/component/content/article?id=25>`_
onwards).

.. _netcdf:

NetCDF
------

.. module:: quippy.netcdf
  :synopsis: Python reference implementation of NetCDF I/O

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
^^^^^^^^^^^^^^^^^

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

See the :mod:`quippy.netcdf` module for a reference implementation of the
NetCDF reading and writing routines, in pure Python.

CInOutput objects
-----------------

.. automodule:: quippy.cinoutput
   :synopsis: XYZ and NetCDF I/O using fast C routines
   :members:
   :inherited-members:
   :show-inheritance:

.. _asap:

ASAP file format
----------------

.. automodule:: quippy.asap
   :synopsis: ASAP file reader and writer
   :members:

.. _atomeyewriter:

AtomEye Image Writer
--------------------

.. automodule:: quippy.atomeyewriter
   :synopsis: AtomEye image writer
   :members:

.. _castep:

CASTEP
------

.. automodule:: quippy.castep
   :synopsis: CASTEP file readers and writers
   :members:



.. _cube:

Gaussian CUBE
-------------

.. automodule:: quippy.cube
   :synopsis: Gaussian CUBE I/O
   :members:

.. _dan:

DAN visualisation code
----------------------

.. automodule:: quippy.dan
   :synopsis: DAN visualisation code writer
   :members:

.. _imd:

IMD checkpoint
--------------

.. automodule:: quippy.imd
   :synopsis: IMD checkpoint reader
   :members:

.. _povray:

POV-ray
-------

.. automodule:: quippy.povray
   :synopsis: POV-ray script writer
   :members:

.. _vasp:

VASP
----

.. automodule:: quippy.vasp
   :synopsis: VASP POSCAR and OUTCAR readers and writers
   :members:


ASE supported files types
-------------------------

Since :class:`quippy.atoms.Atoms` is a subclass of the ASE
:class:`ase.atoms.Atoms` class, all of the ASE I/O formats can also be
used with quippy: see :func:`ase.io.read` for a list of supported
formats. To convert from :class:`ase.atoms.Atoms` to
:class:`quippy.atoms,Atoms`, simply pass the ASE Atoms object to the
quippy Atoms constructor, e.g.::

   from quippy.atoms import Atoms as QuippyAtoms
   from ase.atoms import Atoms as ASEAtoms
   from ase.io import read

   ase_atoms = read(filename)
   quippy_atoms = QuippyAtoms(ase_atoms)

Similarily, to use one of the quippy file formats with other ASE
tools::

   from quippy.io import read
   quippy_atoms = read(filename)
   ase_atoms = ASEAtoms(quippy_atoms)


Adding a new file type
----------------------

To add support for a new file format, implement routines which
read from or write to files following the templates below. ::

     def sample_reader(filename): 
	 # insert code to open `filename` for reading

	 while True:
	    # determine if more frame are available
	    if more_frames:
	       # read next frame from `filename` into new Atoms object
	       at = Atoms()
	       yield at
	    else:
	       break
	   

     class sample_writer(object):
	 def __init__(self, filename):
	    # insert code to open `filename` for writing
	    pass
	 
	 def write(self, at):
            # insert code to write `at` to `filename`
	    pass	    

         def close(self):
	    # insert code to close `filename`
	    pass

:func:`sample_reader` is a generator which yields a succession of
:class:`~quippy.atoms.Atoms` objects, raising :exc:`StopIteration` when there are no
more available - for a file format which only permits one
configuration per file, a simplified implementation template would be::

   def sample_reader(filename):
      # insert code to open `filename` for reading
      # insert code to read from `filename` into new Atoms object
      yield at
      

To register the new file format, you just need to set entries in
:attr:`AtomsReaders` and :attr:`AtomsWriters`::
     
     from quippy import AtomsReaders, AtomsWriters
     AtomsReaders['new_format'] = sample_reader
     AtomsWriters['new_format'] = sameple_writer

For the case of reading, there is a generator :func:`atoms_reader`
which can be used to simplify the registration process::

   @atoms_reader('new_format')
   def sample_reader(filename):
      ...

See the code in :mod:`quippy.xyz` :mod:`quippy.netcdf` and
:mod:`quippy.castep` for full examples.
