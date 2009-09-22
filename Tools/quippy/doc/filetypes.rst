.. _fileformats:

Supported File Formats
======================

.. currentmodule:: quippy

The :attr:`AtomsReaders` and :attr:`AtomsWriters` dictionaries are
used by the :class:`Atoms` constructor, :meth:`Atoms.write`, the
:class:`AtomsList` constructor and :meth:`AtomsList.write` to work out
how to read or write a particular file type, based on the filename
extension.

The quippy native formats are `extended XYZ
<http://www.jrkermode.co.uk/ExtendedXYZ>`_ and `NetCDF
<https://camtools.cam.ac.uk/access/wiki/site/5b59f819-0806-4a4d-0046-bcad6b9ac70f/netcdf.html>`_.

.. attribute:: AtomsReaders

   Supported file formats for reading :class:`Atoms` objects from files.

   +----------------+-----------------------+
   | File extension |  Description          |
   +================+=======================+
   | ``castep`` or  | CASTEP output file    |
   | ``castep_log`` |                       |
   +----------------+-----------------------+
   | ``cell``       | CASTEP cell file      |
   +----------------+-----------------------+
   | ``geom``       | CASTEP geom file      |
   +----------------+-----------------------+
   | ``nc``         | NetCDF  file          |
   +----------------+-----------------------+
   | ``pos``        | ASAP coordinate file  |
   +----------------+-----------------------+
   | ``stdin``      | Read from stdin in    |
   |                | extended XYZ format   |
   +----------------+-----------------------+
   | ``xyz``        | Extended XYZ format   |
   +----------------+-----------------------+

.. attribute:: AtomsWriters

   Supported file formats for writing :class:`Atoms` objects to files.

   +----------------+-----------------------+
   | File extension |  Description          |
   +================+=======================+
   | ``cell``       | CASTEP cell file      |
   +----------------+-----------------------+
   | ``nc``         | NetCDF  file          |
   +----------------+-----------------------+
   | ``pos``        | ASAP coordinate file  |
   +----------------+-----------------------+
   | ``pov``        | POV-ray script        |
   +----------------+-----------------------+
   | ``stdout``     | Write to stdout in    |
   |                | extended XYZ format   |
   +----------------+-----------------------+
   | ``xyz``        | Extended XYZ format   |
   +----------------+-----------------------+

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
:class:`Atoms` objects, raising :exc:`StopIteration` when there are no
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

See the code in :mod:`quippy.xyz_netcdf` and :mod:`quippy.castep` for
full examples.


