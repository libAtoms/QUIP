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

.. _installation:

Installation of QUIP and quippy
*******************************

These intructions provide more details on the compilation and
installation of `QUIP` (Fortran library and main programs) and
`quippy` (Python interface).

Please start by following the instructions in the main `README
<https://github.com/libAtoms/QUIP/blob/public/README.md>`_, which is
the most up-to-date source of information.

Custom settings
---------------

:makevar:`MATHS_LINKOPTS`
   Library options needed to link to BLAS and LAPACK libraries. Any working
   BLAS/LAPACK installation is fine. If you are using Linux, ATLAS is
   a good option, and you should use something like the following::

     -L/usr/local/atlas -llapack -lf77blas -lcblas -latlas

   On Mac OS X, there are build in LAPACK libraries in the Accelerate
   framework, which you can use by entering

     -framework Accelerate

:makevar:`FOX_LIBDIR`, :makevar:`FOX_INCDIR` and :makevar:`FOX_LIBS`
  Directories containing FoX libraries and header files, and required link options.
  Should be read automatically from QUIP Makefiles.

:makevar:`QUIPPY_FCOMPILER`
   Fortran compiler to use. The shell command::

     $ f2py -c --help-fcompiler

   will print a list of detected compilers on your system. Use ``gnu95`` for gfortran,
   ``intel`` for ifort on 32-bit platforms and ``intelem`` for ifort on 64-bit platforms.

:makevar:`QUIPPY_DEFINES` Preprocessor macros which should be defined
   when compiling quippy. Note that since the Fortran source files are
   preprocessed *before* being scanned by :mod:`f90doc`, it's
   important to put all the `-D` options needed here and not in
   :makevar:`QUIPPY_F90FLAGS`.


:makevar:`QUIPPY_F90FLAGS` and :makevar:`QUIPPY_F77FLAGS`
   Extra flags to pass to Fortran 90 and 77 compilers

:makevar:`QUIPPY_OPT`
   Optimisation settings for Fortran compiler

:makevar:`QUIPPY_DEBUG`
   Set this to `1` to include debugging information in the compiled extension code.
   This also disables optimisation.

:makevar:`QUIPPY_CPP`
   Fortran preprocessor to use. Default is system `cpp`.

:makevar:`QUIPPY_INSTALL_OPTS`
   Installation options, e.g. specify ``--home=${HOME}``
   or ``--prefix=${PREFIX}`` to install in a non-default location.

:makevar:`QUIPPY_NO_TOOLS`
   If set to 1, omit compilation of extra tools such as the elasticity module.

:makevar:`QUIPPY_NO_CRACK`
  If set to 1, omit compilation of crack utilities.

:makevar:`HAVE_NETCDF`
  Should be set to 1 to enable NetCDF support. Should be read automatically from QUIP.

:makevar:`NETCDF4`
  If set to 1, use version 4 of NetCDF. Should be read automatically from QUIP.

:makevar:`NETCDF_LIBDIR`, :makevar:`NETCDF_INCDIR`, :makevar:`NETCDF_LIBS` and :makevar:`NETCDF4_LIBS`
  Directories containing NetCDF libraries and header files, and required link options.
  Should be read automatically from QUIP.


.. _install_faq:

Common Problems
---------------

Permission errors when installing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are installing as root, you may need to make sure the value of
the :envvar:`QUIP_ARCH` gets through to the install script, e.g. ::

   sudo QUIP_ARCH=darwin_x86_64_gfortran make install-quippy


Installating on Mac OS X with macports
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Macports requires various packages to be installed to compile
everything, and may require extra linking arguments. See the
:file:`README.macports` for the latest details.

RuntimeError when importing
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If, after installing quippy, you get the error shown below when you
try to import it for the first time, then you are a victim of a bug in
early versions of Python 2.6.

::

   >>> import quippy
   Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home/ab686/QUIP/Tools/quippy/quippy/__init__.py", line 31, in
   <module>
      _quippy.system.verbosity_push(0)
   RuntimeError: more argument specifiers than keyword list entries
   (remaining format:'|:_quippy.system.verbosity_push')

The solution is either to compile your own Python from the current svn
snapshot, or to update numpy to workaround the fix. This can be done
either by compiling numpy from source from an up-to-date svn snapshot,
or by applying `the patch manually
<http://projects.scipy.org/numpy/changeset/6193>`_.

ImportError when importing
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get an :exc:`ImportError` with a message about unresolved
dependancies then something went wrong with the linking process -
check that all the libraries you're linking against are correct. You
can used `ldd` on Linux of `otool -L` on Mac OS X to check which
libraries the :file:`_quippy.so` Python extension is linked against.

Possible problems installing atomeye module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get an :exc:`ImportError` with a message ::
   >>> import atomeye
   ImportError: dlopen(/Users/silvia/lib/python/_atomeye.so, 2): Symbol not found: _Config_load_libatoms
   Referenced from: /Users/silvia/lib/python/_atomeye.so
   Expected in: flat namespace
   in /Users/silvia/lib/python/_atomeye.so

be sure that you have set :envvar:`QUIP_ROOT` variable before starting the compilation.
If not make clean and recompile again

If you get an :exc:`ImportError` with a message ::
   >>> import atomeye
   ImportError: dlopen(/Users/silvia/lib/python/_atomeye.so, 2): Symbol not found: __gfortran_adjustl
   Referenced from: /Users/silvia/lib/python/_atomeye.so
   Expected in: flat namespace
   in /Users/silvia/lib/python/_atomeye.so

be sure that the gfortran libraries are properly set in :makevar:`ATOMEYE_LIBS` in Makefile.atomeye

Error compiling IPModel_GAP
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get the following error during compilation::

   /QUIP/QUIP_Core/IPModel_GAP.f95:51.22:

   use descriptors_module
                         1
   Fatal Error: Can't open module file 'descriptors_module.mod' for reading at (1): No such file or directory

The `GAP_predict` module is not publicly available, so the
:file:`Makefile.inc` must contain :makevar:`HAVE_GP_PREDICT` = 0, and
:makevar:`HAVE_GP_TEACH` = 0.


Warning about :mod:`quippy.castep` when importing quippy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get the following warning message when importin quippy::

   $ python
   >>> from quippy import *
   WARNING:root:quippy.castep import quippy.failed.

then don't worry, the quippy.castep module is not redistributed with
the main code. The rest of quippy works fine without it.


Internal compiler error with `ifort`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you see an error like the following when using the Intel fortran compiler::

   fortcom: Severe: **Internal compiler error: internal abort** Please
   report this error along with the circumstances in which it occurred
   in a Software Problem Report.
    Note: File and line given may not be explicit cause of this error.

   ifort: error #10014: problem during multi-file optimization compilation (code 3)
   backend signals

Then the problem is due to bugs in the compiler. As a workaround,
setting :makevar:`QUIPPY_NO_CRACK` =1 in Makefile.inc should solve the
problem, at the cost of excluding the fracture utilities from quippy.


Linking error on Mac OS X
^^^^^^^^^^^^^^^^^^^^^^^^^

When recompiling quippy on top of a previous compilation, you may see
errors like this::

   collect2: ld returned 1 exit status ld: in
   /QUIP/build.darwin_x86_64_gfortran/libquiputils.a, malformed
   archive TOC entry for  ___elasticity_module_MOD_einstein_frequencies,
   offset 1769103734 is beyond end of file 1254096

This seems to be a Mac OS X Lion problem with rebuilding static
libraries (.a files). Removing the static libraries with `rm
../../build.${QUIP_ARCH}/*.a` and recompiling should solve the
problem.





