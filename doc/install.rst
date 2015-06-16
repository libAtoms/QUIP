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

Installation
************

The aim of quippy is to make all the functions, subroutines and types
defined in QUIP available from Python scripts. We do this by creating
a Python extension which wraps all of these routines automatically
using `f2py <http://www.scipy.org/F2py>`_, part of 
`numpy <http://numpy.scipy.org>`_. f2py is a very useful tool, but it does
not yet support Fortran 90 extensions such as derived types, so quippy
has to do a little more work to allow us to use them transparently.

Quick start
-----------

For people who don't read manuals::

 $ git clone https://github.com/libAtoms/QUIP.git
 $ cd QUIP
 $ export QUIP_ARCH=linux_x86_64_gfortran
 $ make config
 $ make
 $ make quippy
 $ make install-quippy

Then add the installation directory to your :envvar:`PYTHONPATH`.

Requirements
------------

Essential:
 * `Python 2.6 <http://www.python.org>`_ or later
 * `numpy`_  - version 1.2.1 or later
 * A working fortran compiler:

   * ifort 10 or later (tested with 10.1.015 and 11.0.084)
   * gfortran 4.3.3 or later
   * Others that should work but haven't been tested: `pathf95`, `g95`, `pgf95`, `xlf95`

Optional:
 * `ipython <http://ipython.scipy.org>`_ makes using python interactively 
   much more productive.
 * `matplotlib <http://matplotlib.sourceforge.net>`_ is a useful plotting library which integrates well with ipython
 * `AtomEye <http://mt.seas.upenn.edu/Archive/Graphics/A3/A3.html>`_
   atomistic configuration viewer.  A modified version of AtomEye
   which integrates with quippy is available from my `GitHub
   repository <https://www.github.com/jameskermode/AtomEye>`_
 * `scipy <http://www.scipy.org>`_ provides more scientific
   functionality e.g. least squares fitting, optimisation, etc.

Getting quippy
--------------

You can download QUIP and quippy from our public `GitHub repository
<https://www.github.com/libAtoms/QUIP>`_ ::

  git clone https://github.com/libAtoms/QUIP.git

We'll use the environment variable :envvar:`QUIP_ROOT` to refer
to the root of the QUIP tree::

  export QUIP_ROOT=~/QUIP/

The :file:`libAtoms`, :file:`QUIP_Core`, :file:`QUIP_Util` and
:file:`QUIP_Programs` directories contain Fortran 95 source
code. quippy itself lives in the :file:`{QUIP_ROOT}/Tools/quippy`
directory.

Configuring quippy
------------------

Before compiling quippy, you need to set the environment variable
:envvar:`QUIP_ARCH`. This sets the architecture for the QUIP framework which
quippy wraps. Architectures are defined by creating a file
:file:`${QUIP_ROOT}/Makefiles/Makefile.${QUIP_ARCH}` which describes which
compilers and libraries should be used and where they can be found. quippy has
been tested on the following architectures::

  darwin_x86_64_gfortran
  linux_x86_64_gfortran
  linux_x86_64_ifort_gcc_serial
  linux_x86_64_pgi

If you're on one of these platforms then just set :envvar:`QUIP_ARCH`
appropriately, for example on 64-bit Mac OS X with gfortran you would
do::

  $ export QUIP_ARCH=darwin_x86_64_gfortran

Otherwise you'll have to make a new :file:`Makefile.${QUIP_ARCH}`,
containing some of the variables defined below. The next step is to
answer a number of questions about your system::

  $ make config

This generates a configurtion file `build/${QUIP_ARCH}/Makefile.inc`
which you can subsequently edit if you need to customize any of the
options described below.

Custom settings
^^^^^^^^^^^^^^^

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


Compilation
-----------

It's best to compile QUIP before trying to compile quippy. This will
compile the FoX Fortran XML library as well as generating the required
Makefiles. To compile QUIP, run `make` from the :envvar:`QUIP_ROOT`
directory after setting :envvar:`QUIP_ARCH` appropriately, e.g. ::

  cd ${QUIP_ROOT}
  export QUIP_ARCH=linux_x86_64_gfortran
  make

After this, it's time to compile quippy itself ::

  make quippy
  
to compile quippy. The compilation process is quite long; here is an
overview of the various steps that are performed.

- :mod:`patch_f2py` is invoked to patch the :mod:`numpy.f2py`
  package at runtime to make several changes to the f2py-generated
  C code. This will fail if you don't have :mod:`numpy` 1.2.1 or
  later.

- Fortran sources are imported from the :file:`libAtoms`, :file:`QUIP_Core`, 
  :file:`QUIP_Utils` (if :makevar:`QUIPPY_NO_TOOLS` is not set) 
  and :file:`QUIP_Programs` (if :makevar:`QUIPPY_NO_CRACK` is not set)
  directories. At this stage the sources are preprocessed with the
  :makevar:`QUIPPY_CPP` preprocessor. This removes ``#ifdef`` sections
  so that the tools which read the Fortran source do not get confused
  by multiple version of routines.

- The :mod:`f90doc` module is used to parse Fortran sources and
  analyse all the types, subroutines and functions.

- Using the definitions read by :mod:`f90doc`, the
  :mod:`f2py_wrapper_gen` module writes a Fortran wrapper file for
  each source file that we're going to wrap. These files are named
  :file:`quippy_${STEM}_wrap.f90` and are designed to use the
  restricted set of Fortran 90 features understood by f2py.

- The :file:`quippy_${STEM}_wrap.f90` files are passed to f2py, which 
  generates a Python extension module :mod:`_quippy`. This is a low-level
  module which allows all the Fortran functions to be called from Python,
  but doesn't know anything about derived-types. See :ref:`wrapping-fortran-90-code`
  for more details.

- All the Fortran sources - both those imported and the generated
  wrappers - are compiled using the Fortran compiler specified in
  the :makevar:`QUIPPY_FCOMPILER` Makefile variable. The :mod:`_quippy`
  C extension module is also compiled.

- Finally all the object files are linked, together with external
  libraries such as NetCDF and LAPACK, to create
  :file:`_quippy.so`, the Python extension module. 

If the compilation fails with an error message, please send the full
output to me at j.r.kermode@warwick.ac.uk and I'll do my best to work
out what's going wrong.

Testing
-------

Once quippy is successfully compiled, you should run the test suite to 
check everything is working correctly::

   make test

The tests themselves can be found in :file:`${QUIP_ROOT}/tests/test*.py`.
If any of the tests fail please send me (j.r.kermode@warwick.ac.uk) the output.

Installation
------------

Once all the tests have passed, run ::

   make quippy-install

to install in the standard place for Python extension modules on your
system (this will probably be something like
:file:`/usr/local/lib/python-2.{x}/site-packages`). To install
elsewhere, edit the `Makefile.inc` to include the line::

  QUIPPY_INSTALL_OPTIONS=--prefix=PREFIX

to install to a specific location such as :file:`/usr/local` or your
home directory. You will need to add the installation directory
containing the Python modules to your :envvar:`PYTHONPATH` environment
variable, e.g. by adding the following line your your `~/.bashrc` or
`~/.profile` file::

To test if the installation was successful, change out of the source
directory, start Python and try to import everything from the quippy
package::

   $ cd $HOME
   $ python
   >>> from quippy import *

If you have any problems as this stage such as unresolved linking
dependencies, consult the list of :ref:`install_faq` below.

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





