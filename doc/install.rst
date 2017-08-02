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
installation of ``QUIP`` (Fortran library and main programs) and
``quippy`` (Python interface).

Compilation Instructions
------------------------

First try the quickstart below, which should work with most Linux systems.
For Mac systems, have a look at `Installing on Mac OS X with macports`_ first.

Quick start
^^^^^^^^^^^

Install [#]_ the prerequisites: GCC, gfortran, Python, and the linear algebra
libraries.  For example, on Ubuntu, do (in a terminal):

::

    $ sudo apt-get install gcc gfortran python python-pip libblas-dev liblapack-dev

For other systems, replace the ``apt-get`` part with your system package manager.
Beware that the packages might also have slightly different names; these can
usually be found with a quick search.

Don't forget the ``quippy`` prerequisites:

::

    $ pip install numpy
    $ pip install ase


Now you can get the code and compile:

::

    $ git clone --recursive https://github.com/libAtoms/QUIP.git
    $ export QUIP_ARCH=linux_x86_64_gfortran
    $ export QUIPPY_INSTALL_OPTS=--user  # omit for a system-wide installation
    $ make config

Answer all the questions with their defaults (by pressing enter) for now, just
to get things working.

::

    $ make
    $ make install-quippy

And now open a Python terminal and see if it works:

::

    $ python
    >>> import quippy
    >>>

If the import completes successfully (i.e. with no output) then the
installation was successful.  You may want to continue with `Installing the
Jupyter notebook`_ to run the interactive tutorials.

.. [#] If this isn't your machine and you don't have root access, these
   packages might already be installed by the system administrator.  If not,
   ask them.


Step by step
^^^^^^^^^^^^

If that didn't work, try these step-by-step instructions
instructions excerpted from the top-level `README
<https://github.com/libAtoms/QUIP/blob/public/README.md>`_.  The ``README`` file
is the most up-to-date source of installation information.

  .. include:: ../README.md
    :start-after: Compilation Instructions

If that still doesn't work or you're using a nonstandard architecture, try
looking at `Custom settings`_ and `Common Problems`_.  As a last resort you can
consult the `issue tracker on Github`_.


Installing the Jupyter notebook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Jupyter`_ is an environment for interactive computing that makes using Python
much easier and more intuitive.  Especially useful is its notebook environment,
which provides a handy way to experiment with code, see the results, and have a
record of your progress.  The interactive getting-started tutorial is a Jupyter
notebook that you can run and modify yourself.

To get Jupyter up and running, the following should suffice [#]_:

::

    $ pip install jupyter
    $ jupyter notebook

This will open a new window in your browser that you can use to navigate
through your filesystem.  To access the interactive tutorials, you can run the
``jupyter notebook`` command from your ``QUIP/doc/Examples`` directory (or any
enclosing directory) then navigate to the notebooks and open
``Introduction.ipynb`` to get started.

.. [#] This assumes you've already run ``sudo apt-get install python-pip; pip
   install numpy; pip install ase`` as in the `Quick start`_.


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
   Installation options, e.g. specify ``--user`` to install for the current
   user ``--prefix=${PREFIX}`` to install in a non-default location.

:makevar:`QUIPPY_NO_TOOLS`
   If set to 1, omit compilation of extra tools such as the elasticity module.

:makevar:`QUIPPY_NO_CRACK`
  If set to 1, omit compilation of crack utilities.

:makevar:`HAVE_NETCDF4`
  Should be set to 1 to enable NetCDF4 support. Should be read automatically from QUIP.

:makevar:`NETCDF4_LIBS`, :makevar:`NETCDF4_FLAGS`
  Linker flags for compiling with NetCDF4 support, and flags for finding
  header files. Should be read automatically from QUIP.


.. _install_faq:

Common Problems
---------------

Permission errors when installing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are installing as root, you may need to make sure the value of
the :envvar:`QUIP_ARCH` gets through to the install script, e.g. ::

   sudo QUIP_ARCH=darwin_x86_64_gfortran make install-quippy


Installing on Mac OS X with macports
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

   /src/Potentials/IPModel_GAP.f95:51.22:

   use descriptors_module
                         1
   Fatal Error: Can't open module file 'descriptors_module.mod' for reading at (1): No such file or directory

The `GAP_predict` module is not publicly available, so the
:file:`Makefile.inc` must contain :makevar:`HAVE_GP_PREDICT` = 0, and
:makevar:`HAVE_GP_TEACH` = 0.


Warning about :mod:`quippy.castep` when importing quippy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get the following warning message when importing quippy::

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


Segmentation Faults with OpenBLAS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The threading in OpenBLAS can interfere with the OpenMP resulting in
segfaults. Either recompile OpenBLAS with ``USE_OPENMP=1`` or disable
threading with ``export OPENBLAS_NUM_THREADS=1`` at runtime.


.. _`issue tracker on Github`: https://github.com/libAtoms/QUIP/issues
.. _`Jupyter`: http://jupyter.org/

