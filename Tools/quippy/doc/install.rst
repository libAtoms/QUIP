.. _installation:

Installation
************

The aim of quippy is to make all the functions, subroutines and types
defined in QUIP available from Python scripts. We do this by creating
a Python extension which wraps all of these routines automatically
using `f2py <http://www.scipy.org/F2py>`_, part of 
`numpy <http://www.scipy.org/Download>`_. f2py is a very useful tool, but it does
not yet support Fortran 90 extensions such as derived types, so quippy
has to do a little more work to allow us to use them transparently.

Requirements
------------

Essential:
 * `Python 2.4 <http://www.python.org>`_ or later
 * `numpy`_  - version 1.2.1 or later
 * A fortran compiler:

   * ifort 10 or later (tested with 10.1-015)
   * gfortran, svn version of 4.3.3 branch, available by
     anonymous svn from `svn://gcc.gnu.org/svn/gcc/tags/gcc_4_3_3_release 
     <svn://gcc.gnu.org/svn/gcc/tags/gcc_4_3_3_release>`_.
     With 4.3.2 or earlier you run into an 
     `internal compiler error <http://gcc.gnu.org/bugzilla/show_bug.cgi?id=37735>`_, 
     reported by Steve Winfield and now fixed.
   * (PGI Fortran 90 - almost works)

Optional:
 * `ipython <http://ipython.scipy.org>`_ makes using python interactively 
   much more productive.
 * `matplotlib <http://matplotlib.sourceforge.net>`_ is a useful plotting library which integrates well with ipython
 * `AtomEye <http://mt.seas.upenn.edu/Archive/Graphics/A3/A3.html>`_
   atomistic configuration viewer.  A modified version of AtomEye
   which integrates with quippy is available from the QUIP `svn
   repository <http://src.tcm.phy.cam.ac.uk/viewvc/jrk33/repo/trunk/AtomEye>`_

Getting quippy
--------------

If you have an account on the QUIP `svn repository server
<https://camtools.cam.ac.uk/access/wiki/site/5b59f819-0806-4a4d-0046-bcad6b9ac70f/svnrepository.html>`_, 
you should checkout a copy of the QUIP trunk::

  svn checkout svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/QUIP QUIP

If you don't have an account, you can download a `snapshot
<svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/QUIP?view=tar>`_
of the svn repository as a tarball. Extracting this archive will create a
:file:`QUIP` tree.

We'll use the environment variable :envvar:`QUIP_ROOT` to refer
to the root of the QUIP tree. The :file:`libAtoms`, :file:`QUIP_Core`,
:file:`QUIP_Util` and :file:`QUIP_Programs` directories contain
Fortran 95 source code. quippy itself lives in the
:file:`{QUIP_ROOT}/Tools/quippy` directory

If you want to include the AtomEye extension, you should check it out
under :file:`${QUIP_ROOT}/Tools`::

  cd ${QUIP_ROOT}/Tools
  svn checkout svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/AtomEye AtomEye

Again, if you don't have an account you can download a `snapshot
<svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/QUIP?view=tar>`_, 
which you should extract under :file:`${QUIP_ROOT}/Tools`


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
  linux_x86_64_ifort_gcc
  linux_x86_64_pgi

If you're on one of these platforms then just set :envvar:`QUIP_ARCH`
appropriately, for example on 64-bit Mac OS X with gfortran you would
do::

  $ export QUIP_ARCH=darwin_x86_64_gfortan

Otherwise you'll have to make a new :file:`Makefile.${QUIP_ARCH}`. It
should define various quippy-specific variables:

:makevar:`QUIPPY_FCOMPILER`
   Fortran compiler to use. The shell command::

     $ f2py -c --help-fcompiler 

   will print a list of detected compilers on your system. Use ``gnu95`` for gfortran, 
   ``intel`` for ifort on 32-bit platforms and ``intelem`` for ifort on 64-bit platforms.

:makevar:`QUIPPY_INSTALL_OPTS`
   Installation options, e.g. specify ``--home=${HOME}``
   or ``--prefix=${PREFIX}`` to install in a non-default location.

:makevar:`QUIPPY_F90FLAGS`
   Flags to pass to Fortran compiler

:makevar:`QUIPPY_OPT`
   Optimisation settings for Fortran compiler

:makevar:`QUIPPY_CPP`
   Fortran preprocessor to use 

:makevar:`QUIPPY_NO_TOOLS`
   If set to 1, omit compilation of extra tools such as the elasticity module.

:makevar:`QUIPPY_NO_CRACK`
  If set to 1, omit compilation of crack utilities. Currently this is
  necessary with ``ifort`` to avoid an internal compiler error. 

:makevar:`QUIPPY_HAVE_ATOMEYE`
  Set this to 1 if you want to build the AtomEye interface module.

Compilation
-----------

After all the Makefile variables desribed above have been setup, run
``make install`` to compile and install everything. The process is
quite long; here is an overview of the various steps that are
performed.

* ``make`` invokes :file:`presetup.py` which compiles a small C program
  to determine the size of a ``void*`` pointer on your architecture.

* If AtomEye support is enabled, ``make`` compiles the AtomEye C source 
  as a library, suitable for linking into the :mod:`_atomeye` Python
  C extension module.

* ``make`` invokes :file:`setup.py` which does the rest of the work:

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
     by multiple version of 

   - The :mod:`f90doc` module is used to parse Fortran sources and
     analyse all the subroutines and functions defined. Only the files
     listed in the :makevar:`WRAP_SOURCES_LIST` Makefile variable are
     looked at.

   - Using the definitions read by :mod:`f90doc`, the
     :mod:`f2py_wrapper_gen` module writes a Fortran wrapper file for
     each source file that we're going to wrap. These files are named
     :file:`quippy_{STEM}_wrap.f90` and are designed to use the
     restricted set of Fortran 90 features understood by f2py.

   - The :file:`quippy_{STEM}_wrap.f90` files are passed to f2py, which 
     generates a Python extension module :mod:`_quippy`. This is a low-level
     module which allows all the Fortran functions to be called from Python,
     but doesn't know anything about derived-types. See :ref:`wrapping-fortran-90-code`
     for more details.
   
   - All the Fortran sources - both those imported and the generated
     wrappers - are compiled using the Fortran compiler specified in
     the :makevar:`QUIPPY_COMPILER` Makefile variable. The :mod:`_quippy`
     C extension module is also compiled.

   - Finally all the object files are linked, together with external
     libraries such as NetCDF and LAPACK, to create
     :file:`_quippy.so`, the Python extension module. Along with the
     various pure Python modules which make up quippy, this is
     installed in the standard place for Python extension modules on
     your system. This will probably be something like
     :file:`/usr/local/lib/python-2.{x}/site-packages`, unless you
     overrode this by setting :makevar:`QUIPPY_INSTALL_OPTS`.

If the complication fails with an error message, please send the full
output to me at james.kermode@kcl.ac.uk and I'll do my best to work
out what's going wrong.

Once you've compiled quippy successfully, try to import it to verify that
everything is working okay. You should change directory away from the source
tree to avoid importing the uncompiled version::

   $ cd ~
   $ python
   >>> Python 2.6.1+ (release26-maint:71045, Jul 17 2009, 17:47:18) 
   [GCC 4.3.3] on linux2
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import quippy

If you get an :exc:`ImportError` with a message about unresolved dependancies then
something went wrong with the linking process. 


Running the test suite
----------------------

Once quippy is successfully installed, you should run the test suite to 
check everything is working correctly::

   $ cd ${QUIP_ROOT}/Tools/quippy/tests
   $ python run_all.py

If any of the tests fail please send me (james.kermode@kcl.ac.uk) the output.


Installing the ipython profile
------------------------------

If you use `ipython`_ and have installed `matplotlib`_, there's a
special quippy profile you can install. Copy the files
:file:`quippy_load.py` and :file:`ipythonrc-quippy` from
:file:`${QUIP_ROOT}/Tools/quippy` to your :file:`~/.ipython` directory.
Invoking ipython as ``ipython -p quippy`` sets up matplotlib and
imports all the quippy functionality when you start ipython. This is
equivalent to ``ipython -pylab`` followed by ``from quippy import *``.

I use a shell alias which maps ``ipythonq`` to ``ipython -p quippy``
to save typing.
