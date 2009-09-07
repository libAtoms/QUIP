Installation
============

The aim of quippy is to make all the functions, subroutines and types
defined in QUIP available from Python scripts. We do this by creating
a Python extension which wraps all of these routines automatically
using `f2py <http://www.scipy.org/F2py>`_, part of 
`numpy <http://www.scipy.org/numpy>`_. f2py is a very useful tool, but it does
not yet support Fortran 90 extensions such as derived types, so quippy
has to do a little more work to allow us to use them transparently.

Requirements
------------

Essential:
 * `Python 2.4 <http://www.python.org>`_ or later
 * `numpy`_  - tested with version 1.2.1
 * A fortran compiler:

   * ifort 10 or later (tested with 10.1-015)
   * gfortran, svn version of 4.3.3 branch, available by
     anonymous svn from svn://gcc.gnu.org/svn/gcc/tags/gcc_4_3_3_release.
     With 4.3.2 or earlier you run into this 
     `bug <http://gcc.gnu.org/bugzilla/show_bug.cgi?id=37735>`_, 
     reported by Steve Winfield and now fixed.
   * PGI f90 ?

Optional:
 * `ipython <http://ipython.scipy.org>`_ makes using python interactively much more productive.
 * AtomEye atomistic configuration viewer - available from QUIP `svn repository <http://src.tcm.phy.cam.ac.uk/viewvc/jrk33/repo/trunk/AtomEye>`_

Getting quippy
--------------

If you have an account on the QUIP svn repository server, you should
checkout a copy of the QUIP trunk::

  svn checkout svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/QUIP QUIP

If you don't have an account, you can download a `snapshot
<svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/QUIP?view=tar>`_
of the svn repository as a tarball. Extracting this archive will create a
:file:`QUIP` tree.

The structure of the QUIP source tree is::

  ${QUIP_ROOT}
     Makefiles
     libAtoms
     QUIP_Core
     QUIP_Utils
     QUIP_Programs
     Tests
     Tools
         quippy
         AtomEye

Here we've used the environment variable :envvar:`QUIP_ROOT` to refer
to the root of the QUIP tree. The :file:`libAtoms`, :file:`QUIP_Core`,
:file:`QUIP_Util` and :file:`QUIP_Programs` directories contain
Fortran 95 source code. quippy itself lives in the
:file:`{QUIP_ROOT}/Tools/quippy` directory

If you want to include the AtomEye extension, you should check it out
under :file:`${QUIP_ROOT}/Tools`::

  cd ${QUIP_ROOT}/Tools
  svn checkout svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/AtomEye AtomEye

Again, if you don't have an account you can download a `snapshot
<svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/QUIP?view=tar>`_, which you should
extract under :file:`${QUIP_ROOT}/Tools`

Patching numpy
--------------

There's a patch file which needs to be applied to numpy. This allows
optional arguments to be passed correctly between python and fortran.
To apply the patch, you need to find where numpy is installed::

   $ python
   >>> import numpy
   >>> print numpy.__file__

Then go to this directory, and change to the f2py subdirectory, and apply
the patch. For example::

  $ cd ~/lib/python2.6/site-packages/numpy/f2py
  $ patch < ~/Code/QUIP/Tools/quippy/numpy.f2py.patch


Configuring quippy
------------------

Before compiling quippy, you need to set the environment variable
:envvar:`QUIP_ARCH`, for example on 64-bit Mac OS X with gfortran you
would do::

  $ export QUIP_ARCH=darwin_x86_64_gfortan

This sets the architecture for the QUIP framework which
quippy wraps. Architectures are defined by creating a file
:file:`${QUIP_ROOT}/Makefiles/Makefile.${QUIP_ARCH}` which describes which
compilers and libraries should be used and where they can be found. quippy has
been tested on the following architectures::

  darwin_x86_64_gfortran
  linux_x86_64_gfortran
  linux_x86_64_ifort_gcc
  linux_x86_64_pgi

The :file:`Makefile.${QUIP_ARCH}` should define various quippy-specific variables:

:makevar:`QUIPPY_FCOMPILER`
   Fortran compiler to use. The shell command::

     $ f2py -c --help-fcompiler 

   will print a list of detected compilers on your system. Use :value:`gnu95` for gfortran, 
   :var:`intel` for ifort on 32-bit platforms and :var:`intelem` for ifort on 64-bit platforms.

:makevar:`QUIPPY_INSTALL_OPTS`
   Installation options, e.g. specify --home=${HOME} 
   or --prefix=${PREFIX} to install in a non-default location.

:makevar:`QUIPPY_F90FLAGS`
   Flags to pass to Fortran compiler

:makevar:`QUIPPY_OPT`
   Optimisation settings for Fortran compiler

:makevar:`QUIPPY_CPP`
   Fortran preprocessor to use 

:makevar:`QUIPPY_NO_QUIP`
   If set to 1, omit compilation of QUIP modules, including only low-level libAtoms routines

:makevar:`QUIPPY_NO_TOOLS`
   If set to 1, omit compilation of extra tools

:makevar:`QUIPPY_NO_CRACK`
  If set to 1, omit compilation of crack utilities. Currently necessary with ifort to 
  avoid an internal compiler error. 

Compilation
-----------

To compile quippy, run::

   $ export QUIP_ARCH=<<arch>>
   $ make
   $ make install

Important QUIP Makefile variables include

This should build everything and install the Python modules in
the default location for your python installation.

To check everything is working try the following::

   $ python 
   >>> import quippy

If you are using gfortran, the QUIP bindings require a recent svn
version of gfortran to compile (with 4.3.2 or earlier you run into
this bug, reported by Steve Winfield and now fixed:
http://gcc.gnu.org/bugzilla/show_bug.cgi?id=37735)

