# QUIP - QUantum mechanics and Interatomic Potentials

[![Build Status](https://travis-ci.org/libAtoms/QUIP.svg?branch=public)](https://travis-ci.org/libAtoms/QUIP)

The `QUIP` package is a collection of software tools to carry out
molecular dynamics simulations. It implements a variety of interatomic
potentials and tight binding quantum mechanics, and is also able to
call external packages, and serve as plugins to other software such as
[LAMMPS](http://lammps.sandia.gov), [CP2K](http://www.cp2k.org)
and also the python framework [ASE](https://wiki.fysik.dtu.dk/ase).
Various hybrid combinations are also supported in the style of QM/MM,
with a particular focus on materials systems such as metals and
semiconductors.

For more details, see the [online documentation](http://libatoms.github.io/QUIP).

Long term support of the package is ensured by:
 - Noam Bernstein (Naval Research Laboratory)
 - Gabor Csanyi (University of Cambridge)
 - James Kermode (University of Warwick)

Portions of this code were written by: Albert Bartok-Partay, Livia
Bartok-Partay, Federico Bianchini, Anke Butenuth, Marco Caccin,
Silvia Cereda, Gabor Csanyi, Alessio Comisso, Tom Daff, ST John,
Chiara Gattinoni, Gianpietro Moras, James Kermode, Letif Mones,
Alan Nichol, David Packwood, Lars Pastewka, Giovanni Peralta, Ivan
Solt, Oliver Strickson, Wojciech Szlachta, Csilla Varnai, Steven
Winfield.

Copyright 2006-2016.

Most of the publicly available version is released under the GNU
General Public license, version 2, with some portions in the public
domain.

## Features

The following interatomic potentials are presently coded or linked in QUIP:

 - BKS (van Beest, Kremer and van Santen) (silica)
 - EAM (fcc metals)
 - Fanourgakis-Xantheas (water)
 - Finnis-Sinclair (bcc metals)
 - Flikkema-Bromley
 - [GAP](http://www.libatoms.org/gap/gap_download.html) (Gaussian Approximation Potentials)
 - Guggenheim-!McGlashan
 - Brenner (carbon)
 - OpenKIM (general interface)
 - Lennard-Jones
 - Morse
 - Partridge-Schwenke (water monomer)
 - Stillinger-Weber (carbon, silicon, germanium)
 - SiMEAM (silicon)
 - Sutton-Chen
 - Tangney-Scandolo (silica, titania etc)
 - Tersoff (silicon, carbon)

The following tight-binding functional forms and parametrisations are
implemented:

 - Bowler
 - DFTB
 - GSP
 - NRL-TB

The following external packages can be called:

 - CASTEP
 - VASP
 - CP2K
 - ASAP
 - Molpro
 - ASE (required if using ``quippy`` Python interface; latest version 
   recommended)

## Code philosophy and goals

QUIP was born because of the need to efficiently tie together a wide
variety of different models, both empirical and quantum mechanical. It
will not be competitive in terms of performance with codes such as LAMMPS
and Gromacs, but has a number of unique features:

- Deep access to most of the Fortran types and routines from Python via the
  ``quippy`` package
- Support for Gaussian Approximation Potentials [GAP](http://www.libatoms.org/gap/gap_download.html)
- Does not assume minimum image convention, so interatomic potentials can
  have cutoffs that are larger than the unit cell size

## Precompiled Containers

If you have access to [Docker](http://singularity.lbl.gov) or
[Singularity](http://singularity.lbl.gov), you can try one of the
[precompiled images](https://github.com/libAtoms/QUIP/blob/public/docker/README.md)
to get up and running quickly.

## Compilation Instructions

1.  To compile QUIP the minimum requirements are:

    - A working Fortran compiler. QUIP is tested with ``gfortran`` 4.4 and
      later, and ``ifort`` 11.1.

    - Linear algebra libraries BLAS and LAPACK. QUIP is tested with
      reference versions ``libblas-dev`` and ``liblapack-dev`` on Ubuntu
      12.04, and ``mkl`` 11.1 with ``ifort``.

2.  Clone the QUIP repository from GitHub. The ``--recursive`` option
    brings in submodules automatically (If you don't do this, then
    you will need to run ``git submodule update --init``
    from the top-level QUIP directory after cloning) ::

        git clone --recursive https://github.com/libAtoms/QUIP.git

3.  Decide your architecture by looking in the ``arch/`` directory, and
    define an environmental variable ``QUIP_ARCH``, e.g.::

        export QUIP_ARCH=linux_x86_64_gfortran

    for standard gfortran on Linux. You may need to create your own
    ``arch/Makefile.${QUIP_ARCH}`` file based on an existing file for
    more exotic systems.

4.  Customise QUIP, set the maths libraries and provide linking options::

        make config

    Makefile.config will create a build directory, ``build/${QUIP_ARCH}``,
    and all the building happen there. First it will ask you some
    questions about where you keep libraries and other stuff, if you
    don't use something it is asking for, just leave it blank. The
    answers will be stored in ``Makefile.inc`` in the ``build/${QUIP_ARCH}``
    directory, and you can edit them later (e.g. to change optimisation
    or debug options).

    If you later make significant changes to the configuration such as
    enabling or disabling tight-binding support you should force a
    full rebuild by doing a `make deepclean; make`.

5.  Compile all programs, modules and libraries::

        make

    From the top-level ``QUIP`` directory. All programs are built in
    ``build/${QUIP_ARCH}/``. You can also find compiled object files
    and libraries (``libquip.a``) in that directory. Programs can be
    called directly from that directory.

    Other useful make targets include:

    - ``make install`` : copies all compiled programs it can find to
      ``QUIP_INSTALLDIR``, if it's defined and is a directory (full path
      required), and copies bundled structures to ``QUIP_STRUCTS_DIR``
      if it is defined.

    - ``make libquip``:   Compile QUIP as a library and link to it.
      This will make all the various libraries and combine them into one:
      ``build/${QUIP_ARCH}/libquip.a``, which is what you need to link with
      (as well as LAPACK).

6.  A good starting point is to use the ``quip`` program, which can
    calculate the properties of an atomic configuration using a
    variety of models. For example::

        quip at_file=test.xyz init_args='IP LJ' \
            param_file=share/Parameters/ip.parms.LJ.xml E

    assuming that you have a file called ``test.xyz`` with the following
    data in it representing Cu atoms in a cubic fcc lattice::

        4
        Lattice="3.61 0 0 0 3.61 0 0 0 3.61" Properties=species:S:1:pos:R:3
        Cu     0.000 0.000 0.000
        Cu     0.000 1.805 1.805
        Cu     1.805 0.000 1.805
        Cu     1.805 1.805 0.000

    The Lennard-Jones parameters in the above example are defined in the
    ``ip.parms.LJ.xml`` file under ``share/Parameters`` (ensure the path
    to this file is correct). The format of the atomic configuration is 
    given in
    [Extended XYZ](http://libatoms.github.io/QUIP/io.html#extendedxyz)
    format, in which the first line is the number of atoms, the second line
    is a series of key=value pairs, which must at least contain the Lattice
    key giving the periodic bounding box and the Properties key that
    describes the remaining lines. The value of Properties is a sequence of
    triplets separated by a colon (:), that give the name, type and number
    of columns, with the type given by I for integers, R for reals, S for
    strings.

    Most string arguments can be replaced by ``--help`` and QUIP programs
    will  then print  a  list  of allowable  keywords  with brief  help
    messages as to their usage,  so e.g. ``init_args=--help`` will give a
    list of potential model  types (and some combinations). The parsing
    is recursive,  so ``init_args="IP --help"``  will then proceed  to list
    the types of interatomic potentials (IP) that are available.

7.  To compile the Python wrappers (``quippy``), the minimum requirements
    are:

    - Python 2.7. (Python 3 support is in progress).

    - [NumPy](http://www.numpy.org)

    - [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) (``ase``)

    - (optional) [SciPy](http://www.scipy.org)

    - (optional) [matscipy](https://github.com/libAtoms/matscipy).

8.  If you are using a Python virtual environment (virtualenv) and would like
    to install ``quippy`` into it, ensure the environment is activated
    (``source <env_dir>/bin/activate``, where ``<env_dir>`` is the root of
    your virtual environment) _before_ building ``quippy`` (otherwise library
    versions may cause unexpected conflicts).

9.  To compile the Python wrappers (``quippy``), run::

        make quippy

    Quippy can be used by adding the ``lib`` directory in
    ``quippy/build/${QUIP_ARCH}`` to your ``$PYTHONPATH``, however it can be
    more convenient to install into a specific Python distribution::

        make install-quippy

    will either install into the current virtualenv or attempt to install
    systemwide (usually fails without ``sudo``). To install only for the
    current user (into ``~/.local``), execute the command
    ``QUIPPY_INSTALL_OPTS=--user make install-quippy``,
    or use ``QUIPPY_INSTALL_OPTS=--prefix=<directory>`` to install into a
    specific directory. ``QUIPPY_INSTALL_OPTS`` can also be set in the file
    ``build/${QUIP_ARCH}/Makefile.inc``.

10.  More details on the quippy installation process and troubleshooting for
    common build problems are available in the 
    [online documentation](http://libatoms.github.io/QUIP/).

11.  To run the unit and regression tests, which depend on ``quippy``::

        make test

12.  To get back to a state near to a fresh clone, use

        make distclean

13. Some functionality is only available if you check out other
    modules within the ``QUIP/src/`` directories, e.g. the ``ThirdParty``
    (DFTB parameters, TTM3f water model), ``GAP`` (Gaussian
    Approximation Potential models) and ``GAP-filler`` (Gaussian
    Approximation Potential model training). These packages are
    not distributed with QUIP because they come with different licensing
    restrictions, but you can get them [here](http://www.libatoms.org/gap/gap_download.html)

    GAP is a machine learning method that uses Gaussian process
    regression, and needs large data files to run. You can find
    potentials that have been published as well as training data in
    our [data repository](http://www.libatoms.org/Home/DataRepository).

14. In order to run QUIP potentials via LAMMPS, ``make libquip`` to get QUIP
    into library form, and then follow the instructions in the
    [LAMMPS documentation](http://lammps.sandia.gov/doc/pair_quip.html).

