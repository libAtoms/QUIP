# QUIP - QUantum mechanics and Interatomic Potentials
[![Build Status](https://img.shields.io/badge/docs-public-brightgreen)](https://libatoms.github.io/QUIP/)
![build](https://github.com/libAtoms/QUIP/actions/workflows/Build.yml/badge.svg)
[![Docker Pulls](https://img.shields.io/docker/pulls/libatomsquip/quip)](https://hub.docker.com/r/libatomsquip/quip)

The `QUIP` package is a collection of software tools to carry out
molecular dynamics simulations. It implements a variety of interatomic
potentials and tight binding quantum mechanics, and is also able to
call external packages, and serve as plugins to other software such as
[LAMMPS](http://lammps.sandia.gov), [CP2K](http://www.cp2k.org)
and also the python framework [ASE](https://wiki.fysik.dtu.dk/ase).
Various hybrid combinations are also supported in the style of QM/MM,
with a particular focus on materials systems such as metals and
semiconductors.

For more details, see the [online documentation](http://libatoms.github.io/QUIP). There is separate [documentation for SOAP and GAP](https://libatoms.github.io/GAP).

Long term support of the package is ensured by:
 - Noam Bernstein ([@bernstei](https://github.com/bernstei), Naval Research Laboratory)
 - Gabor Csanyi ([@gabor1](https://github.com/gabor1), University of Cambridge)
 - James Kermode ([@jameskermode](https://github.com/jameskermode), University of Warwick)

Portions of this code were written by: Albert Bartok-Partay, Livia
Bartok-Partay, Federico Bianchini, Anke Butenuth, Marco Caccin,
Silvia Cereda, Gabor Csanyi, Alessio Comisso, Tom Daff, ST John,
Chiara Gattinoni, Gianpietro Moras, James Kermode, Letif Mones,
Alan Nichol, David Packwood, Lars Pastewka, Giovanni Peralta, Ivan
Solt, Oliver Strickson, Wojciech Szlachta, Csilla Varnai, Steven
Winfield, Tamas K Stenczel, Adam Fekete.

Copyright 2006-2021.

Most of the publicly available version is released under the GNU
General Public license, version 2, with some portions in the public
domain. The GAP code, included as a submodule, is distributed under
a non-commercial [academic source license](https://github.com/libAtoms/GAP/blob/main/LICENSE.md)

## Citing QUIP, quippy and GAP

Please cite the following publication if you use QUIP:

```bibtex
@ARTICLE{Csanyi2007-py,
  title   = "Expressive Programming for Computational Physics in Fortran 95+",
  author  = "Cs{\'a}nyi, G{\'a}bor and Winfield, Steven and Kermode, J R and De
             Vita, A and Comisso, Alessio and Bernstein, Noam and Payne,
             Michael C",
  journal = "IoP Comput. Phys. Newsletter",
  pages   = "Spring 2007",
  year    =  2007
}
```

If you use the `quippy` Python interface, please cite:

```bibtex
@ARTICLE{Kermode2020-wu,
  title    = "f90wrap: an automated tool for constructing deep Python
              interfaces to modern Fortran codes",
  author   = "Kermode, James R",
  journal  = "J. Phys. Condens. Matter",
  month    =  mar,
  year     =  2020,
  keywords = "Fortran; Interfacing; Interoperability; Python; Wrapping codes;
              f2py",
  language = "en",
  issn     = "0953-8984, 1361-648X",
  pmid     = "32209737",
  doi      = "10.1088/1361-648X/ab82d2"
}
```

If you use the GAP code please cite

```bibtex

@ARTICLE{Bartok2010-pw,
  title    = "Gaussian approximation potentials: the accuracy of quantum
              mechanics, without the electrons",
  author   = "Bart{\'o}k, Albert P and Payne, Mike C and Kondor, Risi and
              Cs{\'a}nyi, G{\'a}bor",
  journal  = "Phys. Rev. Lett.",
  volume   =  104,
  number   =  13,
  pages    = "136403",
  month    =  apr,
  year     =  2010,
  issn     = "0031-9007, 1079-7114",
  pmid     = "20481899",
  doi      = "10.1103/PhysRevLett.104.136403"
}
```

## Features

The following interatomic potentials are presently coded or linked in QUIP:

 - BKS (van Beest, Kremer and van Santen) (silica)
 - EAM (fcc metals)
 - Fanourgakis-Xantheas (water)
 - Finnis-Sinclair (bcc metals)
 - Flikkema-Bromley
 - GAP (Gaussian Approximation Potentials) with (growing...) [online documentation](https://libatoms.github.io/GAP)
 - Guggenheim-McGlashan
 - Brenner (carbon)
 - OpenKIM (general interface)
 - Lennard-Jones
 - MBD (many-body dispersion correction)
 - Morse
 - Partridge-Schwenke (water monomer)
 - Stillinger-Weber (carbon, silicon, germanium)
 - SiMEAM (silicon)
 - Sutton-Chen
 - Tangney-Scandolo (silica, titania etc)
 - Tersoff (silicon, carbon)
 - Tkatchenko-Sheffler pairwise dispersion correction

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
 - ASE (required if using `quippy` Python interface; latest version
   recommended)

## Code philosophy and goals

QUIP was born because of the need to efficiently tie together a wide
variety of different models, both empirical and quantum mechanical. It
will not be competitive in terms of performance with codes such as LAMMPS
and Gromacs. The Atomic Simulation Environment also does this, and is
much more widely used, but QUIP has a number of unique features:

- Access to Fortran types and routines from Python via the
  `quippy` package
- Support for Gaussian Approximation Potentials (GAP) - [online docs](https://libatoms.github.io/GAP)
- Does not assume minimum image convention, so interatomic potentials can
  have cutoffs that are larger than the periodic unit cell size

## Binary Installation of QUIP and quippy

Binary 
for QUIP and the associated quippy Python bindings
that provide interopability with the Atomic Simulation Environment (ASE) are
available from the [Python package index](https://pypi.org/project/quippy-ase/)
(PyPI) under the package name `quippy-ase`.
This means you can install the latest release with:

```bash
pip install quippy-ase
```

Installing via `pip` also makes the `quip` and `gap_fit` command line
programs available (providing the [directory that pip installs scripts
to](https://stackoverflow.com/questions/62162970/programmatically-determine-pip-user-install-location-scripts-directory/62167797#62167797) is on your `PATH`).

Currently, wheels are available for `x86_64` architectures
with Python 3.9+ on macOS and glibc-based Linux distributions
(e.g. Ubuntu, CentOS) and for macOS arm64 (Apple Silicon). The wheels are updated periodically
using GitHub Actions CI. Please open [issues](https://github.com/libAtoms/QUIP/issues)
here if you have problems installing with `pip`.

## Precompiled Containers

If you have access to [Docker](https://hub.docker.com) or
[Singularity](http://singularity.lbl.gov), you can try one of the
[precompiled images](https://github.com/libAtoms/quip-docker)
to get up and running quickly.

## Compilation Instructions

1.  To compile QUIP the minimum requirements are:

    - A working Fortran compiler. QUIP is tested with `gfortran` 4.4 and
      later, and `ifort` 11.1.

    - Linear algebra libraries BLAS and LAPACK. QUIP works with
      reference versions, OpenBLAS (recommended), or Intel MKL.

    - [Meson build system](https://mesonbuild.com/) version 1.1 or later.
      Install with `pip install meson` or via your package manager.

    - [Ninja build tool](https://ninja-build.org/).
      Install with `pip install ninja` or via your package manager.

    - MPI (optional): To use MPI parallelization of `gap_fit`, you need a
      ScaLAPACK library, e.g. `libscalapack-openmpi` on Ubuntu, or
      as part of MKL.

2.  Clone the QUIP repository from GitHub. The `--recursive` option
    brings in submodules automatically (If you don't do this, then
    you will need to run `git submodule update --init --recursive`
    from the top-level QUIP directory after cloning):
    ```bash
    git clone --recursive https://github.com/libAtoms/QUIP.git
    ```

    One submodule is the GAP code, which can be found in `src/GAP`.
    Note that GAP is distributed under a different
    [license](https://github.com/libAtoms/GAP/blob/main/LICENSE.md).

    GAP is a machine learning method that uses Gaussian process
    regression, and needs large data files to run. You can find
    potentials that have been published as well as training data in
    our [data repository](https://libatoms.github.io/GAP/data.html), see also the [online docs](https://libatoms.github.io/GAP).

3.  **Important:** Ensure the submodules are on the correct branches with meson support:
    ```bash
    cd src/fox && git checkout master && git pull
    cd ../GAP && git checkout main && git pull
    cd ../..
    ```

4.  Configure the build with meson:
    ```bash
    meson setup builddir
    ```

    This creates a `builddir` directory where all build artifacts will be placed.
    Meson will automatically detect your compiler, BLAS/LAPACK libraries, and configure
    the build appropriately.

    **Build options**: You can customize the build using `-D` flags:
    ```bash
    meson setup builddir -Dgap=true -Dmpi=false
    ```

    Available options:
    - `gap` (default: `true`): Enable GAP (Gaussian Approximation Potentials) support
    - `mpi` (default: `false`): Enable MPI parallelization

    To reconfigure an existing build directory:
    ```bash
    meson configure builddir -Dmpi=true
    ```

5.  Compile all programs, modules and libraries:
    ```bash
    meson compile -C builddir
    ```

    All programs are built in `builddir/src/Programs/`. You can also find compiled
    shared libraries in `builddir/src/*/lib*.so`. Programs can be called directly from
    the build directory:
    ```bash
    ./builddir/src/Programs/quip --help
    ./builddir/src/Programs/gap_fit --help
    ```

    Other useful meson commands:

    - `meson install -C builddir` : Install compiled programs and libraries to the
      system (or use `--destdir` to specify an installation prefix)

    - `meson test -C builddir` : Run the test suite (requires quippy to be built)

    - Clean and rebuild:
      ```bash
      rm -rf builddir
      meson setup builddir
      meson compile -C builddir
      ```

6.  A good starting point is to use the `quip` program, which can
    calculate the properties of an atomic configuration using a
    variety of models. For example:
    ```bash
    ./builddir/src/Programs/quip atoms_filename=test.xyz init_args='IP LJ' \
        param_filename=share/Parameters/ip.parms.LJ.xml E
    ```
    assuming that you have a file called `test.xyz` with the following
    data in it representing Cu atoms in a cubic fcc lattice::
    ```
    4
    Lattice="3.61 0 0 0 3.61 0 0 0 3.61" Properties=species:S:1:pos:R:3
    Cu     0.000 0.000 0.000
    Cu     0.000 1.805 1.805
    Cu     1.805 0.000 1.805
    Cu     1.805 1.805 0.000
    ```
    The Lennard-Jones parameters in the above example are defined in the
    `ip.parms.LJ.xml` file under `share/Parameters` (ensure the path
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

    Most string arguments can be replaced by `--help` and QUIP programs
    will  then print  a  list  of allowable  keywords  with brief  help
    messages as to their usage,  so e.g. `init_args=--help` will give a
    list of potential model  types (and some combinations). The parsing
    is recursive,  so `init_args="IP --help"`  will then proceed  to list
    the types of interatomic potentials (IP) that are available.

7.  To compile the Python wrappers (`quippy`), the minimum requirements
    are as follows:
    - Python 3.9+
    - [Meson build system](https://mesonbuild.com/) version 1.1+
    - [NumPy](http://www.numpy.org) (`numpy>=1.5.0`)
    - [Atomic Simulation Environment ](https://wiki.fysik.dtu.dk/ase/) (`ase>=3.17.0`)
    - [f90wrap](https://github.com/jameskermode/f90wrap) development version - install with:
      ```bash
      pip install git+https://github.com/jameskermode/f90wrap.git@master
      ```
    - (optional) [SciPy](http://www.scipy.org)
    - (optional) [matscipy](https://github.com/libAtoms/matscipy).

    **Recommended**: Use [uv](https://docs.astral.sh/uv/) for fast, isolated environment management:
    ```bash
    # Install uv if not already installed
    pip install uv

    # Create isolated environment
    uv venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate

    # Install f90wrap development version
    uv pip install git+https://github.com/jameskermode/f90wrap.git@master

    # Install other dependencies
    uv pip install numpy ase meson ninja
    ```

    Note: If you are using a Python virtual environment (virtualenv or venv) and would like
    to install `quippy` into it, ensure the environment is activated
    (`source <env_dir>/bin/activate`, where `<env_dir>` is the root of
    your virtual environment) _before_ building `quippy`.

8.  To compile the Python wrappers (`quippy`):

    **Important:** You must first build the main QUIP libraries (steps 4-5 above) before building quippy.

    **Method 1: Using pip (recommended)**

    The simplest way to build and install quippy is with pip:
    ```bash
    cd quippy
    pip install .
    ```

    **Note:** Do NOT use `pip install -e .` (editable install) as this is problematic with meson-based builds.

    **Method 2: Using meson directly**

    Alternatively, you can build with meson and then install:
    ```bash
    cd quippy
    meson setup builddir
    meson compile -C builddir
    ```

    This will:
    - Preprocess the Fortran source files
    - Generate Python wrappers using f90wrap
    - Compile the Python extension module
    - Prepare the quippy package for installation

    Then install with pip:
    ```bash
    pip install .
    ```

9.  More details on the quippy installation process and troubleshooting for
    common build problems are available in the
    [online documentation](http://libatoms.github.io/QUIP/).

10.  To run the unit and regression tests (requires quippy to be installed):
    ```bash
    cd tests
    python3 run_all.py -v
    ```
    Or if using meson:
    ```bash
    meson test -C builddir
    ```

11. Some functionality is only available if you check out other
    modules within the `QUIP/src/` directories, e.g. the `ThirdParty`
    (DFTB parameters, TTM3f water model).

12. In order to run QUIP potentials via LAMMPS, build the QUIP libraries with meson,
    and then follow the instructions in the
    [LAMMPS documentation](http://lammps.sandia.gov/doc/pair_quip.html). You need at least 11 Aug 2017 version or later.
    The required libraries will be in `builddir/src/*/lib*.so`.

## Migrating from Make to Meson

If you have previously built QUIP using the old Makefile-based system, note the following differences:

**Key Changes:**
- **No `QUIP_ARCH` environment variable**: Meson automatically detects your system configuration
- **Build directory**: All build artifacts go to `builddir/` (or whatever name you choose) instead of `build/${QUIP_ARCH}/`
- **Configuration**: Use `meson configure` with `-D` options instead of editing `Makefile.inc`
- **Submodules**: Ensure submodules are updated to versions with meson support (fox: master, GAP: main)
- **Programs location**: Executables are in `builddir/src/Programs/` instead of `build/${QUIP_ARCH}/`
- **Libraries**: Shared libraries (`.so`) are built by default instead of static libraries (`.a`)

**Workflow Comparison:**

| Task | Old (Make) | New (Meson) |
|------|-----------|-------------|
| Setup | `export QUIP_ARCH=...` + `make config` | `meson setup builddir` |
| Configure | Edit `build/${QUIP_ARCH}/Makefile.inc` | `meson configure builddir -Doption=value` |
| Build | `make` | `meson compile -C builddir` |
| Install | `make install` | `meson install -C builddir` |
| Clean | `make clean` | `rm -rf builddir` |
| Test | `make test` | `meson test -C builddir` |

**To clean your old Make builds:**
```bash
# If you still have the old Makefile setup
make distclean

# Remove old build directories
rm -rf build/
```

# Developer notes:

## Fixing/updating the version of GAP:

  ```bash
  cd src/GAP
  ```
  ```bash
  git checkout <commit>
  ```
  OR
  ```bash
  git checkout main

  ```
  Updating the version in the `QUIP` repository:
  ```
  cd ../..
  git add src/GAP
  git commit -m "updating the version of GAP"
  ```

## Mac OS

We do not recommend Apple-shipped compilers and python, and we do not test
compatibility with them. Either use MacPorts or Homebrew to obtain GNU compilers,
and also use the python from there or Anaconda. As of this edit, gcc-8.1 produces as
internal compiler error, but gcc-4.6 through to gcc-7 is fine.

## Triggering the wheel build

Wheels are built on push and pull requests to `public` using [cibuildwheel](https://cibuildwheel.readthedocs.io/en/stable/)
with [this workflow](.github/workflows/build-wheels.yml).

To make a release candidate create a tag with a suffix such as `-rc1` for the first attempt,
push to trigger the build:

```bash
git commit -m 'release v0.x.z-rc1'
git tag v0.x.y-rc1
git push --tags
```

If all goes well, the `.whl` files will show up as assets within a new GitHub
release. The installation process can now be tested locally.

## Release wheels to PyPI

Once everything works correctly, make a full release (i.e. create a tag named
just `v0.x.y` without the `-rc1` suffix). This will trigger the upload of wheels
and source distribution to PyPI.
