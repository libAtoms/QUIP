# QUIP - QUantum mechanics and Interatomic Potentials

The `QUIP` package is a software library written in Fortran 95+ for the
purposes of carrying out molecular dynamics simulations. The code
implements a wide variety of interatomic potentials and tight binding
quantum mechanics, and is also able to call external packages. Various
hybrid combinations are also supported in the style of QM/MM.

For more details, see the [online documentation](http://libatoms.github.io/QUIP).

Portions of this code were written by: Albert Bartok-Partay, Livia
Bartok-Partay, Federico Bianchini, Anke Butenuth, Marco Caccin,
Silvia Cereda, Gabor Csanyi, Alessio Comisso, Tom Daff, Sebastian
John, Chiara Gattinoni, Gianpietro Moras, James Kermode, Letif
Mones, Alan Nichol, David Packwood, Lars Pastewka, Giovanni
Peralta, Ivan Solt, Oliver Strickson, Wojciech Szlachta, Csilla
Varnai, Steven Winfield.

Copyright 2006-2015.

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
 - GAP (Gaussian Approximation Potentials: general many-body)
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

The following tight-binding functional forms and parametrisations are implemented:

 - Bowler
 - DFTB
 - GSP
 - NRL-TB

The following external packages can be called:

 - CASTEP
 - VASP
 - CP2K
 - ASAP
 - ASE (latest svn trunk recommended)
 - Molpro

## Code Philosophy

We try to strike a compromise between readability of code and
efficiency, and think of QUIP/libAtoms as a "developer's code": nice
when you want to try new ideas quickly, but not competitive in
efficiency with other major md codes such as LAMMPS, Gromacs etc. We
use several extensions to the Fortran 95 standard in order to make the
coding style more object oriented. Several compilers support all the
necessary extensions in their recent versions, e.g. GNU v4.4 and
later. Support in the Intel compiler suite is there in principle, but
not every recent version has correct implementation, although we have
not encountered many problems past version 11.

## Compilation Instructions

1.  Decide your architecture by looking at the Makefiles/README, and
    define an environmental variable QUIP_ARCH, e.g.
    
		export QUIP_ARCH=linux_x86_64_ifort_icc
    
    You may well need to create your own
    `Makefiles/Makefile.${QUIP_ARCH}` file based on an existing file.
    
2.  Ensure that you have sufficiently up-to-date compilers. If you are
    using GNU compiler suite, you need version 4.4 or later. From
    Intel, you need version > 11.0.084.
    
3.  Customise QUIP, set the maths libraries and provide linking options:
    
        make config
    
    Makefile.config will create a build directory, `build.${QUIP_ARCH}`,
    and all the building happen there. First it will ask you some
    questions about where you keep libraries and other stuff, if you
    don't use something it is asking for, just leave it blank. The
    answers will be stored in `Makefile.inc` in the `build.${QUIP_ARCH}`
    directory, and you can edit them later (e.g. to change optimisation
    or debug options).  Note that the default state is usually with
    rather heavy debugging on, including bounds checking, which makes
    the code quite slow.
    
    If you later make significant changes to the configuration such as
    enabling or disabling tight-binding support you should force a
    full rebuild by doing a `make deepclean; make`.
    
4.  Compile all modules, libraries and programs with:
    
		make

	From the top-level `QUIP` directory. All object files and programs
    are built under `build/${QUIP_ARCH}/libquip.a`.

    Other useful make targets include:

	- `make install` : copies all compiled programs it can find to
		`QUIP_INSTALLDIR`, if it's defined and is a directory
	- `make libquip`:   Compile QUIP as a library and link to it. 
	  This will make all the various libraries and combine them into one:
	  `build/${QUIP_ARCH}/libquip.a`, which is what you need to link (and
	  of course LAPACK).
    
5.  QUIP is a "developer's" code, and is not for the faint
    hearted. A good start is to use the `quip` program, which is under
    `src/Programs`, and allows the evaluation of properties of an atomic
    configuration using a variety of models. For example:
    
    		quip at_file=test.xyz init_args='IP LJ' \
    			param_file=QUIP_Core/parameters/ip.parms.LJ.xml E
    
    assuming that you have a file called `test.xyz` with the following
    data in it representing Cu atoms in a simple cubic lattice:
    
    		1
    		Lattice="4 0 0 0 4 0 0 0 4" Properties=Z:I:1:pos:R:3
    		29 0 0 0 
    		
    The Lennard-Jones parameters in the above example is defined in the
    `ip.parms.LJ.xml` file. 
    
    Most string arguments can be replaced by `--help` and QUIP programs
    will  then print  a  list  of allowable  keywords  with brief  help
    messages as to their usage,  so e.g. `init_args=--help` will give a
    list of potential model  types (and some combinations). The parsing
    is recursive,  so `init_args=IP --help`  will then proceed  to list
    the types of interatomic potentials (IP) that are available.

6.  To compile the Python wrappers (`quippy`), run

		make quippy

	`quippy` depends on Python and [numpy](http://www.numpy.org), the
	[Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
	(ASE) and some features require [scipy](http://www.scipy.org)
	and/or [matscipy](https://github.com/libAtoms/matscipy).
	
	More details on the quippy installation process are available in
	the [online documentation](http://libatoms.github.io/QUIP/).

7.  To run the unit and regression tests, which depend on `quippy`:

		make test
    
8.  Some functionality is only available if you check out other
	modules within the `QUIP/src/` directories, e.g. the `ThirdParty`
	(DFTB parameters, TTM3f water model), `GAP` (Gaussian
	Approximation Potential models) and `GAP-filler` (Gaussian
	Approximation Potential model training).
