# QUIP - Quantum Mechanics and Interatomic Potentials

This is the top level directory for QUIP, libAtoms and friends. This
package is a software library written in Fortran 95+ for the purposes
of carrying out molecular dynamics simulations. The QUIP package,
built on top of libAtoms, implements a wide variety of interatomic
potentials and tight binding quantum mechanics, and is also able to
call external packages. Various hybrid combinations are also supported
in the style of QM/MM.

For more details, see the [online documentation](|http://libatoms.github.io/QUIP).

Most of the publicly available version is released under the GNU
General Public license, version 2, with some portions in the public
domain.

The main libAtoms/QUIP contributors are:

 - University of Cambridge: Albert P. Bartók, Gábor Csányi, Alan
   Nichol, Letif Mones, Wojciech Szlachta
 - University of Warwick: James Kermode
 - King's College London: Alessandro De Vita
 - Naval Research Laboratory, Washington DC: Noam Bernstein
 - Fraunhofer IWM, Freiburg: Lars Pastewka

## Potentials

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

1) decide your architecture by looking at the Makefiles/README, and
   define an environmental variable QUIP_ARCH, e.g. 

   export QUIP_ARCH=linux_x86_64_ifort_icc

   You may well need to create your own
   Makefiles/Makefile.${QUIP_ARCH} file based on an existing file.

2) Ensure that you have sufficiently up-to-date compilers. If you are
   using GNU compiler suite, you need version 4.4 or later. From
   Intel, you need version > 11.0.084.

3) customise QUIP, set the maths libraries and provide linking options:

   make config

   Makefile.config will create a build directory, build.${QUIP_ARCH},
   and all the building happen there. First it will ask you some
   questions about where you keep libraries and other stuff, if you
   don't use something it is asking for, just leave it blank. The
   answers will be stored in Makefile.inc in the build.${QUIP_ARCH}
   directory, and you can edit them later (e.g. to change optimisation
   or debug options).  Note that the default state is usually with
   rather heavy debugging on, including bounds checking, which makes
   the code quite slow.

   If you later make significant changes to the configuration such as
   enabling or disabling tight-binding support you should do forcee a
   full rebuild by doing a 'make deepclean; make'.

4) make something:

   make libAtoms
 
   or (usually)

   make QUIP_Programs/<progname>

   Note that the `make' command has to be executed from the top level
   directory, even for targets in subdirectories.

   Most useful make targets include
    all : pretty much every vaguely useful top level program
    QUIP_Programs/eval : evaluate energies, forces, minimize energy etc
    QUIP_Programs/md : basic md program
    QUIP_Core : generates everything that is necessary to use QUIP as a library
    install : copies all compiled programs it can find to QUIP_INSTDIR,
      if it's defined and is a directory

   You can also use QUIP and libAtoms as a library and link to it. To
   make the library version only, execute

   make libquip

   This will make all the various libraries and combine them into one:
   build.${QUIP_ARCH}/libquip.a, which is what you need to link (and
   of course LAPACK).
 
5) QUIP/libatoms is a "developer's" code, and is not for the faint
   hearted. A good start is to use the 'eval' program, which is under
   QUIP_Programs, and allows the evaluation of properties of an atomic
   configuration using a variety of models. For example:

   eval at_file=test.xyz init_args='IP LJ' \
      param_file=QUIP_Core/parameters/ip.parms.LJ.xml E

   assuming that you have a file called test.xyz with the following
   data in it (without the dashes, obviously) representing Cu atoms in
   a simple cubic lattice:

   ----
   1
   Lattice="4 0 0 0 4 0 0 0 4" Properties=Z:I:1:pos:R:3
   29 0 0 0 
   ----

   The Lennard-Jones parameters in the above example is defined in the
   ip.parms.LJ.xml file. 

   Most string arguments can be replaced by '--help' and QUIP programs
   will  then print  a  list  of allowable  keywords  with brief  help
   messages as to their usage,  so e.g. 'init_args=--help' will give a
   list of potential model  types (and some combinations). The parsing
   is recursive,  so 'init_args=IP --help'  will then proceed  to list
   the types of interatomic potentials (IP) that are available.

Some functionality is only available if you check out other packages
within the QUIP tree, e.g. the ThirdParty (DFTB parameters, TTM3f
water model), GAP (Gaussian Approximation Potential models) and
GAP-filler (Gaussian Approximation Potential model training).

