#!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#!X
#!X     libAtoms: atomistic simulation library
#!X     
#!X     Copyright 2006-2007.
#!X
#!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
#!X     Contributors: Noam Bernstein, Alessio Comisso
#!X
#!X     The source code is released under the GNU General Public License,
#!X     version 2, http://www.gnu.org/copyleft/gpl.html
#!X
#!X     If you would like to license the source code under different terms,
#!X     please contact Gabor Csanyi, gabor@csanyi.net
#!X
#!X     When using this software, please cite the following reference:
#!X
#!X     http://www.libatoms.org
#!X
#!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

all: libatoms.a


include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif


LIBATOMS_F95_FILES = \
  mpi \
  System \
  ExtendableStr\
  Units \
  linearalgebra \
  Dictionary \
  Table \
  PeriodicTable \
  minimization \
  Atoms \
  Quaternions \
  RigidBody \
  Group \
  Constraints \
  Thermostat \
  DynamicalSystem \
  ParamReader \
  Spline \
  Sparse \
  CInOutput \
  clusters \
  Structures \
  frametools\
  nye_tensor \
  Topology \
  Atoms_ll \
  libAtoms \
  libAtoms_misc_utils

LIBATOMS_F77_FILES= \
   lbfgs	

LIBATOMS_C_FILES = \
  cutil\
  xyz_netcdf
#  malloc

LIBATOMS_F77_SOURCES = ${addsuffix .f, ${LIBATOMS_F77_FILES}}
LIBATOMS_F77_OBJS = ${addsuffix .o, ${LIBATOMS_F77_FILES}}
LIBATOMS_F95_SOURCES = ${addsuffix .f95, ${LIBATOMS_F95_FILES}}
LIBATOMS_F95_OBJS = ${addsuffix .o, ${LIBATOMS_F95_FILES}}
LIBATOMS_C_SOURCES = ${addsuffix .c, ${LIBATOMS_C_FILES}}
LIBATOMS_C_OBJS = ${addsuffix .o, ${LIBATOMS_C_FILES}}

LIBATOMS_OBJS = ${LIBATOMS_F77_OBJS} ${LIBATOMS_F95_OBJS} ${LIBATOMS_C_OBJS}


libatoms.a: ${LIBATOMS_OBJS}
	ar ${AR_ADD} libatoms.a $?

clean:
	rm -f libatoms.a *.o *.mod

allclean: clean docclean


${LIBATOMS_C_OBJS}  : %.o : %.c
	${CC} ${CFLAGS} -c $< 

md:  libatoms.a md.o
	$(LINKER)  -o $@ ${LINKFLAGS} $@.o  libatoms.a ${LINKOPTS}

xyz2pdb:  libatoms.a xyz2pdb.o
	$(LINKER)  -o $@ ${LINKFLAGS} $@.o  libatoms.a ${LINKOPTS}

xyz2xyz: xyz_netcdf.c
	${CC} ${CFLAGS} $< -o $@ ${LINKOPTS} -lm -DMAIN_PROGRAM

install-converters: xyz2xyz
	@echo ; echo "Where would you like the xyz and NetCDF converters installed to?"
	@read NCBINDIR && cp xyz2xyz $$NCBINDIR ;\
	ln -sf $$NCBINDIR/xyz2xyz $$NCBINDIR/nc2nc ;\
	ln -sf $$NCBINDIR/xyz2xyz $$NCBINDIR/ncstat ;\
	ln -sf $$NCBINDIR/xyz2xyz $$NCBINDIR/xyz2nc ;\
	ln -sf $$NCBINDIR/xyz2xyz $$NCBINDIR/nc2xyz ;\
	ln -sf $$NCBINDIR/xyz2xyz $$NCBINDIR/xyzstat

-include depend

