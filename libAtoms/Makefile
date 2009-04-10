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

include Makefile.common
ifeq (${MAKEFILE_ARCH_SUFFIX},)
  include Makefile.arch
else
  include Makefile.${MAKEFILE_ARCH_SUFFIX}
endif

ifneq (${BUILDDIR},)
	ifeq (${BUILDPRE},)
		export BUILDPRE=..
	endif
	ifeq (${VPATH},)
		export VPATH=..
	endif
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
  clusters \
  Structures \
  frametools\
  nye_tensor \
  Topology \
  CInOutput

ifeq (${HAVE_NETCDF},1)
  COMMON_F95FLAGS += -DHAVE_NETCDF ${ARCH_NETCDF_INCLUDEDIR} ${ARCH_NETCDF4}
  CFLAGS += -DHAVE_NETCDF ${ARCH_NETCDF_INCLUDEDIR}  ${ARCH_NETCDF4}
  ARCH_LIBS += ${ARCH_NETCDF_LIBDIR} ${ARCH_NETCDF_LIBS}
  LIBATOMS_F95_FILES +=  netcdf \
	  libAtoms \
	  libAtoms_misc_utils
else
  LIBATOMS_F95_FILES += libAtoms \
	  libAtoms_misc_utils
endif		

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

SVN_VERSION := -D'SVN_VERSION="$(shell svnversion -n ${BUILDPRE})"'
ARCH_F95FLAGS += ${SVN_VERSION}
CFLAGS += ${SVN_VERSION}

all: libatoms.a

libatoms.a: ${LIBATOMS_OBJS}
	ar ${AR_ADD} libatoms.a $?

clean:
	rm -f libatoms.a *.o *.mod

allclean: clean docclean


${LIBATOMS_C_OBJS}  : %.o : %.c
	${CC} ${CFLAGS} ${ARCH_CCFLAGS} -c $< 

md:  libatoms.a md.o
	$(LINKER)  -o $@ ${LINKFLAGS} $@.o  libatoms.a ${LIBDIRS} ${ARCH_LAPACK_LIBS} ${ARCH_LIBS}

nc2xyz: xyz_netcdf.c
	${CC} ${CFLAGS}  ${ARCH_CCFLAGS} $< -o $@ ${ARCH_LIBS} ${SVN_VERSION} -lm -DMAIN_PROGRAM

converters: nc2xyz
	cp $< ${HOME}/bin
	ln -sf ${HOME}/bin/nc2xyz ${HOME}/bin/nc2nc
	ln -sf ${HOME}/bin/nc2xyz ${HOME}/bin/ncstat
	ln -sf ${HOME}/bin/nc2xyz ${HOME}/bin/xyz2nc
	ln -sf ${HOME}/bin/nc2xyz ${HOME}/bin/xyz2xyz
	ln -sf ${HOME}/bin/nc2xyz ${HOME}/bin/ncstat

-include depend

libatoms.tex: ${LIBATOMS_F95_SOURCES} f90doc intro.tex md.f95
	@if (! ./checkcode ${LIBATOMS_F95_SOURCES}); then exit 1; fi
	./f90doc -t 'libAtoms Reference Manual' -a "James Kermode, G\\'abor Cs\\'anyi, Steven Winfield" -l ${LIBATOMS_F95_SOURCES}  > libatoms.tex

libatoms.dvi: libatoms.tex f90doc
	latex libatoms.tex
	latex libatoms.tex
	makeindex general
	latex libatoms.tex

libatoms.ps: libatoms.dvi
	dvips libatoms.dvi -o libatoms.ps

libatoms.pdf: libatoms.dvi
	dvipdfm libatoms.dvi

htmldoc: libatoms.tex
	latex2html libatoms.tex

doc: libatoms.ps libatoms.pdf

docclean: 
	-rm -f libatoms.dvi libatoms.tex libatoms.aux libatoms.toc libatoms.pdf libatoms.ps libatoms.log libatoms.idx libatoms.out general.i* 
	-rm -rf libatoms

