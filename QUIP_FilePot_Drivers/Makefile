# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2010.
# H0 X
# H0 X   These portions of the source code are released under the GNU General
# H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# H0 X
# H0 X   If you would like to license the source code under different terms,
# H0 X   please contact Gabor Csanyi, gabor@csanyi.net
# H0 X
# H0 X   Portions of this code were written by Noam Bernstein as part of
# H0 X   his employment for the U.S. Government, and are not subject
# H0 X   to copyright in the USA.
# H0 X
# H0 X
# H0 X   When using this software, please cite the following reference:
# H0 X
# H0 X   http://www.libatoms.org
# H0 X
# H0 X  Additional contributions by
# H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# H0 X
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.inc
include Makefile.rules

PROGRAMS = evb_driver
EVB_OBJS =

evb_driver : % :  %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${EVB_OBJS} ${LIBS}  ${LINKOPTS}

ifeq ($(HAVE_CP2K), 1)
   PROGRAMS += cp2k_driver
   CP2K_OBJS = 
cp2k_driver : % :  ${CP2K_OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${CP2K_OBJS} ${LIBS}  ${LINKOPTS}
else
cp2k_driver:
	echo "CP2K support not active" 1>&2 && exit 1
endif
ifeq ($(HAVE_VASP), 1)
   PROGRAMS += vasp_driver
   VASP_OBJS = 
vasp_driver : % :  ${VASP_OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${VASP_OBJS} ${LIBS}  ${LINKOPTS}
else
vasp_driver:
	echo "VASP support not active" 1>&2 && exit 1
endif

all: ${PROGRAMS}

LIBS =  -L. -lquiputils -lquip_core ${GP_LIB} -latoms -lthirdparty
LIBFILES = libatoms.a ${GP_LIBFILE} libquip_core.a libquiputils.a

.PHONY: clean allclean test

clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

