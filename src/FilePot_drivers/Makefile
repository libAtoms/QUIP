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

include Makefile.${QUIP_ARCH}
include Makefile.inc
include Makefile.rules

default: all

LIBS =  -L. -lquiputils -lquip_core ${GP_LIB} -latoms 
ifeq (${HAVE_THIRDPARTY},1)
LIBS += -lthirdparty
endif
LIBFILES = libatoms.a ${GP_LIBFILE} libquip_core.a libquiputils.a

PROGRAMS = evb_driver

evb_driver : % :  %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${LIBS}  ${LINKOPTS}

ifeq ($(HAVE_CP2K), 1)
   PROGRAMS += cp2k_driver
   CP2K_LIB_OBJS = cp2k_driver_module.o
   CLEAN_OBJS += ${CP2K_LIB_OBJS}
   CLEAN_LIBS += libcp2k_driver.a
libcp2k_driver.a: ${CP2K_LIB_OBJS}
ifneq (${LIBTOOL},)
	${LIBTOOL} -o libcp2k_driver.a ${CP2K_LIB_OBJS}
else
	${AR} ${AR_ADD} libcp2k_driver.a $?
endif
cp2k_driver : % :  libcp2k_driver.a %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o -lcp2k_driver ${LIBS} ${LINKOPTS}
else
cp2k_driver:
	echo "CP2K support not active" 1>&2 && exit 1
endif

ifeq ($(HAVE_VASP), 1)
   PROGRAMS += vasp_driver
vasp_driver : % :  %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${LIBS}  ${LINKOPTS}
else
vasp_driver:
	echo "VASP support not active" 1>&2 && exit 1
endif

.PHONY: clean allclean test install

ifeq ($(HAVE_CP2K), 1)
all: libcp2k_driver.a ${PROGRAMS}
else
all: ${PROGRAMS}
endif

install:
	@if [ ! -d ${QUIP_INSTALLDIR} ]; then \
	  echo "make install: QUIP_INSTALLDIR '${QUIP_INSTALLDIR}' doesn't exist or isn't a directory"; \
	  exit 1; \
	else	\
	  for f in ${PROGRAMS} ; do \
	    echo "Copying $$f to ${QUIP_INSTALLDIR}" ; \
	    cp $$f ${QUIP_INSTALLDIR} ; \
	  done ;\
	fi

clean:	
	-rm -f ${CLEAN_OBJS} ${CLEAN_LIBS} ${PROGRAMS} ${addsuffix .o, ${PROGRAMS}}

allclean : clean

