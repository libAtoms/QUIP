# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2015.
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


.PHONY: config doc clean deepclean distclean install test quippy doc install-structures install-dtds install-Tools install-build.QUIP_ARCH install-quippy libquip

ifndef QUIP_ARCH
$(error "You need to define the architecture using the QUIP_ARCH variable. Check out the arch/ subdirectory.")
else
include arch/Makefile.${QUIP_ARCH}
endif

# include other makefiles and export env variables
export QUIP_ARCH

ifeq (${QUIP_ROOT},)
QUIP_ROOT=${PWD}
endif

export QUIP_ROOT
export SCRIPT_PATH=${QUIP_ROOT}/bin
export BUILDDIR=${QUIP_ROOT}/build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX}

-include ${BUILDDIR}/Makefile.inc

# create modules list

MODULES=
# add any third party packages first
ifeq (${HAVE_THIRDPARTY},1)
   THIRDPARTY = ThirdParty
   MODULES += ThirdParty
   THIRDPARTY_LIBS := libthirdparty.a
ifeq (${HAVE_FX},1)
   THIRDPARTY_LIBS += libfx.a
endif
ifeq (${HAVE_SCME},1)
   THIRDPARTY_LIBS += libscme.a
endif
ifeq (${HAVE_MTP},1)
   THIRDPARTY_LIBS += libmtp.a
endif
endif

MODULES += libAtoms

# add GAP modules if we have them - they need to come before other modules, except for libAtoms
ifeq (${HAVE_GAP},1)
MODULES += GAP
GAP += GAP
GAP_LIBS += libgap_predict.a
else
GAP =
endif

ifeq (${HAVE_GAP_FILLER},1)
MODULES += GAP-filler
endif

# now add the rest of the modules
MODULES += Potentials Utils Programs FilePot_drivers Structure_processors


# diagnostic
$(info Using QUIP_ARCH=${QUIP_ARCH}, MODULES=${MODULES}, QUIP_ROOT=${QUIP_ROOT})

# the first target
all: ${MODULES}

FOX = fox
export FOX_LIBDIR=${QUIP_ROOT}/src/${FOX}/objs.${QUIP_ARCH}/lib
export FOX_INCDIR=${QUIP_ROOT}/src/${FOX}/objs.${QUIP_ARCH}/finclude

# now we can include the config makefile, it needs to come after the default target
include Makefile.config
include Makefile.rules

${BUILDDIR}/Makefile.inc:
	@if [ "$(MAKECMDGOALS)" != config ]; then\
		echo ;\
		echo "${BUILDDIR}/Makefile.inc not found. Perhaps you forgot to run \`make config'?" ;\
		echo ;\
		exit 1 ;\
		fi

# Automatically pull the submodules if the user didn't.
# Remove the Makefile.QUIP if it has failed previously.
src/${FOX}/configure:
	@echo "Attempting to automatically clone fox submodule"
	rm -f src/${FOX}/Makefile.QUIP
	git submodule update --init src/${FOX} || \
	    { echo "fox clone failed. Download it manually" ; exit 1 ; }

${FOX}: src/${FOX}/configure src/${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a
src/${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a:
	cp Makefile.fox src/${FOX}/Makefile.QUIP
	make -C src/${FOX} -I${PWD} -I${PWD}/arch -I${BUILDDIR} -f Makefile.QUIP

FOX_STATIC_LIBFILES = $(patsubst -l%,${FOX_LIBDIR}/lib%.a,${FOX_LIBS})
FOX_STATIC_LIBFILE_OBJS = $(shell for i in ${FOX_STATIC_LIBFILES}; do ar -t $$i; done | grep \.o)

# general rule to make a module

${MODULES}:  ${BUILDDIR}/Makefile.inc ${BUILDDIR} ${FOX}
	rm -f ${BUILDDIR}/Makefile
	cp ${PWD}/src/$@/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/$@ -I${PWD} -I${PWD}/arch
	rm ${BUILDDIR}/Makefile

# general rule to make a program in the Programs
# src directory, makes sure everything else is
# built first

Programs/% src/Programs/% : ${MODULES}
	rm -f ${BUILDDIR}/Makefile
	cp ${PWD}/src/Programs/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/Programs -I${PWD} -I${PWD}/arch $(lastword $(subst /, ,$@))
	rm ${BUILDDIR}/Makefile

# dependencies between modules

ifeq (${HAVE_GAP},1)
GAP: libAtoms ${FOX}
endif

ifeq (${HAVE_GAP_FILLER},1)
GAP-filler: libAtoms GAP Potentials Utils
endif

Potentials: libAtoms  ${GAP}
Utils:  libAtoms ${GAP} Potentials
FilePot_drivers:  libAtoms  Potentials Utils
Programs: libAtoms ${GAP} Potentials Utils
Tests: libAtoms  ${GAP} Potentials Utils
libatoms: libAtoms

libquip: libquip.a

libquip.a: ${MODULES}
	LIBQUIP_OBJS="$(shell for i in ${BUILDDIR}/libquiputils.a ${BUILDDIR}/libquip_core.a $(addprefix ${BUILDDIR}/,${GAP_LIBS}) ${BUILDDIR}/libatoms.a $(addprefix ${BUILDDIR}/,${THIRDPARTY_LIBS}) ${FOX_STATIC_LIBFILES}; do ar -t $$i; done | grep \.o)" && \
		     cd ${BUILDDIR} && for i in ${FOX_STATIC_LIBFILES}; do ar -x $$i; done && ar -rcs $@ $$LIBQUIP_OBJS

${BUILDDIR}:
	@if [ ! -d build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ] ; then mkdir -p build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ; fi

quippy: libquip.a
	rm -f ${BUILDDIR}/Makefile
	cp ${PWD}/quippy/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/Programs -I${PWD} -I${PWD}/arch build
	rm ${BUILDDIR}/Makefile

clean: ${BUILDDIR}
	-${MAKE} clean-quippy
	for mods in ${MODULES} ; do \
	  echo "clean in $$mods"; \
	  rm -f ${BUILDDIR}/Makefile ; \
	  cp ${PWD}/src/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	  ${MAKE} -C ${BUILDDIR} USE_MAKEDEP=0 QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/$$mods -I${PWD} -I${PWD}/arch clean ; \
	done
	rm -f ${BUILDDIR}/libquip.a
	rm -rf src/${FOX}/objs.${QUIP_ARCH}

deepclean: clean

distclean: clean
	rm -rf build

install-structures:
ifeq (${QUIP_STRUCTS_DIR},)
	@echo
	@echo "QUIP_STRUCTS_DIR must be defined to install structures"
else
	${MAKE} -C share/Structures QUIP_STRUCTS_DIR=$(QUIP_STRUCTS_DIR) install
endif

install: ${MODULES} install-structures
ifeq (${QUIP_INSTALLDIR},)
	@echo
	@echo "'make install' needs QUIP_INSTALLDIR to be defined to install "
	@echo "programs"
else
	@if [ ! -d ${QUIP_INSTALLDIR} ]; then \
	  echo "make install: QUIP_INSTALLDIR '${QUIP_INSTALLDIR}' doesn't exist or isn't a directory"; \
	  exit 1; \
	else	 \
	  for mods in ${MODULES} ; do \
	    rm -f ${BUILDDIR}/Makefile ;\
	    cp ${PWD}/src/$$mods/Makefile ${BUILDDIR}/Makefile ;\
	    ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} -I${PWD} -I${PWD}/arch install ;\
	    rm ${BUILDDIR}/Makefile ;\
	  done ;\
	fi
endif

test: quippy
	${MAKE} -C tests -I${PWD} -I${PWD}/arch -I${BUILDDIR}

GIT_SUBDIRS=src/GAP src/GAP-filler src/ThirdParty

git_pull_all:
	git pull
	@for d in ${GIT_SUBDIRS}; do if [ -d $$d ]; then pushd $$d; git pull; popd; fi; done

distribution:
	./bin/gitversion > GIT_VERSION
	./bin/gapversion.sh > GAP_VERSION
	git archive HEAD > ../QUIP.distribution.`date +%Y-%m-%d`.tar
	tar rvf ../QUIP.distribution.`date +%Y-%m-%d`.tar GIT_VERSION GAP_VERSION
	bzip2 ../QUIP.distribution.`date +%Y-%m-%d`.tar
	rm GIT_VERSION GAP_VERSION
	@for d in ${GIT_SUBDIRS}; do if [ -d $$d ]; then pushd $$d; git archive HEAD | bzip2 > ../../$$d.distribution.`date +%Y-%m-%d`.tar.bz2; popd; fi; done



