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

ifeq (${QUIP_ROOT},)
   QUIP_ROOT=${PWD}
endif
export QUIP_ROOT

ifneq (${QUIP_ARCH},)
   export BUILDDIR=build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX}
   export QUIP_ARCH
   $(info Configuring with QUIP_ARCH=$(QUIP_ARCH))
   include arch/Makefile.${QUIP_ARCH}
   include Makefile.rules
   ifneq ("$(wildcard $(BUILDDIR)/Makefile.inc)","")
      -include ${BUILDDIR}/Makefile.inc   
   endif
   include Makefile.config
else
   BUILDDIR=dummy
endif

export SCRIPT_PATH=${QUIP_ROOT}/bin

MODULES =

ifeq (${HAVE_THIRDPARTY},1)
   MODULES += ThirdParty
   THIRDPARTY_LIBS := libthirdparty.a
ifeq (${HAVE_FX},1)
   THIRDPARTY_LIBS += libfx.a
endif
endif

MODULES += libAtoms Potentials Utils Programs FilePot_drivers Structure_processors 

ifeq (${HAVE_GAP},1)
MODULES += GAP 
GAP += GAP/libgap_predict.a
else
GAP = 
endif

ifeq (${HAVE_GAP_FILLER},1)
MODULES += GAP-filler
endif

FOX = FoX-4.0.3
EXTRA_CLEAN_DIRS = quippy

.DEFAULT_GOAL=all

default: ${MODULES}

all: default
	@for f in ${MODULES} ; do ${MAKE} $$f/all || exit 1; done



.PHONY: arch ${MODULES} config doc clean deepclean distclean install test quippy doc install-structures install-dtds install-Tools install-build.QUIP_ARCH libquip

arch: 
ifeq (${QUIP_ARCH},)
	@echo
	@echo "You need to define the architecture using the QUIP_ARCH variable"
	@echo "Check out the `arch' subdirectory for example architectures."
	@echo
	@exit 1
endif


${BUILDDIR}/Makefile.inc: 
	@echo
	@echo "Perhaps you forgot to run \`make config'?"
	@echo
	@exit 1


${FOX}: ${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a
${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a:
	make -C ${FOX} -I${PWD} -I${PWD}/Makefiles -I${PWD}/${BUILDDIR} -f Makefile.QUIP 


FOX_STATIC_LIBFILES = $(patsubst -l%,${FOX_LIBDIR}/lib%.a,${FOX_LIBS})
FOX_STATIC_LIBFILE_OBJS = $(shell for i in ${FOX_STATIC_LIBFILES}; do ar -t $$i; done | grep \.o)
# general rule to make a module

${MODULES}:  ${BUILDDIR}/Makefile.inc ${BUILDDIR}
	rm -f ${BUILDDIR}/Makefile
	cp ${PWD}/src/$@/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/$@ -I${PWD} -I${PWD}/arch
	rm ${BUILDDIR}/Makefile

# dependencies between modules

ifeq (${HAVE_GAP},1)
GAP: libAtoms ${FOX}
endif

ifeq (${HAVE_GAP_FILLER},1)
GAP-filler: libAtoms/libatoms.a ${FOX} GAP ${GAP} Potentials Utils
endif

Potentials: libAtoms ${FOX} ${GAP}
Utils: libAtoms ${FOX} ${GAP} Potentials
FilePot_drivers: libAtoms ${FOX} ${GAP} Potentials Utils
Programs: libAtoms ${FOX} ${GAP} Potentials Utils FilePot_drivers
Tests: libAtoms ${FOX} ${GAP} Potentials Utils
libatoms: libAtoms

libquip: libquip.a

libquip.a: ThirdParty libAtoms ${FOX} ${GAP} Potentials Utils
	LIBQUIP_OBJS="$(shell for i in ${BUILDDIR}/libquiputils.a ${BUILDDIR}/libquip_core.a $(subst GAP,${BUILDDIR},${GAP}) ${BUILDDIR}/libatoms.a $(addprefix ${BUILDDIR}/,${THIRDPARTY_LIBS}) ${FOX_STATIC_LIBFILES}; do ar -t $$i; done | grep \.o)" && \
		     cd ${BUILDDIR} && for i in ${FOX_STATIC_LIBFILES}; do ar -x $$i; done && ar -rcs $@ $$LIBQUIP_OBJS

${BUILDDIR}: arch
	@if [ ! -d build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ] ; then mkdir -p build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ; fi


clean: ${BUILDDIR}
	for mods in  ${MODULES} ; do \
	  echo "clean in $$mods"; \
	  rm -f ${BUILDDIR}/Makefile ; \
	  cp ${PWD}/src/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	  ${MAKE} -C ${BUILDDIR} USE_MAKEDEP=0 QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/$$mods -I${PWD} -I${PWD}/Makefiles clean ; \
	done

deepclean: clean
	-for mods in  ${MODULES} ; do \
	  echo "deepclean in $$mods"; \
          rm -f ${BUILDDIR}/Makefile ; \
	  cp ${PWD}/src/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	  ${MAKE} -C ${BUILDDIR} USE_MAKEDEP=0 QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/src/$$mods -I${PWD} -I${PWD}/Makefiles deepclean ; \
	done
	-for dir in ${EXTRA_CLEAN_DIRS}; do \
	  cd $$dir; make clean; \
	done
	-if [[ -d ${FOX}/objs.${QUIP_ARCH} ]]; then \
	  rm -rf ${FOX}/objs.${QUIP_ARCH} ; \
	fi

distclean: deepclean
	rm -rf ${BUILDDIR}

install:
	@if [ "x${QUIP_INSTDIR}" == "x" ]; then \
	  echo "make install needs QUIP_INSTDIR defined"; \
	  exit 1; \
	fi; \
	if [ ! -d ${QUIP_INSTDIR} ]; then \
	  echo "make install QUIP_INSTDIR '${QUIP_INSTDIR}' doesn't exist or isn't a directory"; \
	  exit 1; \
	fi
	${MAKE} install-build.QUIP_ARCH install-Tools install-structures install-dtds

install-build.QUIP_ARCH:
	@echo "installing from build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX}"; \
	for f in `/bin/ls build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX} | egrep -v '\.o|\.a|\.mod|Makefile*|^test$$'`; do \
	  if [ -x build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f ]; then \
	    if [ build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f -nt ${QUIP_INSTDIR}/$$f ]; then \
	       echo "copying f $$f to ${QUIP_INSTDIR}"; \
	       cp build/${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f ${QUIP_INSTDIR}; \
	    fi; \
	    if [ $$f == eval ]; then \
	       rm -f ${QUIP_INSTDIR}/quip_eval; \
	       ln -s ${QUIP_INSTDIR}/eval ${QUIP_INSTDIR}/quip_eval; \
	    fi; \
	  fi; \
	done

install-Tools:
	@echo "installing from Tools"; \
	for f in `/bin/ls Tools | egrep -v '\.f95|Makefile*'`; do \
	  if [ -f Tools/$$f ] && [ -x Tools/$$f ]; then \
	    echo "copying f $$f to ${QUIP_INSTDIR}"; \
	    cp Tools/$$f ${QUIP_INSTDIR}; \
	  fi; \
	done

install-structures:
	cd structures; ${MAKE} QUIP_STRUCTS_DIR=$(QUIP_STRUCTS_DIR) install-structures 

install-dtds:
	cd dtds; ${MAKE} QUIP_STRUCTS_DIR=$(QUIP_STRUCTS_DIR) install-dtds 

quippy:
	make -C quippy install QUIP_ROOT=${QUIP_ROOT}

test:
	${MAKE} -C Tests -I${PWD} -I${PWD}/arch -I${PWD}/${BUILDDIR}

GIT_SUBDIRS=src/GAP src/GAP-filler src/ThirdParty

git_pull_all:
	git pull
	@for d in ${GIT_SUBDIRS}; do if [ -d $$d ]; then pushd $$d; git pull; popd; fi; done

distribution:
	./utility_scripts/gitversion > GIT_VERSION
	./utility_scripts/gapversion.sh > GAP_VERSION
	git archive HEAD > ../QUIP.distribution.`date +%Y-%m-%d`.tar
	tar rvf ../QUIP.distribution.`date +%Y-%m-%d`.tar GIT_VERSION GAP_VERSION
	bzip2 ../QUIP.distribution.`date +%Y-%m-%d`.tar
	rm GIT_VERSION GAP_VERSION
	@for d in ${GIT_SUBDIRS}; do if [ -d $$d ]; then pushd $$d; git archive HEAD | bzip2 > ../../$$d.distribution.`date +%Y-%m-%d`.tar.bz2; popd; fi; done



