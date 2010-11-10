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

ifneq (${QUIP_ARCH},)
	export BUILDDIR=build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX}
	export QUIP_ARCH
	include Makefile.config
	include Makefiles/Makefile.${QUIP_ARCH}
	include ${BUILDDIR}/Makefile.inc
	include Makefile.rules
else
	BUILDDIR=crap
endif

ifeq (${QUIP_ROOT},)
   QUIP_ROOT=${PWD}
endif

export SCRIPT_PATH=${QUIP_ROOT}/utility_scripts/

MODULES = libAtoms QUIP_Core QUIP_Utils QUIP_Programs # Tests
GP = 

ifeq (${HAVE_GP_PREDICT},1)
MODULES += GAP_predict 
GP += GAP_predict
endif

ifeq (${HAVE_GP_TEACH},1)
MODULES += GAP_teach GAProgs 
GP += GAP_teach
endif

FOX = FoX-4.0.3
EXTRA_CLEAN_DIRS = Tools/quippy

default: ${MODULES}

EXTRA_ALL_DIRS = Tools
all: default
	@for f in ${MODULES} ${EXTRA_ALL_DIRS}; do ${MAKE} $$f/all; done

.PHONY: arch ${MODULES} doc clean install test quippy doc install-structures install-dtds install-Tools install-build.QUIP_ARCH

arch: 
ifeq (${QUIP_ARCH},)
	@echo
	@echo "You need to define the architecture using the QUIP_ARCH variable"
	@echo
	@exit 1
endif

${FOX}: ${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a
${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a:
	make -C ${FOX} -I${PWD} -I${PWD}/Makefiles -I${PWD}/${BUILDDIR} -f Makefile.QUIP 


${MODULES}:  ${BUILDDIR}
	ln -sf ${PWD}/$@/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/$@ -I${PWD} -I${PWD}/Makefiles
	rm ${BUILDDIR}/Makefile

ifeq (${HAVE_GP_PREDICT},1)
GAP_predict: libAtoms ${FOX}
endif

ifeq (${HAVE_GP_TEACH},1)
GAP_teach: libAtoms ${FOX} GAP_predict 
GAProgs: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils
endif

QUIP_Core: libAtoms ${FOX} ${GP}
QUIP_Util: libAtoms ${FOX} ${GP} QUIP_Core
QUIP_Programs: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils 
Tests: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils

ifeq (${HAVE_GP_PREDICT},1)
GAP_predict/%: libAtoms ${FOX}
	ln -sf ${PWD}/GAP_predict/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/GAP_predict -I${PWD} -I${PWD}/Makefiles $${targ#GAP_predict/}
	rm ${BUILDDIR}/Makefile
endif

ifeq (${HAVE_GP_TEACH},1)
GAP_teach/%: libAtoms ${FOX} GAP_predict
	ln -sf ${PWD}/GAP_teach/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/GAP_teach -I${PWD} -I${PWD}/Makefiles $${targ#GAP_teach/}
	rm ${BUILDDIR}/Makefile

GAProgs/%: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/GAProgs/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/GAProgs -I${PWD} -I${PWD}/Makefiles $${targ#GAProgs/}
	rm ${BUILDDIR}/Makefile
endif

QUIP_Programs/%: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/QUIP_Programs/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/QUIP_Programs -I${PWD} -I${PWD}/Makefiles $${targ#QUIP_Programs/}
	rm ${BUILDDIR}/Makefile

QUIP_Core/%: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/QUIP_Core/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/QUIP_Core -I${PWD} -I${PWD}/Makefiles $${targ#QUIP_Core/}
	rm ${BUILDDIR}/Makefile

libAtoms/%: libAtoms 
	ln -sf ${PWD}/libAtoms/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/libAtoms -I${PWD} -I${PWD}/Makefiles $${targ#libAtoms/}
	rm ${BUILDDIR}/Makefile

Tools/%: libAtoms ${FOX} ${GP} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/Tools/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/Tools -I${PWD} -I${PWD}/Makefiles $${targ#Tools/}
	rm ${BUILDDIR}/Makefile


${BUILDDIR}: arch
	@if [ ! -d build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ] ; then mkdir build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ; fi


clean: ${BUILDDIR}
	for mods in  ${MODULES} ; do \
	  ln -sf ${PWD}/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	  ${MAKE} -C ${BUILDDIR} USE_MAKEDEP=0 QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/$$mods -I${PWD} -I${PWD}/Makefiles clean ; \
	done ; \

deepclean: clean
	for dir in ${EXTRA_CLEAN_DIRS}; do \
	  cd $$dir; make clean; \
	done

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
	@echo "installing from build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX}"; \
	for f in `/bin/ls build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX} | egrep -v '\.o|\.a|\.mod|Makefile*|^test$$'`; do \
	  if [ -x build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f ]; then \
	    echo "copying f $$f to ${QUIP_INSTDIR}"; \
	    cp build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f ${QUIP_INSTDIR}; \
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

doc: quip-reference-manual.pdf

quip-reference-manual.pdf:
	./Tools/mkdoc

quippy:
	make -C Tools/quippy install QUIP_ROOT=${QUIP_ROOT}

test:
	${MAKE} -C Tests -I${PWD} -I${PWD}/Makefiles -I${PWD}/${BUILDDIR}
