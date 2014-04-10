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

ifeq (${QUIP_ROOT},)
   QUIP_ROOT=${PWD}
endif
export QUIP_ROOT

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

export SCRIPT_PATH=${QUIP_ROOT}/utility_scripts/

MODULES =

ifeq (${HAVE_THIRDPARTY},1)
  MODULES += ThirdParty
endif

MODULES += libAtoms QUIP_Core QUIP_Utils QUIP_Programs QUIP_FilePot_Drivers # Tests
GAP = 

ifeq (${HAVE_GAP},1)
MODULES += GAP 
GAP += GAP/libgap_predict.a
endif

ifeq (${HAVE_GAP_FILLER},1)
MODULES += GAP-filler
endif

FOX = FoX-4.0.3
EXTRA_CLEAN_DIRS = Tools/quippy

default: ${MODULES}

EXTRA_ALL_DIRS = Tools
all: default
	@for f in ${MODULES} ${EXTRA_ALL_DIRS}; do ${MAKE} $$f/all || exit 1; done

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

ifeq (${HAVE_GAP},1)
GAP: libAtoms/libatoms.a ${FOX}
endif

ifeq (${HAVE_GAP_FILLER},1)
GAP-filler: libAtoms/libatoms.a ${FOX} GAP/libgap_predict.a ${GAP} QUIP_Core/libquip_core.a QUIP_Utils
endif

QUIP_Core: libAtoms/libatoms.a ${FOX} ${GAP}
QUIP_Util: libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a
QUIP_FilePot_Drivers: libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils 
QUIP_Programs: libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils QUIP_FilePot_Drivers
Tests: libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils

ifeq (${HAVE_GAP},1)
GAP/%: libAtoms/libatoms.a ${FOX}
	ln -sf ${PWD}/GAP/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/GAP -I${PWD} -I${PWD}/Makefiles $${targ#GAP/}
	rm ${BUILDDIR}/Makefile
endif

ifeq (${HAVE_GAP_FILLER},1)
GAP-filler/%: libAtoms/libatoms.a ${FOX} GAP/libgap_predict.a QUIP_Core/libquip_core.a QUIP_Utils
	ln -sf ${PWD}/GAP-filler/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/GAP-filler -I${PWD} -I${PWD}/Makefiles $${targ#GAP-filler/}
	rm ${BUILDDIR}/Makefile
endif

QUIP_FilePot_Drivers/%: libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils
	ln -sf ${PWD}/QUIP_FilePot_Drivers/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/QUIP_FilePot_Drivers -I${PWD} -I${PWD}/Makefiles $${targ#QUIP_FilePot_Drivers/}
	rm ${BUILDDIR}/Makefile

QUIP_Programs/%: ThirdParty libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils QUIP_FilePot_Drivers
	ln -sf ${PWD}/QUIP_Programs/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/QUIP_Programs -I${PWD} -I${PWD}/Makefiles $${targ#QUIP_Programs/}
	rm ${BUILDDIR}/Makefile

QUIP_Core/%: ThirdParty libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core
	ln -sf ${PWD}/QUIP_Core/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/QUIP_Core -I${PWD} -I${PWD}/Makefiles $${targ#QUIP_Core/}
	rm ${BUILDDIR}/Makefile

QUIP_Utils/%: ThirdParty libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils
	ln -sf ${PWD}/QUIP_Utils/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/QUIP_Utils -I${PWD} -I${PWD}/Makefiles $${targ#QUIP_Utils/}
	rm ${BUILDDIR}/Makefile

libatoms: libAtoms
libAtoms/%: libAtoms 
	ln -sf ${PWD}/libAtoms/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/libAtoms -I${PWD} -I${PWD}/Makefiles $${targ#libAtoms/}
	rm ${BUILDDIR}/Makefile

Tools/%: libAtoms/libatoms.a ${FOX} ${GAP} QUIP_Core/libquip_core.a QUIP_Utils/libquiputils.a Tools
	ln -sf ${PWD}/Tools/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/Tools -I${PWD} -I${PWD}/Makefiles $${targ#Tools/}
	rm ${BUILDDIR}/Makefile

ifeq (${HAVE_THIRDPARTY},1)
ThirdParty/%: 
	ln -sf ${PWD}/ThirdParty/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/ThirdParty -I${PWD} -I${PWD}/Makefiles $${targ#ThirdParty/}
	rm ${BUILDDIR}/Makefile
else
ThirdParty:
	@echo "Placeholder ThirdParty rule"
endif

${BUILDDIR}: arch
	@if [ ! -d build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ] ; then mkdir build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX} ; fi


clean: ${BUILDDIR}
	for mods in  ${MODULES} ; do \
	  echo "clean in $$mods"; \
	  ln -sf ${PWD}/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	  ${MAKE} -C ${BUILDDIR} USE_MAKEDEP=0 QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/$$mods -I${PWD} -I${PWD}/Makefiles clean ; \
	done

deepclean: clean
	-for mods in  ${MODULES} ; do \
	  echo "deepclean in $$mods"; \
	  ln -sf ${PWD}/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	  ${MAKE} -C ${BUILDDIR} USE_MAKEDEP=0 QUIP_ROOT=${QUIP_ROOT} VPATH=${PWD}/$$mods -I${PWD} -I${PWD}/Makefiles deepclean ; \
	done
	-for dir in ${EXTRA_CLEAN_DIRS}; do \
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
	    if [ build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f -nt ${QUIP_INSTDIR}/$$f ]; then \
	       echo "copying f $$f to ${QUIP_INSTDIR}"; \
	       cp build.${QUIP_ARCH}${QUIP_ARCH_SUFFIX}/$$f ${QUIP_INSTDIR}; \
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

doc: quip-reference-manual.pdf

quip-reference-manual.pdf:
	./Tools/mkdoc

quippy:
	make -C Tools/quippy install QUIP_ROOT=${QUIP_ROOT}

test:
	${MAKE} -C Tests -I${PWD} -I${PWD}/Makefiles -I${PWD}/${BUILDDIR}

GIT_SUBDIRS=GAP GAP-filler

git_pull_all:
	git pull
	@for d in ${GIT_SUBDIRS}; do if [ -d $$d ]; then pushd $$d; git pull; popd; fi; done

distribution:
	/bin/bash ./utility_scripts/gitversion > GIT_VERSION
	/bin/bash ./utility_scripts/gapversion.sh > GAP_VERSION
	git archive HEAD > ../QUIP.distribution.`date +%Y-%m-%d`.tar
	tar rvf ../QUIP.distribution.`date +%Y-%m-%d`.tar GIT_VERSION GAP_VERSION
	bzip2 ../QUIP.distribution.`date +%Y-%m-%d`.tar
	rm GIT_VERSION GAP_VERSION
	@for d in ${GIT_SUBDIRS}; do if [ -d $$d ]; then pushd $$d; git archive HEAD | bzip2 > ../../$$d.distribution.`date +%Y-%m-%d`.tar.bz2; popd; fi; done
