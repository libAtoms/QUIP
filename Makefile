
ifneq (${QUIP_ARCH},)
	export BUILDDIR=build.${QUIP_ARCH}
	export QUIP_ARCH
	include Makefiles/Makefile.${QUIP_ARCH}
else
	BUILDDIR=crap
endif


FOX = FoX-4.0.3
MODULES = libAtoms QUIP_Core QUIP_Utils QUIP_Programs # Tests

all: ${MODULES}

.PHONY: arch ${MODULES} doc

arch: 
ifeq (${QUIP_ARCH},)
	@echo
	@echo "You need to define the architecture using the QUIP_ARCH variable"
	@echo
	@exit 1
endif

${FOX}: ${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a
${FOX}/objs.${QUIP_ARCH}/lib/libFoX_common.a:
	make -C ${FOX} -I${PWD}/Makefiles -I${PWD}/${BUILDDIR} -f Makefile.QUIP 


${MODULES}: ${BUILDDIR}
	ln -sf ${PWD}/$@/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} VPATH=${PWD}/$@ -I${PWD}/Makefiles
	rm ${BUILDDIR}/Makefile

QUIP_Core: libAtoms ${FOX}
QUIP_Util: libAtoms ${FOX} QUIP_Core
QUIP_Programs: libAtoms ${FOX} QUIP_Core QUIP_Utils 
Tests: libAtoms ${FOX} QUIP_Core QUIP_Utils

QUIP_Programs/Examples/%: libAtoms ${FOX} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/QUIP_Programs/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} VPATH=${PWD}/QUIP_Programs/Examples -I${PWD}/Makefiles $${targ#QUIP_Programs/Examples/}
	rm ${BUILDDIR}/Makefile

QUIP_Programs/%: libAtoms ${FOX} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/QUIP_Programs/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} VPATH=${PWD}/QUIP_Programs -I${PWD}/Makefiles $${targ#QUIP_Programs/}
	rm ${BUILDDIR}/Makefile

libAtoms/%: libAtoms 
	ln -sf ${PWD}/libAtoms/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} VPATH=${PWD}/libAtoms -I${PWD}/Makefiles $${targ#libAtoms/}
	rm ${BUILDDIR}/Makefile

Tools/%: libAtoms ${FOX} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/Tools/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} VPATH=${PWD}/Tools -I${PWD}/Makefiles $${targ#Tools/}
	rm ${BUILDDIR}/Makefile


${BUILDDIR}: arch
	@if [ ! -d build.${QUIP_ARCH} ] ; then mkdir build.${QUIP_ARCH} ; fi


clean:
	for mods in  ${MODULES}  ; do \
	ln -sf ${PWD}/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	${MAKE} -C ${BUILDDIR} -I${PWD}/Makefiles clean ; \
	done


doc: quip-reference-manual.pdf

quip-reference-manual.pdf:
	./Tools/mkdoc

atomeye:
	if [[ ! -d Tools/AtomEye ]]; then svn co svn+ssh://cvs.tcm.phy.cam.ac.uk/home/jrk33/repo/trunk/AtomEye Tools/AtomEye; fi
	make -C Tools/AtomEye QUIP_ROOT=${PWD}

quippy:
	make -C Tools/quippy install QUIP_ROOT=${PWD}

test:
	${MAKE} -C Tests -I${PWD}/Makefiles -I${PWD}/${BUILDDIR}
