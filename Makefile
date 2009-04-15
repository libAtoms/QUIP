
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

${FOX}: ${FOX}/objs/lib/libFoX_common.a
${FOX}/objs/lib/libFoX_common.a:
	make -C ${FOX} -I${PWD}/Makefiles -I${PWD}/${BUILDDIR} -f Makefile.QUIP


${MODULES}: ${BUILDDIR} ${BUILDDIR}/Makefile.inc
	ln -sf ${PWD}/$@/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} VPATH=${PWD}/$@ -I${PWD}/Makefiles
	rm ${BUILDDIR}/Makefile

QUIP_Core: libAtoms ${FOX}
QUIP_Util: libAtoms ${FOX} QUIP_Core
QUIP_Programs: libAtoms ${FOX} QUIP_Core QUIP_Utils 
Tests: libAtoms ${FOX} QUIP_Core QUIP_Utils

QUIP_Programs/%: libAtoms ${FOX} QUIP_Core QUIP_Utils
	ln -sf ${PWD}/QUIP_Programs/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} VPATH=${PWD}/QUIP_Programs -I${PWD}/Makefiles $${targ#QUIP_Programs/}
	rm ${BUILDDIR}/Makefile

libAtoms/%: libAtoms 
	ln -sf ${PWD}/libAtoms/Makefile ${BUILDDIR}/Makefile
	targ=$@ ; ${MAKE} -C ${BUILDDIR} VPATH=${PWD}/libAtoms -I${PWD}/Makefiles $${targ#libAtoms/}
	rm ${BUILDDIR}/Makefile


${BUILDDIR}: arch
	@if [ ! -d build.${QUIP_ARCH} ] ; then mkdir build.${QUIP_ARCH} ; fi

${BUILDDIR}/Makefile.inc:
	-rm -f ${BUILDDIR}/Makefile.inc
ifndef MATH_LINKOPTS
	@echo ; echo "In the following, hit enter to accept the defaults."
	@echo
	@echo "Please enter the linking options for LAPACK and BLAS libraries:"
	@echo "   Default: ${DEFAULT_MATH_LINKOPTS}"
	@read MATH_LINKOPTS && if [[ $$MATH_LINKOPTS == "" ]] ; then \
	echo "MATH_LINKOPTS=${DEFAULT_MATH_LINKOPTS}" >> ${BUILDDIR}/Makefile.inc ; \
	else echo "MATH_LINKOPTS=$$MATH_LINKOPTS" >> ${BUILDDIR}/Makefile.inc ; fi
endif
ifndef FOX_LIBDIR
	@echo ; \
        echo "Please enter directory where FoX libraries are kept:" ; \
	echo "   Default: use included version ${FOX}" ; \
        read FOX_LIBDIR && if [[ $$FOX_LIBDIR ]] ; then \
	  echo "FOX_LIBDIR=$$FOX_LIBDIR" >> ${BUILDDIR}/Makefile.inc ; echo ; \
	  echo "Please enter directory where FoX include files are kept:" ; \
	  read FOX_INCDIR && echo "FOX_INCDIR=$$FOX_INCDIR" >> ${BUILDDIR}/Makefile.inc ; \
	  echo "HAVE_EXTERNAL_FOX=1" >> ${BUILDDIR}/Makefile.inc ; \
	else echo "FOX_LIBDIR=$${PWD}/FoX-4.0.3/objs/lib" >> ${BUILDDIR}/Makefile.inc; \
	  echo "FOX_INCDIR=$${PWD}/FoX-4.0.3/objs/finclude" >> ${BUILDDIR}/Makefile.inc; \
	  echo "HAVE_EXTERNAL_FOX=0" >> ${BUILDDIR}/Makefile.inc ; fi
endif
ifndef NETCDF_LIBDIR
	@echo ; \
        echo "Please enter directory where NetCDF libraries are kept:" ; \
	echo "   Default: no NetCDF present" ; \
        read NETCDF_LIBDIR && if [[ $$NETCDF_LIBDIR ]] ; then \
	  echo "NETCDF_LIBDIR=$$NETCDF_LIBDIR" >> ${BUILDDIR}/Makefile.inc ; echo ; \
	  echo "Please enter directory where NetCDF include files are kept:" ; \
	  read NETCDF_INCDIR && echo "NETCDF_INCDIR=$$NETCDF_INCDIR" >> ${BUILDDIR}/Makefile.inc ; \
	  echo "HAVE_NETCDF=1" >> ${BUILDDIR}/Makefile.inc ; \
	  if nm $$NETCDF_LIBDIR/libnetcdf.a | grep -q deflate; then \
            echo "NetCDF version 4 found."; echo "NETCDF4=1" >> ${BUILDDIR}/Makefile.inc; \
	  else echo "NetCDF older than version 4 found."; echo "NETCDF4=0" >> ${BUILDDIR}/Makefile.inc; fi \
	else echo "HAVE_NETCDF=0" >> ${BUILDDIR}/Makefile.inc ; fi
endif
ifndef LARSPOT_LIBDIR
	@echo ; \
        echo "Please enter directory where the Lars Potential libraries are kept:" ; \
	echo "   Default: no Lars potential present" ; \
        read LARSPOT_LIBDIR && if [[ $$LARSPOT_LIBDIR ]] ; then \
	  echo "LARSPOT_LIBDIR=$$LARSPOT_LIBDIR" >> ${BUILDDIR}/Makefile.inc ; echo ; \
	  echo "Please enter directory where Lars Potential include files are kept:" ; \
	  read LARSPOT_INCDIR && echo "LARSPOT_INCDIR=$$LARSPOT_INCDIR" >> ${BUILDDIR}/Makefile.inc ; \
	  echo "HAVE_LARSPOT=1" >> ${BUILDDIR}/Makefile.inc ; \
	else echo "HAVE_LARSPOT=0" >> ${BUILDDIR}/Makefile.inc ; fi
endif
ifndef EXTRA_LINKOPTS
	@echo
	@echo "Please enter any other extra inking options:"
	@echo "   Default: none"
	@read EXTRA_LINKOPTS && echo "EXTRA_LINKOPTS=$$EXTRA_LINKOPTS" >> ${BUILDDIR}/Makefile.inc
endif
	echo "HAVE_LOTF=0" >> ${BUILDDIR}/Makefile.inc
	echo "OPTIM=${DEFAULT_OPTIM}" >> ${BUILDDIR}/Makefile.inc
	echo "DEBUG=${DEFAULT_DEBUG}" >> ${BUILDDIR}/Makefile.inc


clean:
	for mods in  ${MODULES}  ; do \
	ln -sf ${PWD}/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	${MAKE} -C ${BUILDDIR} -I${PWD}/Makefiles clean ; \
	done


doc: quip-reference-manual.pdf

quip-reference-manual.pdf:
	./Tools/mkdoc

