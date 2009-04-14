
ifneq (${ARCH},)
	export BUILDDIR=build.${ARCH}
	export ARCH
	include Makefiles/Makefile.${ARCH}
else
	BUILDDIR=crap
endif

MODULES = libAtoms QUIP_Core QUIP_Utils QUIP_Programs # Tests

all: ${MODULES}

${MODULES}: ${BUILDDIR} ${BUILDDIR}/Makefile.inc
	ln -sf ${PWD}/$@/Makefile ${BUILDDIR}/Makefile
	${MAKE} -C ${BUILDDIR} VPATH=${PWD}/$@ -I${PWD}/Makefiles
	rm ${BUILDDIR}/Makefile

QUIP_Core: libAtoms
QUIP_Util: libAtoms QUIP_Core
QUIP_Programs: libAtoms QUIP_Core QUIP_Utils 
Tests: libAtoms QUIP_Core QUIP_Utils


${BUILDDIR}:
ifneq (${ARCH},)
	if [ ! -d build.${ARCH} ] ; then mkdir build.${ARCH} ; fi
else
	@echo
	@echo "You need to define the architecture using the ARCH variable"
	@echo
	@exit 1
endif

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
	echo "   Default: no FoX present" ; \
        read FOX_LIBDIR && if [[ $$FOX_LIBDIR ]] ; then \
	echo "FOX_LIBDIR=$$FOX_LIBDIR" >> ${BUILDDIR}/Makefile.inc ; echo ; \
	echo "Please enter directory where FoX include files are kept:" ; \
	read FOX_INCDIR && echo "FOX_INCDIR=$$FOX_INCDIR" >> ${BUILDDIR}/Makefile.inc ; \
	echo "HAVE_FOX=1" >> ${BUILDDIR}/Makefile.inc ; \
	else echo "HAVE_FOX=0" >> ${BUILDDIR}/Makefile.inc ; fi
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
ifndef ARCH_LINKOPTS
	@echo
	@echo "Please enter any other linking options:"
	@echo "   Default: none"
	@read ARCH_LINKOPTS && echo "ARCH_LINKOPTS=$$ARCH_LINKOPTS" >> ${BUILDDIR}/Makefile.inc
endif
	echo "HAVE_LOTF=0" >> ${BUILDDIR}/Makefile.inc
	echo "HAVE_CP2K=0" >> ${BUILDDIR}/Makefile.inc


clean:
	for mods in  ${MODULES}  ; do \
	ln -sf ${PWD}/$$mods/Makefile ${BUILDDIR}/Makefile ; \
	${MAKE} -C ${BUILDDIR} -I${PWD}/Makefiles clean ; \
	done


