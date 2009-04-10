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
		export VPATH+=..
	endif
endif

DEPLIBDIRS = ../${BUILDPRE}/libAtoms ../${BUILDPRE}/QUIP 

ifeq (${HAVE_LOTF},1)
	DEPLIBDIRS += ../${BUILDPRE}/LOTF95
endif

COMMON_F95FLAGS += $(addsuffix /${BUILDDIR}, $(addprefix -I,${DEPLIBDIRS}))
COMMON_F95FLAGS += ${ARCH_FOX_INCLUDEDIR} ${ARCH_MDCORE_INCLUDEDIR}

FOX_LIBS = -lFoX_sax -lFoX_wxml -lFoX_common -lfsys
COMMON_LIBS += -latoms $(FOX_LIBS)

COMMON_LIBDIRS = ${ARCH_FOX_LIBDIR} ${LOTF95LIBDIR} ${DEPLIBDIRS}
LIBS = ${FOX_LIBS} ${COMMON_LIBS} ${ARCH_LIBS}


F95_FILES =  transition_state elasticity phonons potential_test
OBJS = ${addsuffix .o, ${F95_FILES}}

ifeq (${HAVE_NETCDF},1)
	COMMON_F95FLAGS += -DHAVE_NETCDF ${ARCH_NETCDF_INCLUDEDIR}
endif

all: libquiputils.a

.PHONY : clean allclean  ${DEPLIBDIRS}

${DEPLIBDIRS}:
	-mkdir -p $@/${BUILDDIR} && ln -s $@/Makefile $@/${BUILDDIR}/Makefile
	${MAKE} -C $@/${BUILDDIR} 

libquiputils.a: ${DEPLIBDIRS} ${OBJS}
	ar ${AR_ADD} libquiputils.a ${OBJS}

clean:
	rm -f *.o *.mod *.mod.save libquip.a

allclean: clean

deepclean: allclean
	for f in ${DEPLIBDIRS} ; do if [ -d $$f/${BUILDDIR} ] ; then ${MAKE} -C $$f/${BUILDDIR} clean; fi ; done


