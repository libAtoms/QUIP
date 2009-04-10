include Makefile.common
ifeq (${MAKEFILE_ARCH_SUFFIX},)
  include Makefile.arch
else
  include Makefile.${MAKEFILE_ARCH_SUFFIX}
endif

include Makefile.programs

ifneq (${BUILDDIR},)
	export BUILDDIR
	export BUILDPRE = ../
	export VPATH += ..
endif

ifeq (${HAVE_CP2K},1)
  COMMON_F95FLAGS += -DHAVE_CP2K
endif
ifeq (${HAVE_LOTF},1)
  LOTFLIB_DIR = ../${BUILDPRE}/LOTF95
  LOTFLIB = -llotf
  COMMON_F95FLAGS += -DHAVE_LOTF
endif
ifeq (${HAVE_QUIP},1)
  QUIPLIB_DIR = ../${BUILDPRE}/QUIP ../${BUILDPRE}/QUIP_Utils
  QUIPLIB = -lquiputils -lquip
endif
ifeq (${HAVE_HYBRID},1)
  COMMON_F95FLAGS += -DHAVE_HYBRID
endif
ifeq (${HAVE_NETCDF},1)
  COMMON_F95FLAGS += -DHAVE_NETCDF ${ARCH_NETCDF_INCLUDEDIR}
  ARCH_LIBS += ${ARCH_NETCDF_LIBDIR} ${ARCH_NETCDF_LIBS}
endif

ifeq (${HAVE_MDCORE},1)
  COMMON_F95FLAGS += ${ARCH_MDCORE_INCLUDEDIRS}
  ARCH_LIBS += ${ARCH_MDCORE_LIBDIRS} ${ARCH_MDCORE_LIBS}
endif

DEPLIBDIRS = ../${BUILDPRE}/libAtoms ${QUIPLIB_DIR} ${LOTFLIB_DIR}
DEPLIBS =  ${QUIPLIB} $(LOTFLIB) -latoms 

COMMON_F95FLAGS += $(addprefix -I,$(addsuffix /${BUILDDIR}, ${DEPLIBDIRS}))
COMMON_F95FLAGS += $(addprefix -I,${DEPLIBDIRS})
COMMON_F95FLAGS += ${ARCH_FOX_INCLUDEDIR}
ARCH_LIBS += ${ARCH_FOX_LIBDIR} ${ARCH_FOX_LIBS}

LIBDIRS += $(addprefix -L,$(addsuffix /${BUILDDIR}, ${DEPLIBDIRS}))

.PHONY: clean allclean ${DEPLIBDIRS} test

all: ${PROGRAMS}

${DEPLIBDIRS}:
	-@if [ -d $@ ] ; then mkdir -p $@/${BUILDDIR} && cd $@/${BUILDDIR} && ln -s ../Makefile . ; else echo ; echo ; echo "Directory $@ not found, skipping" ; echo ; echo ; fi
	${MAKE} -C $@/${BUILDDIR} 

${PROGRAMS}: % : ${DEPLIBDIRS} ${OBJS} %.o
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBDIRS} ${DEPLIBS}  ${ARCH_LAPACK_LIBS} ${ARCH_LIBS}

${TEST_PROGRAMS}: % : ${DEPLIBDIRS} ${OBJS} %.o
	mkdir -p tests
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBDIRS} ${DEPLIBS} ${ARCH_LAPACK_LIBS} ${ARCH_LIBS}


clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

deepclean: clean
	for f in ${DEPLIBDIRS}; do if [ -d $$f/${BUILDDIR} ] ; then ${MAKE} -C $$f/${BUILDDIR} clean; fi ; done

test: ${TEST_PROGRAMS}
	cd ${BUILDPRE}; ./test.py ${MAKEFILE_ARCH_SUFFIX} -r tests/*.test; cd -
