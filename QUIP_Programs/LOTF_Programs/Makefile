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

ifeq (${HAVE_NETCDF},1)
  COMMON_F95FLAGS += -DHAVE_NETCDF ${ARCH_NETCDF_INCLUDEDIR} ${ARCH_NETCDF4}
  ARCH_LIBS += ${ARCH_NETCDF_LIBDIR} ${ARCH_NETCDF_LIBS}
endif

DEPLIBDIRS = ../${BUILDPRE}/libAtoms ../${BUILDPRE}/LOTF95 
DEPLIBS = -llotf -latoms


ifeq (${HAVE_QUIP},1)
  DEPLIBDIRS += ../${BUILDPRE}/QUIP ../${BUILDPRE}/QUIP_Utils
  DEPLIBS = -llotf -lquip -llotf -lquiputils -latoms
  ARCH_LIBS += ${ARCH_FOX_LIBDIR} ${ARCH_FOX_LIBS}
endif

ifeq (${HAVE_MDCORE},1)
  COMMON_F95FLAGS += ${ARCH_MDCORE_INCLUDEDIR}
  ARCH_LIBS += ${ARCH_MDCORE_LIBDIR} ${ARCH_MDCORE_LIBS}
endif

COMMON_F95FLAGS += $(addprefix -I,$(addsuffix /${BUILDDIR}, ${DEPLIBDIRS}))
LIBDIRS += $(addprefix -L,$(addsuffix /${BUILDDIR}, ${DEPLIBDIRS}))

MODULES = crackparams cracktools

MOD_F95FILES = ${addsuffix .f95, ${MODULES}}
MOD_OBJS = ${addsuffix .o, ${MODULES}}


.PHONY: clean allclean ${DEPLIBDIRS}


all: ${PROGRAMS}

${DEPLIBDIRS}:
	-@if [ -d $@ ] ; then mkdir -p $@/${BUILDDIR} && ln -s $@/Makefile $@/${BUILDDIR}/Makefile ; else echo ; echo ; echo "Directory $@ not found, skipping" ; echo ; echo ; fi
	${MAKE} -C $@/${BUILDDIR} 

${PROGRAMS}: % : ${DEPLIBDIRS} ${MOD_OBJS} %.o 
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${MOD_OBJS} ${LIBDIRS}  ${DEPLIBS} ${ARCH_LIBS} ${ARCH_LAPACK_LIBS}

${TEST_PROGRAMS}: % : ${DEPLIBDIRS} ${OBJS} %.o
	mkdir -p tests
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBDIRS} ${DEPLIBS} ${ARCH_LIBS} ${ARCH_LAPACK_LIBS} 


clean:	
	-rm -f *.mod ${PROGRAMS} $(addsuffix .o,${PROGRAMS}) ${MOD_OBJS}

allclean : clean

deepclean: clean
	for f in ${DEPLIBDIRS}; do if [ -d $$f/${BUILDDIR} ] ; then ${MAKE} -C $$f/${BUILDDIR} clean; fi ; done

test: ${TEST_PROGRAMS}
	cd ${BUILDPRE}; ../QUIP_Programs/test.py ${MAKEFILE_ARCH_SUFFIX} -r tests/*.test; cd -
