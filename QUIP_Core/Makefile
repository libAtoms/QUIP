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
	SRCDIR=..
else
	SRCDIR=.
endif

COMMON_F95FLAGS += -I$(SRCDIR)

DEPLIBDIRS = ../${BUILDPRE}/libAtoms
ifeq (${HAVE_LOTF},1)
	LOTF95LIBDIR = ../${BUILDPRE}/LOTF95
	COMMON_F95FLAGS += -DHAVE_LOTF -I${LOTF95LIBDIR}/${BUILDDIR}
endif

ifeq (${HAVE_CP2K},1)
	COMMON_F95FLAGS += -DHAVE_CP2K
endif

ifeq (${HAVE_NETCDF},1)
	COMMON_F95FLAGS += -DHAVE_NETCDF ${ARCH_NETCDF_INCLUDEDIR}
endif

ifeq (${HAVE_LOCAL_E_MIX},1)
	COMMON_F95FLAGS += -DHAVE_LOCAL_E_MIX
endif

ifeq (${HAVE_ONIOM},1)
	COMMON_F95FLAGS += -DHAVE_ONIOM
endif

ifeq (${HAVE_MDCORE},1)
	COMMON_F95FLAGS += ${ARCH_MDCORE_INCLUDEDIR} -DHAVE_MDCORE
endif

COMMON_F95FLAGS += $(addsuffix /${BUILDDIR}, $(addprefix -I,${DEPLIBDIRS}))
COMMON_F95FLAGS += ${ARCH_FOX_INCLUDEDIR}

FOX_LIBS = -lFoX_sax -lFoX_wxml -lFoX_common -lfsys
COMMON_LIBS += -latoms $(FOX_LIBS)
ifneq (${HAVE_LOTF},1)
	COMMON_LIBS += -llotf
endif

COMMON_LIBDIRS = ${ARCH_FOX_LIBDIR} ${LOTF95LIBDIR} ${DEPLIBDIRS}
LIBS = ${FOX_LIBS} ${COMMON_LIBS} ${ARCH_LIBS}


BASE_F95_FILES = MPI_context ScaLAPACK Matrix RS_SparseMatrix QUIP_Common
BASE_F95_SOURCES = ${addsuffix .f95, ${BASE_F95_FILES}}
BASE_F95_OBJS = ${addsuffix .o, ${BASE_F95_FILES}}

TB_F95_FILES = TB_Common  TB_Kpoints TBModel_NRL_TB TBModel_Bowler TBModel_DFTB TBModel_GSP \
  TBModel TBMatrix TB_Mixing Functions Ewald ginted TBSystem  ApproxFermi TB_GreensFunctions TB FilePot
TB_F95_SOURCES = ${addsuffix .f95, ${TB_F95_FILES}}
TB_F95_OBJS = ${addsuffix .o, ${TB_F95_FILES}}

IP_F95_FILES = IPModel_GAP IPModel_LJ IPModel_SW IPModel_Tersoff IPModel_EAM_Ercolessi_Adams IPModel_Brenner IPModel_FS IPModel_BOP IPEwald IPModel_FB IPModel_Si_MEAM IPModel_Brenner_Screened IPModel_Brenner_2002 IP
IP_F95_SOURCES = ${addsuffix .f95, ${IP_F95_FILES}}
IP_F95_OBJS = ${addsuffix .o, ${IP_F95_FILES}}

QUIP_OBJS_PASS1 = ${BASE_F95_OBJS} ${TB_F95_OBJS} ${IP_F95_OBJS} 
ifeq (${HAVE_CP2K},1)
  QUIP_OBJS_PASS1 += cp2k_driver.o 
endif
QUIP_OBJS_PASS1 += Potential.o
QUIP_OBJS_PASS2 = MetaPotential.o QUIP_module.o 
ifeq (${HAVE_QC},1)
  QUIP_OBJS_PASS2 += QC_QUIP_Wrapper.o 
endif

ALL_F95_FILES = ${BASE_F95_SOURCES} ${TB_F95_SOURCES} ${IP_F95_SOURCES}
ifeq (${HAVE_CP2K},1)
  ALL_F95_FILES += cp2k_driver.f95
endif
ALL_F95_FILES += Potential.f95 MetaPotential.f95 QUIP_module.f95 
ifeq (${HAVE_QC},1)
  ALL_F95_FILES += QC_QUIP_Wrapper.f95
endif

TEST_TARGETS = test_KP test_Matrix test_TBMatrix test_TBModel test_TBSystem test_TBCalculate \
  test_RS_SparseMatrix test_DFTB test_Potential test_Self_Consistency test_GreensFunc test_AF \
  test_ScaLAPACK test_SelfEnergy_passivation test_cat test_QC_QUIP_Wrapper

# UTIL_TARGETS = NRL_TB_to_xml DFTB_to_xml
# 
# QUIP_TARGETS = eval

TARGETS = ${TEST_TARGETS}

all: libquip.a

.PHONY : clean allclean depend doc ${DEPLIBDIRS} ${LOTF95LIBDIR}

${DEPLIBDIRS}:
	-mkdir -p $@/${BUILDDIR} && cd $@/${BUILDDIR} && ln -s ../Makefile .
	${MAKE} -C $@/${BUILDDIR} 

${LOTF95LIBDIR}:
	-mkdir -p $@/${BUILDDIR} && cd $@/${BUILDDIR} && ln -s ../Makefile .
	${MAKE} -C $@/${BUILDDIR} 

pass1: ${DEPLIBDIRS} ${QUIP_OBJS_PASS1}

pass2: pass1 ${QUIP_OBJS_PASS2}

libquip.a: pass1 ${LOTF95LIBDIR} pass2
	ar ${AR_ADD} libquip.a ${QUIP_OBJS_PASS1} ${QUIP_OBJS_PASS2}

MetaPotential.o: \
	MetaPotential_Hybrid_utils.f95 \
	MetaPotential_Local_E_Mix_header.f95 MetaPotential_Local_E_Mix_routines.f95 \
	MetaPotential_ONIOM_header.f95 MetaPotential_ONIOM_routines.f95 \
	MetaPotential_ForceMixing_header.f95 MetaPotential_ForceMixing_routines.f95

distribution: ${DEPLIBDIRS} libquip.a
	mkdir -p distribution/lib distribution/include
	rm -f distribution/lib/* distribution/include/*
	cp extern${EXTERN_SUFFIX}/FoX/lib/* libAtoms/libatoms.a libquip.a distribution/lib
	cp extern${EXTERN_SUFFIX}/FoX/include/* libAtoms/*.mod *.mod distribution/include


${TARGETS}: % : ${DEPLIBDIRS} ${LOTF95LIBDIR} libquip.a %.o 
	${F95} ${LINKFLAGS} -o $@ $@.o libquip.a ${LIBDIRS} ${LIBS}

clean:
	rm -f *.o *.mod *.mod.save libquip.a
	-rm -r doctmp

allclean: clean
	rm -f depend 

deepclean: allclean
	for f in ${DEPLIBDIRS} ${LOTF95LIBDIR} ; do if [ -d $$f/${BUILDDIR} ] ; then ${MAKE} -C $$f/${BUILDDIR} clean; fi ; done


doc:
	./mkdoc ${ALL_F95_FILES}

quip.dtd: ${IP_F95_SOURCES} ${TB_F95_SOURCES} mkdtd.pl
	./mkdtd.pl ${IP_F95_SOURCES} ${TB_F95_SOURCES} > $@

# depend:
# 	${MAKE} -f Makefile.depend

-include depend

