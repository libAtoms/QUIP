all: libquip.a

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

BASE_F95_FILES = MPI_context ScaLAPACK Matrix RS_SparseMatrix QUIP_Common
BASE_F95_SOURCES = ${addsuffix .f95, ${BASE_F95_FILES}}
BASE_F95_OBJS = ${addsuffix .o, ${BASE_F95_FILES}}

TB_F77_FILES = ginted
TB_F95_FILES = TB_Common  TB_Kpoints TBModel_NRL_TB_defs TBModel_NRL_TB TBModel_Bowler TBModel_DFTB TBModel_GSP \
  TBModel TBMatrix TB_Mixing Functions Ewald TBSystem  ApproxFermi TB_GreensFunctions TB FilePot
TB_F77_SOURCES = ${addsuffix .f, ${TB_F77_FILES}}
TB_F95_SOURCES = ${addsuffix .f95, ${TB_F95_FILES}}
TB_F95_OBJS = ${addsuffix .o, ${TB_F77_FILES} ${TB_F95_FILES}}

IP_F95_FILES = IPModel_GAP IPModel_LJ IPModel_SW IPModel_Tersoff IPModel_EAM_Ercolessi_Adams IPModel_Brenner IPModel_FS IPModel_BOP IPEwald IPModel_FB IPModel_Si_MEAM IPModel_Brenner_Screened IPModel_Brenner_2002 IPModel_ASAP IPModel_ASAP2 IPModel_Template IP
IP_F95_SOURCES = ${addsuffix .f95, ${IP_F95_FILES}}
IP_F95_OBJS = ${addsuffix .o, ${IP_F95_FILES}}

POT_F95_FILES = Potential AdjustablePotential MetaPotential QUIP_module
POT_F95_SOURCES = ${addsuffix .f95, ${POT_F95_FILES}}
POT_F95_OBJS = ${addsuffix .o, ${POT_F95_FILES}}

ALL_F95_FILES = ${BASE_F95_SOURCES} ${TB_F95_SOURCES} ${IP_F95_SOURCES}
QUIP_OBJS = ${BASE_F95_OBJS} ${TB_F95_OBJS} ${IP_F95_OBJS} 
ifeq (${HAVE_CP2K},1)
  QUIP_OBJS += cp2k_driver_old.o cp2k_driver.o
  ALL_F95_FILES += cp2k_driver_old.f95 cp2k_driver.f95
endif
QUIP_OBJS += ${POT_F95_OBJS}
ALL_F95_FILES += ${POT_F95_SOURCES}
ifeq (${HAVE_QC},1)
  QUIP_OBJS += QC_QUIP_Wrapper.o 
  ALL_F95_FILES += QC_QUIP_Wrapper.f95
endif


TEST_TARGETS = test_KP test_Matrix test_TBMatrix test_TBModel test_TBSystem test_TBCalculate  test_RS_SparseMatrix test_DFTB test_Potential test_Self_Consistency test_GreensFunc test_AF   test_ScaLAPACK test_SelfEnergy_passivation test_cat test_QC_QUIP_Wrapper

# UTIL_TARGETS = NRL_TB_to_xml DFTB_to_xml
# 
# QUIP_TARGETS = eval

TARGETS = ${TEST_TARGETS}


.PHONY : clean allclean depend doc

libquip.a: ${QUIP_OBJS}
	ar ${AR_ADD} libquip.a ${QUIP_OBJS}

MetaPotential.o: \
	MetaPotential_Hybrid_utils.f95 \
	MetaPotential_Local_E_Mix_header.f95 MetaPotential_Local_E_Mix_routines.f95 \
	MetaPotential_ONIOM_header.f95 MetaPotential_ONIOM_routines.f95 \
	MetaPotential_ForceMixing_header.f95 MetaPotential_ForceMixing_routines.f95

distribution: libquip.a
	mkdir -p distribution/lib distribution/include
	rm -f distribution/lib/* distribution/include/*
	cp extern${EXTERN_SUFFIX}/FoX/lib/* libAtoms/libatoms.a libquip.a distribution/lib
	cp extern${EXTERN_SUFFIX}/FoX/include/* libAtoms/*.mod *.mod distribution/include


${TARGETS}: % :  libquip.a %.o 
	${F95} ${LINKFLAGS} -o $@ $@.o libquip.a libatoms.a ${LINKOPTS}

clean:
	rm -f *.o *.mod *.mod.save libquip.a
	-rm -r doctmp

allclean: clean
	rm -f depend 


quip.dtd: ${IP_F95_SOURCES} ${TB_F95_SOURCES} mkdtd.pl
	${VPATH}/mkdtd.pl ${IP_F95_SOURCES} ${TB_F95_SOURCES} > $@

LIBS = -L. -lquip -latoms

cp2k_filepot_old: % : %.o libatoms.a libquip.a
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${LIBS} ${LINKOPTS}

cp2k_filepot: % : %.o libatoms.a libquip.a
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${LIBS} ${LINKOPTS}

# depend:
# 	${MAKE} -f Makefile.depend

-include depend

