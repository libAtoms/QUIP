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

default: libquip_core.a 

ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.inc
include Makefile.rules

ALL_TARGETS = libquip_core.a

all: ${ALL_TARGETS}

BASE_F95_FILES = ScaLAPACK Matrix RS_SparseMatrix QUIP_Common Functions Ewald Multipole_Interactions  Partridge_Schwenke_Dipole
BASE_F95_SOURCES = ${addsuffix .f95, ${BASE_F95_FILES}}
BASE_F95_OBJS = ${addsuffix .o, ${BASE_F95_FILES}}

TB_F77_FILES = ginted
TB_F95_FILES = TB_Common  TB_Kpoints TBModel_NRL_TB_defs TBModel_NRL_TB TBModel_Bowler TBModel_DFTB TBModel_GSP \
  TBModel TBMatrix TB_Mixing  TBSystem  ApproxFermi TB_GreensFunctions TB 
TB_F77_SOURCES = ${addsuffix .f, ${TB_F77_FILES}}
TB_F95_SOURCES = ${addsuffix .f95, ${TB_F95_FILES}}
TB_F95_OBJS = ${addsuffix .o, ${TB_F77_FILES} ${TB_F95_FILES}}

IP_F95_FILES = IPEwald Yukawa Multipoles IPModel_GAP IPModel_LJ IPModel_Morse IPModel_FC IPModel_SW IPModel_Tersoff IPModel_EAM_Ercolessi_Adams IPModel_Brenner IPModel_FS IPModel_BOP IPModel_FB IPModel_Si_MEAM IPModel_Brenner_Screened IPModel_Brenner_2002 IPModel_TS IPModel_Glue IPModel_PartridgeSchwenke IPModel_Einstein IPModel_Coulomb IPModel_Sutton_Chen IPModel_FX IPModel_HFdimer IPModel_BornMayer IPModel_Custom IPModel_Template IPModel_SW_VP IPModel_WaterDimer_Gillan IPModel_WaterTrimer_Gillan IPModel_Tether IPModel_LMTO_TBE IPModel_Multipoles IPModel_FC4 IPModel_Spring IPModel_SCME IPModel_MTP IPModel_ZBL IPModel_LinearSOAP

ifeq (${HAVE_ASAP},1)
  IP_F95_FILES += IPModel_ASAP
endif
ifeq (${HAVE_KIM},1)
  IP_F95_FILES += IPModel_KIM
endif
IP_F95_FILES += IP
IP_F95_SOURCES = ${addsuffix .f95, ${IP_F95_FILES}}
IP_F95_OBJS = ${addsuffix .o, ${IP_F95_FILES}}

POT_F95_FILES = FilePot CallbackPot SocketPot Potential_simple AdjustablePotential Potential ElectrostaticEmbed quip_unified_wrapper quip_lammps_wrapper
ifeq (${HAVE_PRECON},1)
  POT_F95_FILES += Potential_Precon_Minim
endif
POT_F95_SOURCES = ${addsuffix .f95, ${POT_F95_FILES}}
POT_F95_OBJS = ${addsuffix .o, ${POT_F95_FILES}}

ALL_F95_FILES = ${BASE_F95_SOURCES}  ${IP_F95_SOURCES}
QUIP_OBJS = ${BASE_F95_OBJS}  ${IP_F95_OBJS} 
ifeq (${HAVE_TB},1)
  ALL_F95_FILES += ${TB_F95_SOURCES}
  QUIP_OBJS += ${TB_F95_OBJS}	
endif
QUIP_OBJS += ${POT_F95_OBJS}
ALL_F95_FILES += ${POT_F95_SOURCES}
ifeq (${HAVE_QC},1)
  QUIP_OBJS += QC_QUIP_Wrapper.o 
  ALL_F95_FILES += QC_QUIP_Wrapper.f95
endif


TEST_TARGETS = test_KP test_Matrix test_TBMatrix test_TBModel test_TBSystem test_TBCalculate  test_RS_SparseMatrix test_DFTB test_Potential test_Self_Consistency test_GreensFunc test_AF   test_ScaLAPACK test_SelfEnergy_passivation test_cat test_QC_QUIP_Wrapper

TARGETS = ${TEST_TARGETS}


.PHONY : clean allclean depend doc

libquip_core.a: ${QUIP_OBJS}
ifneq (${LIBTOOL},)
	${LIBTOOL} -o libquip_core.a ${QUIP_OBJS}
else
	${AR} ${AR_ADD} libquip_core.a $?
endif

Potential.o: \
	Potential_Hybrid_utils.f95 \
	Potential_Local_E_Mix_header.f95 Potential_Local_E_Mix_routines.f95 \
	Potential_ONIOM_header.f95 Potential_ONIOM_routines.f95 \
	Potential_ForceMixing_header.f95 Potential_ForceMixing_routines.f95

distribution: libquip_core.a
	mkdir -p distribution/lib distribution/include
	rm -f distribution/lib/* distribution/include/*
	cp extern${EXTERN_SUFFIX}/FoX/lib/* libAtoms/libatoms.a libquip_core.a distribution/lib
	cp extern${EXTERN_SUFFIX}/FoX/include/* libAtoms/*.mod *.mod distribution/include


${TARGETS}: % :  libquip_core.a %.o 
	${F95} ${LINKFLAGS} -o $@ $@.o libquip_core.a libatoms.a ${LINKOPTS}

clean:
	rm -f *.o *.mod *.mod.save *.mod.txt *.mod.save.txt *.fpp *.f90doc libquip_core.a quip.dtd  Potentials.depend

LIBS = -L. -lquip_core ${GP_LIB} -latoms
ifeq (${HAVE_THIRDPARTY},1)
  LIBS += -lthirdparty
endif

ifeq (${USE_MAKEDEP},1)
ALL_F95_FPP_FILES = ${addsuffix .fpp, ${basename ${ALL_F95_FILES}}}

Potentials.depend: ${ALL_F95_FPP_FILES}
	${SCRIPT_PATH}/${MAKEDEP} ${MAKEDEP_ARGS} -- ${addprefix ../../src/Potentials/,${ALL_F95_FILES}} > Potentials.depend

-include Potentials.depend
endif

install: 
