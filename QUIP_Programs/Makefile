all: ${PROGRAMS}

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

OBJS = vacancy_map_mod.o 

PROGRAMS = eval md lj_example md_example \
  NRL_TB_to_xml DFTB_to_xml \
  calc_n_poles \
  analyze_md_phonons \
  xyz_to_SPECSYM_start xyz_to_SPECSYM_forces \
  local_E_fd wrap align local_random_search\
  metapot_test tabletest lotf_metapot lotfpc_metapot ts_calculation\
  lotf_time

LIBS =  -L. -lquiputils -lquip -latoms 
LIBFILES = libatoms.a libquip.a libquiputils.a

.PHONY: clean allclean test


${PROGRAMS}: % :  ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

test: % :  ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

CRACK_OBJS = crackparams.o cracktools.o

crack: % :  ${CRACK_OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${CRACK_OBJS} ${LIBS}  ${LINKOPTS}

makecrack: % :  ${CRACK_OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${CRACK_OBJS} ${LIBS}  ${LINKOPTS}

ts_calculation: % :  %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${LIBS}  ${LINKOPTS}

QMMM_CP2K_md: % :  %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${LIBS}  ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

