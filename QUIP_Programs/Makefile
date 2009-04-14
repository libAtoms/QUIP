include Makefile.inc
include Makefile.common
ifeq (${ARCH},)
  include Makefile.arch
else
  include Makefile.${ARCH}
endif

OBJS = vacancy_map_mod.o

PROGRAMS = eval md lj_example md_example \
  NRL_TB_to_xml DFTB_to_xml \
  calc_n_poles \
  analyze_md_phonons \
  xyz_to_SPECSYM_start xyz_to_SPECSYM_forces \
  local_E_fd wrap align test local_random_search


LIBS =  -L. -lquiputils -lquip -latoms 


.PHONY: clean allclean test

all: ${PROGRAMS}

${PROGRAMS}: % :  ${OBJS} %.o
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

CRACK_OBJS = crackparams.o cracktools.o

crack: % :  ${CRACK_OBJS} %.o
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${CRACK_OBJS} ${LIBS}  ${LINKOPTS}


clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

