all: ${PROGRAMS}

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

OBJS = 

PROGRAMS = elastic_fields density_1d decimate move_displacement_field density_new \
  clean_traj file_rewrite find_space_minim \
  structure_analysis_traj mean_var_correl

LIBS =  -L. -lquiputils -lquip -latoms 
LIBFILES = libatoms.a libquip.a libquiputils.a

.PHONY: clean allclean test


${PROGRAMS}: % :  ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

