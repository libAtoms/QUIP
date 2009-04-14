include Makefile.inc
include Makefile.common
ifeq (${ARCH},)
  include Makefile.arch
else
  include Makefile.${ARCH}
endif

#include Makefile.programs
PROGRAMS = eval

LIBS =  -L. -lquiputils -lquip -latoms 


.PHONY: clean allclean test

all: ${PROGRAMS}

${PROGRAMS}: % :  ${OBJS} %.o
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

