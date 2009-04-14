include Makefile.inc
include Makefile.common
ifeq (${ARCH},)
  include Makefile.arch
else
  include Makefile.${ARCH}
endif


F95_FILES =  transition_state elasticity phonons potential_test
OBJS = ${addsuffix .o, ${F95_FILES}}

all: libquiputils.a

.PHONY : clean allclean 


libquiputils.a:  ${OBJS}
	ar ${AR_ADD} libquiputils.a ${OBJS}

clean:
	rm -f *.o *.mod *.mod.save libquiputils.a

allclean: clean
