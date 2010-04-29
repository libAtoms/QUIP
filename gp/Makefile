all: libgp.a

include Makefile.inc
include Makefile.common

ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

GP_F95_FILES = bispectrum clustering gp gp_sparse qw_continous
GP_F95_SOURCES = ${addsuffix .f95, ${GP_F95_FILES}}
GP_F95_OBJS = ${addsuffix .o, ${GP_F95_FILES}}

.PHONY : clean allclean depend doc

libgp.a: ${GP_F95_OBJS}
	ar ${AR_ADD} libgp.a $?

clean:
	rm -f *.o *.mod *.mod.save libquip.a
	-rm -r doctmp

allclean: clean
	rm -f depend 

-include depend
