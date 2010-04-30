# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   libAtoms+QUIP: atomistic simulation library
# HND X
# HND X   Portions of this code were written by
# HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HND X
# HND X   Copyright 2006-2010.
# HND X
# HND X   Not for distribution
# HND X
# HND X   Portions of this code were written by Noam Bernstein as part of
# HND X   his employment for the U.S. Government, and are not subject
# HND X   to copyright in the USA.
# HND X
# HND X   When using this software, please cite the following reference:
# HND X
# HND X   http://www.libatoms.org
# HND X
# HND X  Additional contributions by
# HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HND X
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

all: ${PROGRAMS}

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

OBJS = 

PROGRAMS = teach_sparse

LIBS = -L. -lquiputils -lquip -lgp -latoms 
LIBFILES = libatoms.a libgp.a libquip.a libquiputils.a

.PHONY: clean allclean

${PROGRAMS}: %: ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS} ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean: clean
