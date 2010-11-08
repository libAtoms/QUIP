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

all: libgap_teach.a

include Makefile.features
include Makefile.libs
include Makefile.rules
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.override_arch

GP_TEACH_F95_FILES = clustering gp_teach 
GP_TEACH_F95_SOURCES = ${addsuffix .f95, ${GP_TEACH_F95_FILES}}
GP_TEACH_F95_OBJS = ${addsuffix .o, ${GP_TEACH_F95_FILES}}

.PHONY : clean allclean depend doc

libgap_teach.a: ${GP_TEACH_F95_OBJS}
	${AR} ${AR_ADD} libgap_teach.a $?

clean:
	rm -f *.o *.mod *.mod.save libgap_teach.a
	-rm -r doctmp

allclean: clean
	rm -f depend 

-include depend
