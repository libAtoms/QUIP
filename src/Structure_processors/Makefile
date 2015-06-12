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

include Makefile.${QUIP_ARCH}
include Makefile.inc
include Makefile.rules

OBJS = 

PROGRAMS = elastic_fields density_1d decimate move_displacement_field density_new \
  clean_traj file_rewrite find_space_minim \
  structure_analysis_traj mean_var_correl angle_distr coordination_number \
  make_bulk_supercell \
  extract_cv extract_EVB mean_var_decorrelated_err rings UI_integrate solvate \
  merge_traj test_ran diffusion xyz2residue_library make_surface_slab_vasp

ifeq (${HAVE_TB},1)
  PROGRAMS += make_k_mesh
endif

default: all

all: ${PROGRAMS}

LIBS =  -L. -lquiputils -lquip_core ${GP_LIB} -latoms 
ifeq (${HAVE_THIRDPARTY},1)
  LIBS += -lthirdparty
endif
LIBFILES = libatoms.a libquip_core.a ${GP_LIBFILE} libquiputils.a

.PHONY: clean allclean test


${PROGRAMS}: % :  ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS}

allclean : clean

install-scripts:
	@echo "installing from Structure_processors"; \
	for f in `/bin/ls src/Structure_processors | egrep -v '\.f95|Makefile*'`; do \
	  if [ -f src/Structure_processors/$$f ] && [ -x src/Structure_processors/$$f ]; then \
	    echo "copying f $$f to ${QUIP_INSTDIR}"; \
	    cp src/Structure_processors/$$f ${QUIP_INSTDIR}; \
	  fi; \
	done

