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

OBJS = vacancy_map_mod.o 

COBJS = vacancy_map_mod.o 

PROGRAMS = quip md \
  FC_to_xml analyze_md_phonons \
  xyz_to_SPECSYM_start xyz_to_SPECSYM_forces \
  local_E_fd wrap align local_random_search \
  metapot_test tabletest ts_calculation \
  fix_traj_latest randomise_calc \
  quip_wrapper_example descriptors_wrapper_example \
  slice_sample order_atoms_as_molecules md_gid \
  quip_wrapper_simple_example  \
  get_qw test_task_manager

CPROGRAMS = quip_wrapper_simple_example_C

ifeq (${HAVE_TB},1)
   PROGRAMS += NRL_TB_to_xml DFTB_to_xml calc_n_poles 
endif

default: all

all: ${PROGRAMS}

LIBS =  -L. -lquiputils -lquip_core ${GAP_LIB} -latoms
ifeq (${HAVE_THIRDPARTY},1)
  LIBS += -lthirdparty
endif
LIBFILES = libatoms.a ${GAP_LIBFILE} libquip_core.a libquiputils.a

.PHONY: clean allclean test

MINIM_PROGRAMS = minimtest1 dimertest1

EXTRA_PROGRAMS = test_CASTEP_MM_buffer_crack test_CASTEP_water_chain test_CASTEP_water_bulk test QMMM_CP2K_md QMMM_md_buf annealing deval_dphonon crack analytical_free_E_UI DTDE_GP_free_E DTDE_grad_GP_free_E ABF_GP_free_E fgp
-include Makefile.extra


${CPROGRAMS} : % :  ${COBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${COBJS} ${LIBS}  ${LINKOPTS}


${PROGRAMS} ${EXTRA_PROGRAMS} ${MINIM_PROGRAMS}: % :  ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS} ${EXTRA_PROGRAMS} ${MINIM_PROGRAMS} crack makecrack

allclean : clean

install:
	@if [ ! -d ${QUIP_INSTALLDIR} ]; then \
	  echo "make install: QUIP_INSTALLDIR '${QUIP_INSTALLDIR}' doesn't exist or isn't a directory"; \
	  exit 1; \
	else	\
	  for f in ${PROGRAMS} ; do \
	    echo "Copying $$f to ${QUIP_INSTALLDIR}/$${f}${QUIP_MPI_SUFFIX}" ; \
	    cp $$f ${QUIP_INSTALLDIR}/$${f}${QUIP_MPI_SUFFIX} ; \
	  done ;\
	  for f in ${CPROGRAMS} ; do \
	    echo "Copying $$f to ${QUIP_INSTALLDIR}/$${f}${QUIP_MPI_SUFFIX}" ; \
	    cp $$f ${QUIP_INSTALLDIR}/$${f}${QUIP_MPI_SUFFIX} ; \
	  done ;\
	fi
