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
  NRL_TB_to_xml DFTB_to_xml FC_to_xml \
  calc_n_poles \
  analyze_md_phonons \
  xyz_to_SPECSYM_start xyz_to_SPECSYM_forces \
  local_E_fd wrap align local_random_search\
  metapot_test tabletest lotf_metapot lotfpc_metapot ts_calculation\
  lotf_time lotf_time_metapot cxyzreader fix_traj_latest

LIBS =  -L. -lquiputils -lquip -latoms 
LIBFILES = libatoms.a libquip.a libquiputils.a

.PHONY: clean allclean test

EXTRA_PROGRAMS = test_CASTEP_MM_buffer_crack test_CASTEP_water_chain test_CASTEP_water_bulk test ts_calculation QMMM_CP2K_md QMMM_md_buf annealing deval_dphonon

${PROGRAMS} ${EXTRA_PROGRAMS}: % :  ${OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${OBJS} ${LIBS}  ${LINKOPTS}

CRACK_OBJS = crackparams.o cracktools.o

crack: % :  ${CRACK_OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${CRACK_OBJS} ${LIBS}  ${LINKOPTS}

makecrack: % :  ${CRACK_OBJS} %.o ${LIBFILES}
	$(LINKER) $(LINKFLAGS) -o $@ ${F90OPTS} $@.o ${CRACK_OBJS} ${LIBS}  ${LINKOPTS}

clean:	
	-rm -f *.mod *.o ${PROGRAMS} ${EXTRA_PROGRAMS} crack makecrack

allclean : clean

