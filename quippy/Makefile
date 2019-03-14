# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

####### NOTES
# original Makefile taken from branch hackathon-f90wrap




# default to be overridden by site-specific config
QUIPPY_INSTALL_DIR = ${HOME}/anaconda2/lib/python2.7/site-packages     #fixme: this was for debug only

# include other makefiles and export env variables
export QUIP_ARCH

ifeq (${QUIP_ROOT},)
QUIP_ROOT=${PWD}/..
endif
include Makefile.${QUIP_ARCH}
include ${BUILDDIR}/Makefile.inc
include Makefile.rules

# working on python3
PYTHON := python

QUIPPY_SRC_DIR := ../../quippy

FPP = gfortran -E -cpp
CFLAGS = -fPIC


LIBATOMS_SOURCES =  Atoms_types.f95 Atoms.f95 System.f95 Dictionary.f95   # Table.f95
POT_SOURCES =  Potential.f95
GAP_SOURCES =  descriptors.f95 descriptors_wrapper.f95

LIBATOMS_FILES = $(addprefix ../../src/libAtoms/,${LIBATOMS_SOURCES})
POT_FILES = $(addprefix ../../src/Potentials/,${POT_SOURCES})
GAP_FILES = $(addprefix ../../src/GAP/,${GAP_SOURCES})

ifeq (${HAVE_GAP},1)
	WRAP_SOURCES = ${LIBATOMS_SOURCES} ${POT_SOURCES} ${GAP_SOURCES}
	WRAP_FILES = ${LIBATOMS_FILES} ${POT_FILES} ${GAP_FILES}
else
	WRAP_SOURCES = ${LIBATOMS_SOURCES} ${POT_SOURCES}
	WRAP_FILES = ${LIBATOMS_FILES} ${POT_FILES}
endif

WRAP_FPP_FILES = $(addsuffix .fpp,$(basename ${WRAP_SOURCES}))

F90WRAP_FILES = $(addsuffix .f90,$(addprefix f90wrap_,$(basename ${WRAP_SOURCES}))) #f90wrap_toplevel.f90

F2PY_LINK_ARGS = $(shell python -c 'import sys; print " ".join([arg for arg in sys.argv if arg.startswith("-l") or arg.startswith("-L")])' ${SYSLIBS} ${QUIPPY_LDFLAGS})


###
.PHONY: build clean install

build: _quippy.so
	echo WRAP_SOURCES__ ${WRAP_SOURCES}
	echo WRAP_FILES___ ${WRAP_FILES}

clean:
	-rm _castep.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -r castep
	-rm -r f90wrap_*.o
	-rm -rf src.*

${WRAP_FPP_FILES}: ${WRAP_FILES}
	echo quippy_stone_1 ### debug
	for f in ${WRAP_FILES}; do \
	    echo ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
	    ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
	done

_quippy.so: ${WRAP_FPP_FILES}
	echo ${WRAP_FILES}
	echo ${WRAP_SOURCES}
	echo quippy_stone_2 ### debug
	-rm -r quippy
	mkdir quippy
	echo quippy_stone_3 ### debug
	cp -r ${QUIPPY_SRC_DIR}/quippy/* quippy
	echo quippy_stone_4 ### debug
	f90wrap -v -m quippy ${WRAP_FPP_FILES} -P \
	   -a c_error_abort --init-file ${QUIPPY_SRC_DIR}/init.py  --move-methods \
		 --shorten-routine-names -k ${QUIPPY_SRC_DIR}/KIND_MAP -s ${QUIPPY_SRC_DIR}/STRING_LENGTHS -S 10240 \
		 --skip atoms_shallowcopy atoms_initialise_ptr atoms_finalise_multi \
		         potential_initialise_inoutput #potential_fm potential_evb potential_cluster
	echo quippy_stone_5 ### debug

	f2py-f90wrap --build-dir . -c -m _quippy ${F90WRAP_FILES} ${F2PY_LINK_ARGS} \
		 -L. -lquip --link-lapack

install: build
	QUIPPY_INSTALL_DIR=/home/tks32/.local/lib/python3.6/site-packages       #fixme: this was for debug only
	cp -r quippy ${QUIPPY_INSTALL_DIR}
	cp _quippy.so ${QUIPPY_INSTALL_DIR}
