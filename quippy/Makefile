# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Portions of this code were written by
# HQ X     Tamas K. Stenczel, James Kermode
# HQ X
# HQ X   Copyright 2019
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   https://warwick.ac.uk/fac/sci/eng/staff/jrk
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# include other makefiles and export env variables
export QUIP_ARCH

ifeq (${QUIP_ROOT},)
QUIP_ROOT=${PWD}/..
endif
include Makefile.${QUIP_ARCH}
include ${BUILDDIR}/Makefile.inc
include Makefile.rules

PY_VERSION=$(shell ${PYTHON} -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)")
HAVE_PYTHON3=$(shell ${PYTHON} -c 'import sys; print(int(sys.version_info[0] >= 3))')

QUIPPY_SRC_DIR := ../../quippy

FPP = ${F90} -E -cpp ${QUIPPY_F90FLAGS}
CFLAGS = -fPIC ${QUIPPY_CFLAGS}


# The following files will be wrapped
LIBATOMS_SOURCES =  Atoms_types.f95 Atoms.f95 System.f95 Dictionary.f95 DynamicalSystem.f95 nye_tensor.f95
POT_SOURCES =  Potential.f95
GAP_SOURCES =  descriptors.f95

LIBATOMS_FILES = $(addprefix ../../src/libAtoms/,${LIBATOMS_SOURCES})
POT_FILES = $(addprefix ../../src/Potentials/,${POT_SOURCES})
GAP_FILES = $(addprefix ../../src/GAP/,${GAP_SOURCES})

WRAP_SOURCES = ${LIBATOMS_SOURCES} ${POT_SOURCES}
WRAP_FILES = ${LIBATOMS_FILES} ${POT_FILES}

PROGRAMS = quip
ifeq (${HAVE_GAP},1)
	WRAP_SOURCES +=  ${GAP_SOURCES}
	WRAP_FILES +=  ${GAP_FILES}
	PROGRAMS += gap_fit
endif

WRAP_FPP_FILES = $(addsuffix .fpp,$(basename ${WRAP_SOURCES}))

F90WRAP_FILES = $(addsuffix .f90,$(addprefix f90wrap_,$(basename ${WRAP_SOURCES})))

EXT_SUFFIX=$(shell ${PYTHON} -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
PLATFORM=$(shell ${PYTHON} -c 'import sysconfig; print(sysconfig.get_platform()+"-"+sysconfig.get_python_version())')

F90WRAP_OBJS= $(addsuffix .o,$(basename ${F90WRAP_FILES})) src.${PLATFORM}/_quippymodule.o src.${PLATFORM}/src.${PLATFORM}/fortranobject.o

F2PY_LINK_ARGS = $(shell ${PYTHON} -c 'import sys; print(" ".join([arg for arg in sys.argv if arg.startswith("-l") or arg.startswith("-L")] or arg == "-fopenmp"))' ${SYSLIBS} ${QUIPPY_LDFLAGS})


#################################
.PHONY: all build install clean f90wrap


all: build

f90wrap:
	echo Building f90wrap with CFLAGS=${CFLAGS}
	CC=${CC} CFLAGS=${CFLAGS} F90=${F90} ${PIP} install ${QUIP_ROOT}/src/f90wrap ${QUIPPY_INSTALL_OPTS}

clean:
	rm -f _quippy${EXT_SUFFIX} ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	rm -f f90wrap_*.o
	rm -rf src.*
	rm -rf build dist quippy quippy.egg-info

${WRAP_FPP_FILES}: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do \
	    echo ${FPP} ${DEFINES} $$f > $$(basename $${f%.*}).fpp ; \
	    ${FPP} ${DEFINES} $$f > $$(basename $${f%.*}).fpp ; \
	done

${F90WRAP_FILES}: ${WRAP_FPP_FILES}
	@if [ ${HAVE_PYTHON3} != 1 ] ; then echo "make quippy: Error: python3 needed. current version: ${PY_VERSION}."; exit 1; fi
	if [ -d "quippy" ] ; then rm -r quippy; fi
	mkdir quippy
	f90wrap -m quippy --f90-mod-name 'quippy._quippy' ${WRAP_FPP_FILES} -P -a c_error_abort \
		--init-file ${QUIPPY_SRC_DIR}/init.py  --move-methods \
		--shorten-routine-names -k ${QUIPPY_SRC_DIR}/KIND_MAP -s ${QUIPPY_SRC_DIR}/STRING_LENGTHS -S 10240 \
		--skip atoms_shallowcopy atoms_initialise_ptr atoms_finalise_multi \
		potential_initialise_inoutput cplx_2d  cplx_2d_array1_finalise \
                potential_local_e_mix_initialise potential_local_e_mix_finalise \
                --skip-types spherical_harmonics_type potential_local_e_mix --force-public set_cutoff \
		--documentation-plugin ${QUIPPY_SRC_DIR}/doc_plugin.py \
		--py-max-line-length 120 --f90-max-line-length 120

_quippy${EXT_SUFFIX} ${F90WRAP_OBJS}: ${F90WRAP_FILES}
	F90=${F90} LDFLAGS="${QUIPPY_LDFLAGS}" f2py-f90wrap --build-dir . -c -m _quippy ${F90WRAP_FILES} \
		-L. -lquip_nostub ${F2PY_LINK_ARGS} --f90flags="${QUIPPY_F90FLAGS}"

build: f90wrap _quippy${EXT_SUFFIX} ${QUIPPY_SRC_DIR}/quippy/*.py
	${QUIP_ROOT}/bin/gitversion --hash-only > VERSION
	cp -r ${QUIPPY_SRC_DIR}/quippy/*.py quippy

install: build
	cp ../../README.md ../../quippy/setup.py .
	cp ${PROGRAMS} quippy/
	${PIP} install ${QUIPPY_INSTALL_OPTS} .
