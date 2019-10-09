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

QUIPPY_SRC_DIR := ../../quippy

FPP = ${F90} -E -cpp ${QUIPPY_F90FLAGS}
CFLAGS = -fPIC ${QUIPPY_CFLAGS}


# The following files will be wrapped
LIBATOMS_SOURCES =  Atoms_types.f95 Atoms.f95 System.f95 Dictionary.f95 DynamicalSystem.f95
POT_SOURCES =  Potential.f95
GAP_SOURCES =  descriptors.f95 gap_fit_module.f95

LIBATOMS_FILES = $(addprefix ../../src/libAtoms/,${LIBATOMS_SOURCES})
POT_FILES = $(addprefix ../../src/Potentials/,${POT_SOURCES})
GAP_FILES = $(addprefix ../../src/GAP/,${GAP_SOURCES})

WRAP_SOURCES = ${LIBATOMS_SOURCES} ${POT_SOURCES}
WRAP_FILES = ${LIBATOMS_FILES} ${POT_FILES}

ifeq (${HAVE_GAP},1)
	WRAP_SOURCES +=  ${GAP_SOURCES}
	WRAP_FILES +=  ${GAP_FILES}
endif

WRAP_FPP_FILES = $(addsuffix .fpp,$(basename ${WRAP_SOURCES}))

F90WRAP_FILES = $(addsuffix .f90,$(addprefix f90wrap_,$(basename ${WRAP_SOURCES})))

F2PY_LINK_ARGS = $(shell python -c 'import sys; print(" ".join([arg for arg in sys.argv if arg.startswith("-l") or arg.startswith("-L")]))' ${SYSLIBS} ${QUIPPY_LDFLAGS})

F90WRAP_PATH = $(shell pip show f90wrap | grep Location | sed -e 's/Location: //' -e 's/lib.*/bin/')
F90WRAP = ${F90WRAP_PATH}/f90wrap
F2PY_F90WRAP = ${F90WRAP_PATH}/f2py-f90wrap

#################################
.PHONY: all build install clean


all: build

clean:
	rm -f _quippy.so _quippy.*.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	rm -f f90wrap_*.o
	rm -rf src.*
	rm -rf build dist quippy quippy.egg-info

${WRAP_FPP_FILES}: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do \
	    echo ${FPP} -D'SIZEOF_FORTRAN_T=2' $$f > $$(basename $${f%.*}).fpp ; \
	    ${FPP} -D'SIZEOF_FORTRAN_T=2' $$f > $$(basename $${f%.*}).fpp ; \
	done

_quippy.so: ${WRAP_FPP_FILES}
	if [ -d "quippy" ] ; then rm -r quippy; fi
	mkdir quippy
	${F90WRAP} -m quippy ${WRAP_FPP_FILES} -P -a c_error_abort \
		--init-file ${QUIPPY_SRC_DIR}/init.py  --move-methods \
		--shorten-routine-names -k ${QUIPPY_SRC_DIR}/KIND_MAP -s ${QUIPPY_SRC_DIR}/STRING_LENGTHS -S 10240 \
		--skip atoms_shallowcopy atoms_initialise_ptr atoms_finalise_multi \
		potential_initialise_inoutput cplx_2d  cplx_2d_array1_finalise \
		--skip-types spherical_harmonics_type --force-public set_cutoff \
		--documentation-plugin ${QUIPPY_SRC_DIR}/doc_plugin.py \
		--py-max-line-length 120 --f90-max-line-length 120

	F90=${F90} ${F2PY_F90WRAP} --build-dir . -c -m _quippy ${F90WRAP_FILES} ${F2PY_LINK_ARGS} \
		-L. -lquip --link-lapack --f90flags="${QUIPPY_F90FLAGS}" -lgomp
	if [ -d "_quippy.so" ] ; then rm _quippy.so; fi
	ln -s _quippy.*.so _quippy.so

build: _quippy.so ${QUIPPY_SRC_DIR}/quippy/*.py
	cp -r ${QUIPPY_SRC_DIR}/quippy/*.py quippy

install: build
	python ../../quippy/setup.py install ${QUIPPY_INSTALL_OPTS}
