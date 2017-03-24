# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2017
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

# default to be overridden by site-specific config
QUIPPY_INSTALL_DIR = ${HOME}/lib/python

# include other makefiles and export env variables
export QUIP_ARCH

ifeq (${QUIP_ROOT},)
QUIP_ROOT=${PWD}/..
endif
include Makefile.${QUIP_ARCH}
include ${BUILDDIR}/Makefile.inc
include Makefile.rules

# Test if the python in the path is a Python 2, if not then
# explicitly call python2
# This is a workaround for old virtualenv versions that dont create the
# python2 symlink.
python_version := $(word 2, $(subst ., , $(shell python -V 2>&1)))
ifeq (${python_version},2)
	PYTHON := python
else
	PYTHON := python2
endif

QUIPPY_SRC_DIR := ../../quippy/

FPP = gfortran -E -cpp
CFLAGS = -fPIC

LIBATOMS_SOURCES = System.f95 Units.f95 linearalgebra.f95 #MPI_context.f95 Quaternions.f95 \
 					 Dictionary.f95 Table.f95 PeriodicTable.f95 Atoms_types.f95 Atoms.f95 Connection.f95 DynamicalSystem.f95 \
 					 	clusters.f95 Structures.f95 DomainDecomposition.f95 CInOutput.f95 ParamReader.f95 Spline.f95 \
 					 frametools.f95 Topology.f95 find_surface_atoms.f95 ringstat.f95 angular_functions.f95 \
 					 steinhardt_nelson_qw.f95 nye_tensor.f95
LIBATOMS_FILES = $(addprefix ../../src/libAtoms/,${LIBATOMS_SOURCES})

WRAP_SOURCES = ${LIBATOMS_SOURCES}
WRAP_FILES = ${LIBATOMS_FILES}

WRAP_FPP_FILES = $(addsuffix .fpp,$(basename ${WRAP_SOURCES}))

F90WRAP_FILES = $(addsuffix .f90,$(addprefix f90wrap_,$(basename ${WRAP_SOURCES}))) #f90wrap_toplevel.f90

F2PY_LINK_ARGS = $(shell python -c 'import sys; print " ".join([arg for arg in sys.argv if arg.startswith("-l") or arg.startswith("-L")])' ${SYSLIBS} ${QUIPPY_LDFLAGS})

###

.PHONY: build clean install

build: _quippy.so

clean:
	-rm _castep.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -r castep
	-rm -r f90wrap_*.o
	-rm -rf src.*

${WRAP_FPP_FILES}: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do \
	  ${FPP} $$f > $$(basename $${f%.*}).fpp ; \
	done

_quippy.so: ${WRAP_FPP_FILES}
	-rm -r quippy
	mkdir quippy
	cp -r ${QUIPPY_SRC_DIR}/quippy/* quippy	
	f90wrap -m quippy ${WRAP_FPP_FILES} -P \
	   -a c_error_abort --init-file $(QUIPPY_SRC_DIR)/init.py --move-methods --shorten-routine-names \
		 --conf-file ${QUIPPY_SRC_DIR}/f90wrap_conf.py
	f2py-f90wrap --build-dir . -c -m _quippy ${F90WRAP_FILES} ${F2PY_LINK_ARGS} -L. -lquip --link-lapack

install:
	cp -r quippy ${QUIPPY_INSTALL_DIR}
	cp _quippy.so ${QUIPPY_INSTALL_DIR}
