# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
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

# include other makefiles and export env variables
export QUIP_ARCH

ifeq (${QUIP_ROOT},)
QUIP_ROOT=${PWD}/..
endif
include Makefile.${QUIP_ARCH}
include ${BUILDDIR}/Makefile.inc
include Makefile.rules

# Test if the 'python' in the path is a Python 2, if not then
# explicitly call python2
# This is a workaround for old virtualenv versions that don't create the
# 'python2' symlink.
python_version := $(word 2, $(subst ., , $(shell python -V 2>&1)))
ifeq (${python_version},2)
	PYTHON := python
else
	PYTHON := python2
endif

.PHONY: build install clean test

build:
	LD=${LINKER} CC=${CC} F90=${F90} F77=${F77} QUIP_ROOT=${QUIP_ROOT} ${PYTHON} setup.py build

install:
	LD=${LINKER} CC=${CC} F90=${F90} F77=${F77} QUIP_ROOT=${QUIP_ROOT} ${PYTHON} setup.py install ${QUIPPY_INSTALL_OPTIONS}

clean:
	LD=${LINKER} CC=${CC} F90=${F90} F77=${F77} QUIP_ROOT=${QUIP_ROOT} ${PYTHON} setup.py clean

test:
	LD=${LINKER} CC=${CC} F90=${F90} F77=${F77} QUIP_ROOT=${QUIP_ROOT} ${PYTHON} setup.py test ${QUIPPY_TEST_OPTIONS}
