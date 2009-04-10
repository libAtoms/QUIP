include Makefile.common
include Makefile.arch

LIBATOMS_DIR=../libAtoms
QUIP_DIR=../QUIP
HAVE_QUIP=1
ATOMEYE_DIR=../AtomEye

SETUP_OPTS=config_fc --f90flags="-ffree-line-length-none -x f95-cpp-input" --cpp="cpp -x f95-cpp-input" build_ext --fcompiler=gnu95 build_src -L/usr/lib64/lam/lib --f2py-opts=--debug-capi
INSTALL_DIR=${HOME}

LIBATOMS=--libatoms-dir=${LIBATOMS_DIR}
ifeq (${HAVE_QUIP},1)
  QUIP=--quip-dir=${QUIP_DIR}
endif
ATOMEYE=--atomeye-dir=${ATOMEYE_DIR}

.PHONY: all clean build install

all: build

clean:
	rm -rf build

build:
	python setup.py ${SETUP_OPTS} ${ARCH_NETCDF_INCLUDEDIR} ${ARCH_NETCDF_LIBDIR} ${ARCH_NETCDF_LIBS} ${ARCH_FOX_INCLUDEDIR} ${ARCH_FOX_LIBDIR} ${ARCH_FOX_LIBS} ${LIBATOMS} ${QUIP} ${ATOMEYE} build

install:
	python setup.py ${SETUP_OPTS} ${ARCH_NETCDF_INCLUDEDIR} ${ARCH_NETCDF_LIBDIR} ${ARCH_NETCDF_LIBS} ${ARCH_FOX_INCLUDEDIR} ${ARCH_FOX_LIBDIR} ${ARCH_FOX_LIBS} ${LIBATOMS} ${QUIP} ${ATOMEYE} install --home=${INSTALL_DIR}
