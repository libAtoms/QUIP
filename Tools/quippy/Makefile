include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

.PHONY: clean build install

quippy: build

SETUP_OPTS=config_fc --f90flags="-ffree-line-length-none -x f95-cpp-input" --cpp="cpp -x f95-cpp-input" build_ext --fcompiler=gnu95 build_src -L/usr/lib64/lam/lib #--f2py-opts=--debug-capi
INSTALL_DIR=${HOME}

space = $(empty) $(empty)

LIBATOMS_DIR=../../libAtoms
include ${LIBATOMS_DIR}/Makefile
LIBATOMS=--libatoms-dir=${LIBATOMS_DIR} \
	--libatoms-sources=$(subst $(space),:,$(strip ${LIBATOMS_F77_SOURCES})):$(subst $(space),:,$(strip ${LIBATOMS_F95_SOURCES})):$(subst $(space),:,$(strip ${LIBATOMS_C_SOURCES}))\

QUIP_CORE_DIR=../../QUIP_Core
include ${QUIP_CORE_DIR}/Makefile
QUIP_CORE=--quip-core-dir=${QUIP_CORE_DIR} \
	--quip-core-sources=$(subst $(space),:,$(strip ${TB_F77_SOURCES})):$(subst $(space),:,$(strip ${ALL_F95_FILES}))

ATOMEYE_DIR=../AtomEye
ATOMEYE=--atomeye-dir=${ATOMEYE_DIR}


build:
	python setup.py ${SETUP_OPTS} -I${NETCDF_INCDIR} -L${NETCDF_LIBDIR} ${NETCDF4_LIBS} -I${FOX_INCDIR} -L${FOX_LIBDIR} ${FOX_LIBS} ${LIBATOMS} ${QUIP_CORE} ${ATOMEYE} build --debug

install:
	python setup.py ${SETUP_OPTS} -I${NETCDF_INCDIR} -L${NETCDF_LIBDIR} ${NETCDF4_LIBS} -I${FOX_INCDIR} -L${FOX_LIBDIR} ${FOX_LIBS} ${LIBATOMS} ${QUIP_CORE} ${ATOMEYE} build --debug install --home=${INSTALL_DIR}

clean:
	rm -rf build
