# HQ

ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.inc
include Makefile.rules

QUIPPY_SRC_DIR=${QUIP_ROOT}/quippy
QUIPPY_INSTALL_DIR = ${HOME}/lib/python

LIBATOMS_SOURCES =  System ExtendableStr MPI_context Units linearalgebra Quaternions Dictionary Table PeriodicTable Atoms_types Atoms Connection DynamicalSystem clusters Structures DomainDecomposition CInOutput ParamReader Spline frametools Topology find_surface_atoms ringstat angular_functions steinhardt_nelson_qw nye_tensor

QUIPCORE_SOURCES = Functions Potential AdjustablePotential

QUIPUTILS_SOURCES = elasticity real_space_covariance ts_params transition_state crackparams cracktools

WRAP_SOURCES = ${LIBATOMS_SOURCES} ${QUIPCORE_SOURCES} ${QUIPUTILS_SOURCES}

ifeq (${HAVE_GAP},1)
WRAP_SOURCES += descriptors
endif

ifeq (${HAVE_CP2K},1)
WRAP_SOURCES += cp2k_driver_module
endif

WRAP_FPP_FILES = $(addsuffix .fpp,${WRAP_SOURCES})
F90WRAP_FILES = $(addprefix f90wrap_,$(addsuffix .f90,${WRAP_SOURCES}))
###

.PHONY: python clean install

python: _quippy.so

clean:
	-rm _quippy.so ${F90WRAP_FILES}
	-rm -r quippy
	-rm -r f90wrap_*.o
	-rm -rf src.*

LIBS =  -L. -lquiputils -lquip_core ${GP_LIB} -latoms 
ifeq (${HAVE_CP2K},1)
  LIBS += -lcp2k_driver
endif
ifeq (${HAVE_THIRDPARTY},1)
  LIBS += -lthirdparty
endif

_quippy.so: ${WRAP_FPP_FILES}
	f90wrap -v -m quippy ${WRAP_FPP_FILES} -k ${QUIPPY_SRC_DIR}/kind_map \
	   -s ${QUIPPY_SRC_DIR}/string_lengths --short-names ${QUIPPY_SRC_DIR}/short_names -P \
	   -a f90wrap_abort --init-file ${QUIPPY_SRC_DIR}/init.py
	f2py-f90wrap --build-dir . -c -m _quippy ${F90WRAP_FILES} \
	   ${LINKFLAGS} ${LIBS}

install:
	cp -r quippy/ ${QUIPPY_INSTALL_DIR}
	cp _quippy.so ${QUIPPY_INSTALL_DIR}
