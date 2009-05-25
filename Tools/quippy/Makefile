ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
endif

.PHONY: build

build:
	@echo "QUIP_ROOT=${QUIP_ROOT}"
	$(MAKE) -f Makefile.quippy -I ${QUIP_ROOT}/Makefiles -I ${QUIP_ROOT}/build.${QUIP_ARCH} build QUIP_ROOT=${QUIP_ROOT}

.DEFAULT:
	$(MAKE) -f Makefile.quippy -I ${QUIP_ROOT}/Makefiles -I ${QUIP_ROOT}/build.${QUIP_ARCH} $@ QUIP_ROOT=${QUIP_ROOT}
