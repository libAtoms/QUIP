ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
endif

.PHONY: build

build:
	$(MAKE) -f Makefile.quippy -I ${QUIP_ROOT}/Makefiles -I ${QUIP_ROOT}/build.${QUIP_ARCH} build

.DEFAULT:
	$(MAKE) -f Makefile.quippy -I ${QUIP_ROOT}/Makefiles -I ${QUIP_ROOT}/build.${QUIP_ARCH} $@
