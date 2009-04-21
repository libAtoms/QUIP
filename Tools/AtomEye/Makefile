ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
endif

all:
	$(MAKE) -f Makefile.atomeye -I ${QUIP_ROOT}/Makefiles -I ${QUIP_ROOT}/build.${QUIP_ARCH} all

.DEFAULT:
	$(MAKE) -f Makefile.atomeye -I ${QUIP_ROOT}/Makefiles -I ${QUIP_ROOT}/build.${QUIP_ARCH} $@
