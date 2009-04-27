all: alltests

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

alltests:
ifeq (${QUIP_ARCH},)
	@echo
	@echo "You need to define the architecture using the ARCH variable"
	@echo
	@exit 1
endif
	QUIP_ARCH=${QUIP_ARCH} ./test.py -r *.test
