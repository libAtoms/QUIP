.PHONY: all alltests clean

all: alltests

alltests:
ifeq (${QUIP_ARCH},)
	@echo
	@echo "You need to define the architecture using the ARCH variable"
	@echo
	@exit 1
endif
	QUIP_ARCH=${QUIP_ARCH} ./test.py -r *.test

clean:
	-rm *.candidate *.reference *.stdin *.stdout *.xyz *.xml *.new
