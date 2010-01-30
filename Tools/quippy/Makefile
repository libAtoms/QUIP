ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
endif

.PHONY: build install clean test

build:
	@echo "QUIP_ROOT=${QUIP_ROOT}"
	cd ${QUIP_ROOT}/Tools/quippy; python setup.py build

install:
	@echo "QUIP_ROOT=${QUIP_ROOT}"
	cd ${QUIP_ROOT}/Tools/quippy; python setup.py install

clean:
	@echo "QUIP_ROOT=${QUIP_ROOT}"
	cd ${QUIP_ROOT}/Tools/quippy; python setup.py clean

test:
	@echo "QUIP_ROOT=${QUIP_ROOT}"
	cd ${QUIP_ROOT}/Tools/quippy; python setup.py test


