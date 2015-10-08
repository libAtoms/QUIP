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

all:
	QUIP_ARCH=${QUIP_ARCH} ${PYTHON} run_all.py -v
