#! /bin/bash

# EXIT if anything oes wrong so
# we don't end up pushing empty
# repositories

set -e

# builds expect this
export QUIP_ROOT=`pwd`

# packages for building docs
pip install sphinx sphinx_rtd_theme numpydoc

# needed to nbconvert ipynb files and to process the rst files
pip install nbconvert\[execute\] ipython

# quippy is working, install it
make install-quippy

## Install atomeye from src
cd src
git clone https://github.com/jameskermode/AtomEye.git
cd AtomEye/Python

# Install in the virtualenv
python setup.py install

# Work from the docs directory
DOCS_DIR=${QUIP_ROOT}/doc
cd ${DOCS_DIR}

# Put a working copy of the gh-pages where they are expected
PAGES_DIR=${DOCS_DIR}/_build/html
git clone -b gh-pages ${PAGES_URL} ${PAGES_DIR} > /dev/null 2>&1 || echo "Failed to clone docs"

# Clean out previous builds; should keep .git and .nojekll but get rid
# of cruft
cd ${PAGES_DIR}
rm -rf *

# html version goes into _build
cd ${DOCS_DIR}
# Force complete rebuild
O=-Ea make html

# set up git so it can push
git config --global user.name "Travis-CI"
git config --global user.email "build@travis.org"

# add -A will both add and delete files
# These commands are moved from docpush/Makefile
cd ${PAGES_DIR}
git add -A
git commit -m docpush
git push origin gh-pages --quiet

