#! /bin/bash

# builds expect this
export QUIP_ROOT=`pwd`

# packages for building docs
pip install sphinx nbconvert numpydoc

# quippy is working, install it
make install-quippy

## Install atomeye from src
cd src
git clone https://github.com/jameskermode/AtomEye.git
cd AtomEye/Python

# Install in the virtualenv
python setup.py install
cd ../../../

# Work from the docs directory
cd doc

# Put a working copy of the gh-pages where they are expected
PAGES_DIR=../../QUIP-pages
git clone -b gh-pages ${PAGES_URL} ${PAGES_DIR} > /dev/null 2>&1

# set up git so it can push
git config --global user.name "Travis-CI"
git config --global user.email "build@travis.org"

# For some reason, it won't import from the current directory;
export PYTHONPATH=`pwd`:$PYTHONPATH

# html version is fine, push it later
make docpush

