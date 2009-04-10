#!/bin/bash

if [[ $# < 1 ]]; then
  echo "Usage: $0 arch [ files... ]" 1>&2
  exit 1
fi

arch=$1; shift

if [[ ! -e ../Makefiles/Makefile.$arch ]]; then
  echo "Can't find Makefile.$arch" 1>&2
  ls -1 ../Makefiles/Makefile.*
  exit 2
fi

mkdir -p build.$arch

echo $PWD
if [ ! -e build.$arch/Makefile ]; then
  ln -s ../Makefile build.$arch/Makefile
fi
if [ -e Makefile.programs ]; then
  if [ ! -e build.$arch/Makefile.programs ]; then
    ln -s ../Makefile.programs build.$arch/Makefile.programs
  fi
fi

BUILDDIR=build.$arch MAKEFILE_ARCH_SUFFIX="$arch" make -C build.$arch -I../../Makefiles $*
