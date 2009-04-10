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

cd ../Makefiles
rm -f Makefile.arch
ln -s Makefile.$arch Makefile.arch
cd -

mkdir -p build.$arch

if [[ ! -e build.$arch/extern ]]; then
  ln -s ../extern.$arch build.$arch/extern
fi

echo $PWD
ln -s ../Makefile build.$arch/Makefile

#export QUIP_PATH=`pwd`/
BUILDDIR=build.$arch make -C build.$arch -I../../Makefiles $*
