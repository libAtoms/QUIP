if [ -n "$IS_MACOS" ]; then
  if [ "$PLAT" == "arm64" ]; then
    export QUIP_ARCH=darwin_arm64_gfortran_openmp
  else
    export QUIP_ARCH=darwin_x86_64_gfortran_openmp
  fi
else
  export QUIP_ARCH=linux_x86_64_gfortran_openmp
fi

pip install oldest-supported-numpy

[[ -d ${REPO_DIR} ]] || mkdir -p ${REPO_DIR}
cp Makefile.${QUIP_ARCH}.inc ${REPO_DIR}/Makefile.inc

export NPY_DISTUTILS_APPEND_FLAGS=1
(cd ${REPO_DIR}/../.. && make quippy)

# if we're building a release then use tag name as version
if [[ -f GITHUB_TAG ]]; then
    cat GITHUB_TAG > ${REPO_DIR}/VERSION
fi

# get ready to run `pip wheel` in build directory
cp ${REPO_DIR}/../../README.md ${REPO_DIR}
cp ${REPO_DIR}/../../quippy/setup.py ${REPO_DIR}

# include `quip` and `gap_fit` command line tools
cp ${REPO_DIR}/quip ${REPO_DIR}/quippy
cp ${REPO_DIR}/gap_fit ${REPO_DIR}/quippy/
