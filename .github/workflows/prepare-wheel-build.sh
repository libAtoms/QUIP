echo "prepare-build.sh received environment ARCHS=${ARCHS} QUIP_ARCH=${QUIP_ARCH} RUNNER_OS=${RUNNER_OS}"

# Install Openblas -- adapted from https://github.com/numpy/numpy/blob/main/tools/wheels/cibw_before_build.sh
echo "Installing OpenBLAS..."

if [[ "${RUNNER_OS" == "Linux" ]]; then
	basedir=$(python .github/workflows/openblas_support.py)
	cp -r $basedir/lib/* /usr/local/lib
	cp $basedir/include/* /usr/local/include
elif [[ "${RUNNER_OS} == "macOS" ]]; then
	if [[ "$ARCHS" == "arm64" ]]; then
		basedir=$(python .github/workflows/openblas_support.py)
		cp -r $basedir/lib/* /usr/local/lib
		cp $basedir/include/* /usr/local/include
		sudo mkdir -p /opt/arm64-builds/lib /opt/arm64-builds/include
		sudo chown -R $USER /opt/arm64-builds
		cp -r $basedir/lib/* /opt/arm64-builds/lib
		cp $basedir/include/* /opt/arm64-builds/include
	else
		brew install openblas
		brew link --force openblas
	fi
fi

WORK_DIR=$(dirname $0)
BUILDDIR=$PWD/build/${QUIP_ARCH}

[[ -d ${BUILDDIR} ]] || mkdir -p ${BUILDDIR}
cp $WORK_DIR/Makefile.${QUIP_ARCH}.inc ${BUILDDIR}/Makefile.inc

echo Building QUIP
(cd ${BUILDDIR}/../.. && make)

echo Building quippy
(cd ${BUILDDIR}/../.. && make quippy)

# if we're building a release then use tag name as version
if [[ -f GITHUB_TAG ]]; then
    cat GITHUB_TAG > ${BUILDDIR}/VERSION
fi

# get ready to run `pip wheel` in build directory
cp ${BUILDDIR}/../../README.md ${BUILDDIR}
cp ${BUILDDIR}/../../quippy/setup.py ${BUILDDIR}

# include `quip` and `gap_fit` command line tools
cp ${BUILDDIR}/quip ${BUILDDIR}/quippy
cp ${BUILDDIR}/gap_fit ${BUILDDIR}/quippy/

# Python build dependencies
pip install oldest-supported-numpy
