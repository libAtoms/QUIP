echo "prepare-build.sh received environment ARCHS=${ARCHS} QUIP_ARCH=${QUIP_ARCH} RUNNER_OS=${RUNNER_OS}"

if [[ "${RUNNER_OS}" == "macOS" && "$ARCHS" == "arm64" ]] then

	echo "Installing arm64 cross compiler..."

    # taken from https://github.com/MacPython/gfortran-install/blob/master/gfortran_utils.sh#L97
	curl -L -O https://github.com/isuruf/gcc/releases/download/gcc-10-arm-20210228/gfortran-darwin-arm64.tar.gz
	export GFORTRAN_SHA=f26990f6f08e19b2ec150b9da9d59bd0558261dd
	if [[ "$(shasum gfortran-darwin-arm64.tar.gz)" != "${GFORTRAN_SHA}  gfortran-darwin-arm64.tar.gz" ]]; then
            echo "shasum mismatch for gfortran-darwin-arm64"
            exit 1
	fi
	sudo mkdir -p /opt/
	sudo cp "gfortran-darwin-arm64.tar.gz" /opt/gfortran-darwin-arm64.tar.gz
	pushd /opt
        sudo tar -xvf gfortran-darwin-arm64.tar.gz
        sudo rm gfortran-darwin-arm64.tar.gz
	popd
	export FC_ARM64="$(find /opt/gfortran-darwin-arm64/bin -name "*-gfortran")"
	local libgfortran="$(find /opt/gfortran-darwin-arm64/lib -name libgfortran.dylib)"
	local libdir=$(dirname $libgfortran)

	export FC_ARM64_LDFLAGS="-L$libdir -Wl,-rpath,$libdir"
	if [[ "${ARCHS:-}" == "arm64" ]]; then
            export FC=$FC_ARM64
	    export F90=$FC
	    export F95=$FC
	    export F77=$FC
	fi
fi

# Install Openblas -- adapted from https://github.com/numpy/numpy/blob/main/tools/wheels/cibw_before_build.sh
echo "Installing OpenBLAS..."

basedir=$(python .github/workflows/openblas_support.py)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include
if [[ "${RUNNER_OS}" == "macOS" && "$ARCHS" == "arm64" ]]; then
    sudo mkdir -p /opt/arm64-builds/lib /opt/arm64-builds/include
    sudo chown -R $USER /opt/arm64-builds
    cp -r $basedir/lib/* /opt/arm64-builds/lib
    cp $basedir/include/* /opt/arm64-builds/include
fi

WORK_DIR=$(dirname $0)
BUILDDIR=$PWD/build/${QUIP_ARCH}

[[ -d ${BUILDDIR} ]] || mkdir -p ${BUILDDIR}
cp $WORK_DIR/Makefile.${QUIP_ARCH}.inc ${BUILDDIR}/Makefile.inc

export NPY_DISTUTILS_APPEND_FLAGS=1

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
