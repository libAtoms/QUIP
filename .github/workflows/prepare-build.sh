if [ "$(uname)" == "Darwin" ]; then

    # taken from https://github.com/MacPython/gfortran-install/blob/master/gfortran_utils.sh#L97
    function install_arm64_cross_gfortran {
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
	if [[ "${PLAT:-}" == "arm64" ]]; then
            export FC=$FC_ARM64
	    export F90=$FC
	    export F95=$FC
	    export F77=$FC
	fi
    }

    if [ "$PLAT" == "arm64" ]; then
	export QUIP_ARCH=darwin_arm64_gfortran_openmp
	install_arm64_cross_gfortran
    else
	export QUIP_ARCH=darwin_x86_64_gfortran_openmp
    fi
  
else
    export QUIP_ARCH=linux_x86_64_gfortran_openmp
fi

WORK_DIR=$(dirname $0)
BUILD_DIR=$PWD/build/${QUIP_ARCH}

[[ -d ${BUILD_DIR} ]] || mkdir -p ${BUILD_DIR}
cp $WORK_DIR/Makefile.inc ${BUILD_DIR}/Makefile.inc

export NPY_DISTUTILS_APPEND_FLAGS=1
(cd ${BUILD_DIR}/../.. && make)

# if we're building a release then use tag name as version
if [[ -f GITHUB_TAG ]]; then
    cat GITHUB_TAG > ${BUILD_DIR}/VERSION
fi

# get ready to run `pip wheel` in build directory
cp ${BUILD_DIR}/../../README.md ${BUILD_DIR}
cp ${BUILD_DIR}/../../quippy/setup.py ${BUILD_DIR}

# include `quip` and `gap_fit` command line tools
cp ${BUILD_DIR}/quip ${BUILD_DIR}/quippy
cp ${BUILD_DIR}/gap_fit ${BUILD_DIR}/quippy/

# Python build dependencies
pip install oldest-supported-numpy
