echo "prepare-build.sh received environment ARCHS=${ARCHS} QUIP_ARCH=${QUIP_ARCH} RUNNER_OS=${RUNNER_OS}"

# Install Openblas -- adapted from https://github.com/numpy/numpy/blob/main/tools/wheels/cibw_before_build.sh
echo "Installing OpenBLAS..."

if [[ "${RUNNER_OS}" == "Linux" ]]; then
	basedir=$(python .github/workflows/openblas_support.py)
	cp -r $basedir/lib/* /usr/local/lib
	cp $basedir/include/* /usr/local/include
	# Set PKG_CONFIG_PATH for meson to find openblas
	export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:${PKG_CONFIG_PATH}
elif [[ "${RUNNER_OS}" == "macOS" ]]; then
	if [[ "$ARCHS" == "arm64" ]]; then
		basedir=$(python .github/workflows/openblas_support.py)
		cp -r $basedir/lib/* /usr/local/lib
		cp $basedir/include/* /usr/local/include
		sudo mkdir -p /opt/arm64-builds/lib /opt/arm64-builds/include
		sudo chown -R $USER /opt/arm64-builds
		cp -r $basedir/lib/* /opt/arm64-builds/lib
		cp $basedir/include/* /opt/arm64-builds/include
		# Set PKG_CONFIG_PATH for meson
		export PKG_CONFIG_PATH=/opt/arm64-builds/lib/pkgconfig:${PKG_CONFIG_PATH}
	else
		brew install openblas
		brew link --force openblas
		# Set PKG_CONFIG_PATH for meson
		export PKG_CONFIG_PATH=$(brew --prefix openblas)/lib/pkgconfig:${PKG_CONFIG_PATH}
	fi
fi

# Python build dependencies
pip install meson-python meson ninja f90wrap numpy

echo Building QUIP with meson
# Build QUIP libraries first
meson setup builddir --buildtype=release -Dgap=true
meson compile -C builddir

# if we're building a release then use tag name as version
if [[ -f GITHUB_TAG ]]; then
    cat GITHUB_TAG > quippy/VERSION
fi
