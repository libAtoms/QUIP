# Base Python (debian) image
FROM python:3.7
LABEL maintainer="Gabor Csanyi <gc121@cam.ac.uk>"

# make /bin/sh symlink to bash instead of dash:
RUN echo "dash dash/sh boolean false" | debconf-set-selections
ENV DEBIAN_FRONTEND noninteractive
RUN dpkg-reconfigure dash

#################################################################################
## Root environment
#################################################################################

RUN cp /etc/skel/.bash* /etc/skel/.profile /root/ \
 && echo "PS1='docker:\W$ '" >> /root/.bashrc


#################################################################################
## OS level deps
#################################################################################

# Build tools and deps for QUIP. Followed by some useful utilities
RUN apt-get -y update \
 && apt-get upgrade -y \
 && apt-get install -y \
        gfortran \
        cmake \
        openmpi-bin \
        libopenmpi-dev \
        liblapack-dev \
        libblas-dev \
        libnetcdf-dev \
        netcdf-bin \
        # using libzmq3-dev instead of libzmq3, this one works
        libzmq3-dev \
        # Useful tools
        pandoc \
        curl \
        vim \
        emacs-nox \
        less \
        bsdmainutils \
        man-db \
        # AtomEye
        libxpm-dev \
        libgsl0-dev \
        xterm \
        # amber
        csh \
        flex \
        # gpaw
        libxc-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Custom compilation of OpenBLAS with OpenMP enabled (linear algebra is limited to single core in debs)
# NUM_THREADS must be set otherwise docker build non-parallel version.
RUN git clone -b 'v0.3.9' --single-branch --depth 1 https://github.com/xianyi/OpenBLAS.git /tmp/OpenBLAS \
 && cd /tmp/OpenBLAS/ \
 && make NO_AFFINITY=1 USE_OPENMP=1 NUM_THREADS=32 > /dev/null \
 && make NO_AFFINITY=1 USE_OPENMP=1 NUM_THREADS=32 install > /dev/null \
 && rm -rf /tmp/OpenBLAS

RUN ln -s /opt/OpenBLAS/lib/libopenblas.so /usr/local/lib/libopenblas.so
RUN ln -s /opt/OpenBLAS/lib/libopenblas.a /usr/local/lib/libopenblas.a

# Make OpenBLAS the default
#RUN update-alternatives --install /usr/lib/libblas.so libblas.so /opt/OpenBLAS/lib/libopenblas.so 1000 \
# && update-alternatives --install /usr/lib/libblas.so.3 libblas.so.3 /opt/OpenBLAS/lib/libopenblas.so 1000 \
# && update-alternatives --install /usr/lib/liblapack.so liblapack.so /opt/OpenBLAS/lib/libopenblas.so 1000 \
# && update-alternatives --install /usr/lib/liblapack.so.3 liblapack.so.3 /opt/OpenBLAS/lib/libopenblas.so 1000 \
# && ldconfig

#################################################################################
## Python 3 environment
#################################################################################

# Put any Python libraries here
RUN pip install --upgrade pip \
 && pip install --no-cache-dir \
        jupyter \
        numpy \
        scipy \
        scikit-learn\
        scikit-build\
        matplotlib \
        pyamg \
        imolecule \
        sphinx \
        sphinx-rtd-theme \
        nbsphinx \
        numpydoc \
        spglib \
        rdkit-pypi \
        nglview \
        RISE \
        pandas \
        phonopy \
        torch \
        torchani \
        gpaw \
 && jupyter nbextension enable --py --sys-prefix widgetsnbextension \
 && jupyter nbextension enable --py --sys-prefix nglview \
 && jupyter-nbextension install rise --py --sys-prefix \
 && jupyter-nbextension enable rise --py --sys-prefix

# get latest version of ASE
RUN pip install git+https://gitlab.com/ase/ase

# Keep the source for examples
RUN git clone https://github.com/libAtoms/matscipy.git /opt/matscipy \
 && cd /opt/matscipy \
 && pip install --no-cache-dir .


#################################################################################
## Julia v1.5.0
#################################################################################

# specify paths for Julia 1.6
ENV JULIA_PATH /opt/julia/v1.6.1
# PKG_DIR is now replaced with DEPOT_PATH
ENV JULIA_DEPOT_PATH /opt/julia/share/site
# Add to path as current version
ENV PATH ${JULIA_PATH}/bin:$PATH

RUN mkdir -p ${JULIA_PATH} \
 && cd ${JULIA_PATH} \
 && curl --location "https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.1-linux-x86_64.tar.gz" | tar xz --strip-components 1

RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("registry add https://github.com/JuliaRegistries/General.git")'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("registry add https://github.com/JuliaMolSim/MolSim.git")'

RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add ACE@0.8.4"); using ACE'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add JuLIP"); using JuLIP'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add ASE"); using ASE'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add IPFitting"); using IPFitting'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add IJulia"); using IJulia'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add Images"); using Images'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add ImageMagick"); using ImageMagick'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add Plots"); using Plots'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add CSV"); using CSV'
RUN julia -e 'using Pkg;Pkg.REPLMode.pkgstr("add DataFrames"); using DataFrames'


# Relevant for Both Julia Environments:
# Add kernelspecs to global Jupyter
RUN mv /root/.local/share/jupyter/kernels/julia* /usr/local/share/jupyter/kernels/

#################################################################################
#RUN pip install julia
#RUN python -c 'import julia; julia.install()'
#RUN pip install git+https://github.com/casv2/pyjulip.git

#################################################################################
# End of the quip-base image
#
# This was the very base env, now we are building QUIP and GAP
#
#
#
#################################################################################

# All the QUIPs go here; added to path in the end.
ENV QUIP_ROOT /opt/quip
ENV QUIP_INSTALLDIR ${QUIP_ROOT}/bin
ENV PATH ${QUIP_INSTALLDIR}:${PATH}

WORKDIR /opt/quip

RUN git clone --depth 1 --recursive https://github.com/libAtoms/QUIP.git .

#################################################################################
# Quip library for LAMMPS
#################################################################################

# lammps should be linked with SERIAL version of QUIP other configurations are
# untested and too complicated for a user (mixed paralleisms).
ENV QUIP_ARCH linux_x86_64_gfortran
ADD Makefile.inc build/${QUIP_ARCH}/Makefile.inc

# Build only libquip for serial to keep a slim image.
# Makefile.inc is also required to compile lammps.
RUN make libquip > /dev/null \
 && find build/${QUIP_ARCH} -type f ! \( -name 'libquip.a' -o -name 'Makefile.inc' \) -delete


#################################################################################
# Quip with OpenMP
#################################################################################

ENV QUIP_ARCH linux_x86_64_gfortran_openmp
ADD Makefile.inc build/${QUIP_ARCH}/Makefile.inc

RUN make \
 && make install \
 && make install-quippy

#################################################################################
# End of the quip-gap image
#
# Now we are doing the general changes for the final image
#
#
#
#################################################################################

#################################################################################
## Additional python packages
#################################################################################

# Atomistica
#   we need to set the arg mismatch flag for the compiler, their setup will not add it automatically
#   should remove this when this is done: https://github.com/Atomistica/atomistica/pull/40
WORKDIR /opt
RUN git clone --depth 1 https://github.com/Atomistica/atomistica.git \
 && cd atomistica \
 && sed -e 's/f77flags=/f77flags=-fallow-argument-mismatch /g' -i setup.cfg.gnu\
 && python setup.py install

# AtomEye needs to link with QUIP for xyz read-write
WORKDIR ${QUIP_ROOT}/src/AtomEye
ENV PATH ${QUIP_ROOT}/src/AtomEye/bin:${PATH}

RUN git clone --depth 1 https://github.com/jameskermode/AtomEye.git . \
 && make > /dev/null \
 && cd Python \
 && python setup.py install

#################################################################################
## Softwares
#################################################################################

# LAMMPS compilation
WORKDIR /opt/lammps

ENV LAMMPS_PATH /opt/lammps
ENV PATH ${LAMMPS_PATH}/src/:${PATH}

# Build LAMMPS from github - faster and more reliable than from the website
RUN git clone --branch stable_29Sep2021_update3 --depth 1 https://github.com/lammps/lammps.git .

# Build `shlib` objects first so they have `-fPIC` then symlink the directory
# so they can be reused to build the binaries halving the compilation time.
# Clean up Obj files immedaitely to keep image smaller.
RUN cd src \
 && make yes-all \
 && make no-lib \
 && QUIP_ARCH=linux_x86_64_gfortran; make yes-ml-quip \
 && make no-intel \
 && make yes-python \
 && make -j4 mpi mode=shlib \
 && make install-python \
 && ln -s Obj_shared_mpi Obj_mpi \
 && make -j4 mpi \
 && make clean-all

##################################### 
## RUST
#####################################

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs > rust.sh
RUN chmod u+rx rust.sh
RUN ./rust.sh -y
ENV PATH /root/.cargo/bin:$PATH

RUN pip install git+https://github.com/Luthaf/rascaline.git

#################################################################################
## Data
#################################################################################

# Published GAPs
#WORKDIR /opt/share/potentials

# Remote URLs do not get decompressed so pipe through tar
#ENV POTENTIALS_DIR /opt/share/potentials

#ADD files/GAPPotentials.md .

#RUN wget -nv -O- "http://www.libatoms.org/pub/Home/TungstenGAP/GAP_6.tbz2" \
#    | tar xj -P --transform "s,^,${POTENTIALS_DIR}/GAP/Tungsten/," \
# && wget -nv -O- "http://www.libatoms.org/pub/Home/IronGAP/gp33b.tar.gz" \
#    | tar xz -P --transform "s,^,${POTENTIALS_DIR}/GAP/Iron/," \
# && wget -nv -O- "http://www.libatoms.org/pub/Home/DataRepository/gap_dft_corrections_water.tgz" \
#    | tar xz -P --transform "s,^,${POTENTIALS_DIR}/GAP/Water/," \
# && wget -nv -O- "http://www.libatoms.org/pub/Home/DataRepository/gap_dft_corrections_ch4_h2o.tgz" \
#    | tar xz -P --transform "s,^,${POTENTIALS_DIR}/GAP/WaterCH4/," \
# && wget -nv -O- "http://www.libatoms.org/pub/Home/DataRepository/gap_dft_1_2_body_LiH2O.tgz" \
#    | tar xz -P --transform "s,^,${POTENTIALS_DIR}/GAP/WaterLiH2O/," \
# && wget -nv -O- "http://www.libatoms.org/pub/Home/DataRepository/aC_GAP.tar.gz" \
#    | tar xz -P --transform "s,^,${POTENTIALS_DIR}/GAP/Carbon/,"

# GPAW data
WORKDIR /opt/share/gpaw

# Ensure we don't run interactively
ENV GPAW_SETUP_VERSION 0.9.20000
ENV GPAW_SETUP_PATH /opt/share/gpaw/gpaw-setups-${GPAW_SETUP_VERSION}

RUN wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.20000.tar.gz -O - | tar -xz  -C .

#################################################################################
## Entrypoint
#################################################################################

WORKDIR /root/

# Set bash instead of sh as defualt terminal in Jupyter
ENV SHELL /bin/bash

# Launch in the home directory of the user
ADD files/demo.ipynb .

CMD bash -c exit && jupyter notebook --ip=$(hostname -i) --port=8899 --allow-root
EXPOSE 8899
