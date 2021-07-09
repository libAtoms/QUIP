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
        python-rdkit \
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
        nglview \
        RISE \
        pandas \
        phonopy \
        torch\
        torchani\
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

# Relevant for Both Julia Environments:
# Add kernelspecs to global Jupyter
RUN mv /root/.local/share/jupyter/kernels/julia* /usr/local/share/jupyter/kernels/

########################################################################################
#RUN pip install julia
#RUN python -c 'import julia; julia.install()'
#RUN pip install git+https://github.com/casv2/pyjulip.git

CMD bash
