# Base Python image has most up to date Python parts
FROM libatomsquip/quip-base-software

MAINTAINER Tom Daff "tdd20@cam.ac.uk"

# All the QUIPs go here; added to path in the end.
ENV QUIP_ROOT /opt/quip

# To build within the image without additonal libraries use
# the git+VANILLA version
# RUN git clone https://github.com/libAtoms/QUIP.git ${QUIP_ROOT}
# ENV BUILD NOGAP
ENV BUILD GAP
# ENV BUILD ALL
ADD . ${QUIP_ROOT}

# LAMMPS compilation

# lammps should be linked with SERIAL version of QUIP
# other configurations are untested and too complicated
# for a user (mixed paralleisms).
ENV QUIP_ARCH linux_x86_64_gfortran
ENV LAMMPS_PATH /opt/lammps

# Build only libquip for serial to keep a slim image.
# Makefile.inc is also required to compile lammps.
RUN cd ${QUIP_ROOT} \
    && mkdir -p build/${QUIP_ARCH} \
    && cp docker/arch/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make libquip > /dev/null \
    && find build/${QUIP_ARCH} -type f ! \( -name 'libquip.a' -o -name 'Makefile.inc' \) -delete

# TODO: prune any unwanted directories in this command
RUN mkdir -p ${LAMMPS_PATH} \
    && cd ${LAMMPS_PATH} \
    && curl http://lammps.sandia.gov/tars/lammps-stable.tar.gz | tar xz --strip-components 1

# Build `shlib` objects first so they have `-fPIC` then symlink the directory
# so they can be reused to build the binaries halving the compilation time.
# Clean up Obj files immedaitely to keep image smaller.
RUN cd ${LAMMPS_PATH}/src \
    && make yes-all \
    && make no-lib \
    && make yes-user-quip yes-python no-user-intel \
    && make mpi mode=shlib \
    && make install-python \
    && ln -s Obj_shared_mpi Obj_mpi \
    && make mpi \
    && make clean-all

ENV PATH ${LAMMPS_PATH}/src/:${PATH}


# MPI QUIP for parallel within QUIP
# Installs with _mpi suffix, e.g. quip_mpi
ENV QUIP_ARCH linux_x86_64_gfortran_openmpi

RUN cd ${QUIP_ROOT} \
    && mkdir -p build/${QUIP_ARCH} \
    && cp docker/arch/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make > /dev/null \
    && QUIP_INSTALLDIR=${QUIP_ROOT}/bin make install \
    && find build/${QUIP_ARCH} -type f ! \( -name 'libquip.a' -o -name 'Makefile.inc' \) -delete


# QUIP for general use is the OpenMP version.
# Installs with no suffix, e.g. quip
# Also installs quippy
# Keep all libraries for AtomEye
ENV QUIP_ARCH linux_x86_64_gfortran_openmp

RUN cd ${QUIP_ROOT} \
    && mkdir -p build/${QUIP_ARCH} \
    && cp docker/arch/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make > /dev/null \
    && QUIP_INSTALLDIR=${QUIP_ROOT}/bin make install \
    && make install-quippy > /dev/null \
    && find build/${QUIP_ARCH} -type f ! \( -name '*.a' -o -name 'Makefile.inc' \) -delete

# AtomEye needs to link with QUIP for xyz read-write
RUN git clone --depth 1 https://github.com/jameskermode/AtomEye.git ${QUIP_ROOT}/src/AtomEye \
    && cd ${QUIP_ROOT}/src/AtomEye \
    && make \
    && cd ${QUIP_ROOT}/src/AtomEye/Python \
    && python setup.py install

ENV PATH ${QUIP_ROOT}/bin:${QUIP_ROOT}/src/AtomEye/bin:${PATH}

# ENTRYPOINT ["/bin/bash", "-c"]

# Public GAP image requires license agreement
# Replace bash with a license check script
RUN mkdir -p /bin/real/ \
    && mv /bin/bash /bin/real/bash

ADD docker/files/fakebash /bin/bash

CMD bash -c exit && jupyter notebook --port=8899 --ip='*' --allow-root --NotebookApp.token='' --NotebookApp.password=''

EXPOSE 8899
