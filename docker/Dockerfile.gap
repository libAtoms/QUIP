# Base Python image has most up to date Python parts
FROM libatomsquip/quip-base

MAINTAINER Tom Daff "tdd20@cam.ac.uk"

# QUIP compilation - OpenMP version
ENV QUIP_ROOT /opt/QUIP

# To build within the image without additonal libraries use
# the git+VANILLA version
# RUN git clone https://github.com/libAtoms/QUIP.git ${QUIP_ROOT}
# ENV BUILD DOCKER_NOGAP
ENV BUILD DOCKER_ALL
ADD . ${QUIP_ROOT}
ENV QUIP_ARCH linux_x86_64_gfortran_openmp

RUN cd ${QUIP_ROOT} \
    && mkdir -p build/${QUIP_ARCH} \
    && cp docker/arch/rules/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make \
    && make install-quippy

ENV PATH ${QUIP_ROOT}/build/${QUIP_ARCH}:${PATH}

# LAMMPS compilation
# lammps should be linked with SERIAL version of QUIP
# other configurations are untested and too complicated
# for a user
ENV QUIP_ARCH linux_x86_64_gfortran
ENV LAMMPS_PATH /opt/lammps

RUN cd ${QUIP_ROOT} \
    && mkdir -p build/${QUIP_ARCH} \
    && cp docker/arch/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make libquip

RUN mkdir -p ${LAMMPS_PATH} \
    && cd ${LAMMPS_PATH} \
    && curl http://lammps.sandia.gov/tars/lammps-stable.tar.gz | tar xz --strip-components 1

RUN cd ${LAMMPS_PATH}/src \
    && make yes-all \
    && make no-lib \
    && make yes-user-quip yes-python \
    && make mpi \
    && make mpi mode=shlib \
    && make install-python

ENV PATH ${LAMMPS_PATH}/src/:${PATH}

# ENTRYPOINT ["/bin/bash", "-c"]

CMD jupyter notebook --port=8899 --ip='*' --allow-root --NotebookApp.token='' --NotebookApp.password=''

EXPOSE 8899
