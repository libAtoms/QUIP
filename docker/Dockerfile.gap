# Base Python image has most up to date Python parts
FROM libatomsquip/quip-base

MAINTAINER Tom Daff "tdd20@cam.ac.uk"

# To build within the image without additonal libraries use
# the git+VANILLA version
# RUN git clone https://github.com/libAtoms/QUIP.git /opt/QUIP
# ENV BUILD DOCKER_NOGAP
ENV BUILD DOCKER_ALL
ADD . /opt/QUIP
ENV QUIP_ARCH linux_x86_64_gfortran_openmp

RUN cd /opt/QUIP \
    && mkdir -p build/${QUIP_ARCH} \
    && cp tests/rules/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make \
    && make install-quippy

ENV PATH="/opt/QUIP/build/${QUIP_ARCH}:${PATH}"

CMD jupyter notebook --port=8899 --ip="*" --allow-root --NotebookApp.token='' --NotebookApp.password=''

ENTRYPOINT /bin/bash

EXPOSE 8899
