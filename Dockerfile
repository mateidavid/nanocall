FROM debian:stable
MAINTAINER Matei David <matei.david.at.oicr.on.ca>
ARG DEBIAN_FRONTEND=noninteractive

# install prerequisites
RUN for i in 1 2 3; do \
        apt-get update \
        && break; sleep 1; \
    done && \
    for i in 1 2 3; do \
        apt-get install -y \
             build-essential \
             cmake \
             libhdf5-dev \
        && break; sleep 1; \
    done

# if necessary, specify compiler
#RUN apt-get install -y g++-4.9 g++-5 g++-6
#ENV CC=gcc-4.9
#ENV CXX=g++-4.9

# add source
ADD . /src/

# build and install
RUN mkdir -p /src/build && \
    cd /src/build && \
    cmake ../src && \
    make && \
    make install

VOLUME ["/data"]
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/nanocall"]
CMD ["--version"]
