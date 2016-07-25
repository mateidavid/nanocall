FROM debian:unstable
MAINTAINER Matei David <matei.david.at.oicr.on.ca>

# install prerequisites
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        cmake \
        git \
        libhdf5-dev \
        pax-utils

# if necessary, specify compiler
#RUN apt-get install -y g++-4.9 g++-5 g++-6
#ENV CC=gcc-4.9
#ENV CXX=g++-4.9

# build master
RUN cd /tmp && \
    git clone --recursive --branch master --depth 1 https://github.com/mateidavid/nanocall.git && \
    mkdir -p /tmp/nanocall/build && \
    cd /tmp/nanocall/build && \
    cmake ../src && \
    make && \
    make install

ADD get-nanocall-lddtree /usr/local/bin
