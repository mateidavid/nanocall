FROM ubuntu:14.04
WORKDIR /
RUN apt-get update && apt-get install -y build-essential cmake git libhdf5-dev
RUN git clone --recursive https://github.com/jts/nanocall.git
RUN mkdir /build
WORKDIR /build
RUN cmake /nanocall/src && make
ENTRYPOINT /build/nanocall/nanocall
