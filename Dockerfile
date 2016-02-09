FROM ubuntu:14.04
MAINTAINER Matei David <matei@cs.toronto.edu>

RUN apt-get update && apt-get install -y build-essential cmake git libhdf5-dev

#RUN mkdir /src && \
#     cd /src && \
#     git clone --recursive https://github.com/jts/nanocall.git
ADD src /src/nanocall/src

RUN mkdir /src/nanocall/build && \
    cd /src/nanocall/build && \
    cmake ../src && \
    make && \
    make install

VOLUME /data
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/nanocall"]
CMD ["--version"]
