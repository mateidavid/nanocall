FROM debian:stable
MAINTAINER Matei David <matei.david.at.oicr.on.ca>
ARG DEBIAN_FRONTEND=noninteractive

ADD lddtree.tgz /

# use host timezone
ENV TZ=${TZ}
RUN ln -snf /usr/share/zoneinfo/${TZ} /etc/localtime && echo ${TZ} > /etc/timezone

# use host id
RUN groupadd --gid ${GROUP_ID} ${GROUP_NAME}
RUN useradd --create-home --uid ${USER_ID} --gid ${GROUP_ID} ${USER_NAME}
USER ${USER_NAME}

VOLUME ["/data"]
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/nanocall"]
CMD ["--version"]
