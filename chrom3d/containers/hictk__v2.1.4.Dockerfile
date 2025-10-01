# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

ARG CONTAINER_VERSION

FROM paulsengroup/hictk:${CONTAINER_VERSION} AS base

ARG CONTAINER_TITLE
ARG CONTAINER_VERSION

RUN apt-get update \
&& apt-get install -y procps \
&& rm -rf /var/lib/apt/lists/*

RUN hictk --help

ENTRYPOINT []
CMD ["/bin/bash"]
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.source='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hictk}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
