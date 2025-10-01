# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.13-slim AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&& apt-get install -y procps zstd \
&& rm -rf /var/lib/apt/lists/*

RUN pip install \
        'bioframe==0.8.*' \
        'pandas>2' \
        'zstandard'

RUN python3 -c 'import bioframe,pandas,zstandard'

CMD ["/bin/bash"]
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.source='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-py-utils}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
