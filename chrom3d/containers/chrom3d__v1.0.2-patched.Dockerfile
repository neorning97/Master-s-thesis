# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ghcr.io/paulsengroup/ci-docker-images/ubuntu-24.04-cxx-clang-20 AS builder

RUN apt-get update \
&& apt-get install -y curl \
                      patch \
                      libboost-dev \
                      libboost-filesystem-dev \
                      libboost-random-dev

ARG CHROM3D_VER=1.0.2

RUN cd /tmp/ \
&& curl -LO "https://github.com/Chrom3D/Chrom3D/archive/refs/tags/v${CHROM3D_VER}.tar.gz" \
&& echo "f19e046216f17c43369ccddf45cdc114b364cf3ad7173501a56ab87a7d310be9  v${CHROM3D_VER}.tar.gz" > checksum.sha256 \
&& sha256sum -c checksum.sha256 \
&& tar -xf *.tar.gz

COPY containers/patches/chrom3d-*.patch /tmp/

RUN cd "/tmp/Chrom3D-${CHROM3D_VER}" \
&& patch -p1 -i /tmp/chrom3d-cmake.patch \
&& patch -p1 -i /tmp/chrom3d-fix-warnings.patch

RUN cmake -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_C_COMPILER=clang \
          -DCMAKE_CXX_COMPILER=clang++ \
          -DCMAKE_INSTALL_PREFIX=/tmp/chrom3d-staging \
          -S "/tmp/Chrom3D-${CHROM3D_VER}/" \
          -B /tmp/build/ \
&& cmake --build /tmp/build/ -j $(nproc) \
&& cmake --install /tmp/build/

FROM ubuntu:24.04 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

COPY --from=builder "/tmp/chrom3d-staging/" "/usr/local/"

CMD ["/usr/local/bin/Chrom3D"]
WORKDIR /data

RUN Chrom3D --version


LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.documentation='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.source='https://github.com/robomics/chrom3d-nf'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-chrom3d}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
