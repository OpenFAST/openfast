# Copyright 2016-2024 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
ARG BASE=ubuntu:jammy

# Build stage 1: builds openfast.
FROM ${BASE} AS build

# Install dependencies
# For gfortran-8
# RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
# apt-get install gfortran-8

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && \
    apt-get install --no-install-recommends -y \
        software-properties-common \
        build-essential \
        cmake  \
        cmake-curses-gui \
        gcc \
        gfortran \
        git \
        libopenblas-dev \
        make \
    && rm -rf /var/lib/apt/lists/*

# Copy in the checked-out code version.
WORKDIR /openfast
COPY . .

# Build the project.
RUN mkdir build
WORKDIR /openfast/build

# NOTE: building with optimizations on (RELEASE or RELWITHDEBINFO), the virtual machine
# will require about 6GB of memory. Otherwise, the gfortran compiler will exit with an
# "internal error"
ENV FC=/usr/bin/gfortran
ARG CMAKE_OPTIONS="-DBUILD_TESTING=OFF -DBUILD_FASTFARM=ON -DDOUBLE_PRECISION=OFF -DCMAKE_BUILD_TYPE=RELEASE"
RUN cmake .. ${CMAKE_OPTIONS}

ARG BUILD_CORES=4
RUN make -j${BUILD_CORES} install

# Build stage 2: provides built openfast in a small image.
FROM ${BASE} as production
COPY --from=build /openfast/install /openfast/install

ARG TIMEZONE=UTC
ENV DEBIAN_FRONTEND=noninteractive TZ=${TIMEZONE}

RUN apt-get update -qq && \
    apt-get install --no-install-recommends -y \
        libopenblas-dev \
        libgomp1 \
        nano \
    && rm -rf /var/lib/apt/lists/*

# Make `openfast` command work.
ENV PATH=/openfast/install/bin:$PATH
