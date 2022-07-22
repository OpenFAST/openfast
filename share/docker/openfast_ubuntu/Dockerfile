#
# Copyright 2016 National Renewable Energy Laboratory
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

FROM ubuntu:bionic

# Install dependencies

# For gfortran-8
# RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
# apt-get install gfortran-8

ENV DEBIAN_FRONTEND=noninteractive TZ=America/Denver

RUN apt update -qq && \
    apt install -y software-properties-common build-essential && \
    add-apt-repository ppa:git-core/ppa -y && \
    apt install -y python3-pip && \
    apt install -y cmake cmake-curses-gui && \
    apt install -y gcc gfortran make && \
    apt install -y libblas-dev liblapack-dev && \
    apt install -y git && \
    apt install -y nano

RUN pip3 install numpy

# Configure the environment
ENV FC=/usr/bin/gfortran

# Clone the project
RUN git clone --recursive https://github.com/openfast/openfast.git openfast
WORKDIR /openfast

# Build the project
RUN mkdir build
WORKDIR /openfast/build

# NOTE: building with optimizations on (RELEASE or RELWITHDEBINFO), the virtual machine
# will require about 6GB of memoery. Otherwise, the gfortran compiler will exit with an
# "internal error"
RUN cmake .. -DBUILD_TESTING=ON -DDOUBLE_PRECISION=ON -DCMAKE_BUILD_TYPE=RELWITHDEBINFO
RUN make -j4 install
