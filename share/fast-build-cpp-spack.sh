#!/bin/bash

set -ex

COMPILER=gcc
SPACK_ROOT=
SPACK_EXE=${SPACK_ROOT}/bin/spack
module purge
module load gcc/5.2.0
module load python/2.7.8
module use ${SPACK_ROOT}/share/spack/modules/$(${SPACK_EXE} arch)
module load $(${SPACK_EXE} module find cmake %${COMPILER})
module load $(${SPACK_EXE} module find openmpi %${COMPILER})
module load $(${SPACK_EXE} module find hdf5 %${COMPILER})
module load $(${SPACK_EXE} module find zlib %${COMPILER})
module load $(${SPACK_EXE} module find libxml2 %${COMPILER})
module load $(${SPACK_EXE} module find xz %${COMPILER})
module load $(${SPACK_EXE} module find binutils %${COMPILER})
module list
which cmake


OPENFAST_DIR=
yaml_install_dir=`${SPACK_EXE} location -i yaml-cpp %${COMPILER}`
hdf5_install_dir=`${SPACK_EXE} location -i hdf5 %${COMPILER}`
zlib_install_dir=`${SPACK_EXE} location -i zlib %${COMPILER}`
libxml2_install_dir=`${SPACK_EXE} location -i libxml2 %${COMPILER}`
CC=gcc CXX=g++ FC=gfortran cmake \
   -DCMAKE_INSTALL_PREFIX=${OPENFAST_DIR}/install/ \
   -DCMAKE_BUILD_TYPE=DEBUG \
   -DBUILD_OPENFAST_CPP_API=ON \
   -DFPE_TRAP_ENABLED:BOOL=ON \
   -DYAML_ROOT:PATH=$yaml_install_dir \
   -DHDF5_USE_STATIC_LIBRARIES=ON \
   -DHDF5_ROOT:PATH=$hdf5_install_dir \
   -DLIBXML2_ROOT:PATH=$libxml2_install_dir \
   -DLIBXML2_USE_STATIC_LIBRARIES=ON \
   -DHDF5_ROOT:PATH=$hdf5_install_dir \
   $EXTRA_ARGS \
../ &> log.cmake
make VERBOSE=1 &> log.make
