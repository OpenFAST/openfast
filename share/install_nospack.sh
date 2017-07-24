#!/bin/bash

passFail() {
    if [ $1 -eq 0 ]
    then
	echo "... PASSED"
    else
	echo "... FAILED"
    fi
}

prepInstall() {
#Prepare for installation
    echo -n "Prepping install"
    echo -n $PWD > /tmp/fastDir
    OPENFAST_DIR=$(sed 's:/:\\/:g' /tmp/fastDir)
    echo "export OPENFAST_DIR=${OPENFAST_DIR}" > .prepInstall
    source .prepInstall
    passFail $?
}

compileLapack() {
#Registry
    echo -n "Compiling Lapack"
    [ -d modules-ext/lapack ] || mkdir modules-ext/lapack
    cd modules-ext/lapack
    curl -k -o lapack-3.6.0.tgz http://www.netlib.org/lapack/lapack-3.6.0.tgz &> log.wget
    tar -zxf lapack-3.6.0.tgz &> log.untar
    [ -d build ] || mkdir build
    cd build
    make clean &> /dev/null
    # if [ -f CMakeCache.txt] ; then
    #     rm CMakeCache.txt
    # fi
    cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_FORTRAN_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../../../install/ -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON -DBUILD_DEPRECATED=ON -DLAPACKE=ON ../lapack-3.6.0 &> log.config
    make -j 8 &> log.make
    make install &> log.makeInstall
    passFail $?
    cd ${OPENFAST_DIR}
}

compileYAMLcpp() {
#yaml-cpp
    echo "Compiling yaml-cpp"
    echo -n "   Setting up build directory"
    [ -d modules-ext/yaml-cpp ] || mkdir modules-ext/yaml-cpp
    cd modules-ext/yaml-cpp
    git clone https://github.com/jbeder/yaml-cpp.git &> /dev/null
    [ -d build ] || mkdir build
    cd build
    if [ -f CMakeCache.txt ] ; then
        rm -rf CMakeCache.txt CMakeFiles/
    fi
    passFail $?
    echo -n "   Configuring"
    cmake ../yaml-cpp/ -DCMAKE_INSTALL_PREFIX=${OPENFAST_DIR}/install &> log.cmake
    passFail $? 
    echo -n "   Compiling"
    make -j 8 &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ${OPENFAST_DIR}
}

compileHDF5() {
    echo "Compiling hdf5"
    echo -n "   Getting source"
    [ -d modules-ext/hdf5 ] || mkdir modules-ext/hdf5
    cd modules-ext/hdf5
    wget --no-check-certificate https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.19.tar.bz2 &> log.wget
    passFail $?
    echo -n "   Setting up build directory"
    tar -jxf hdf5-1.8.19.tar.bz2 &> log.untar
    cd hdf5-1.8.19
    passFail $?
    echo -n "   Configuring"
    ./configure CC=mpicc FC=mpif90 CXX=mpicxx --enable-parallel --prefix=${OPENFAST_DIR}/install &> log.config
    passFail $? 
    echo -n "   Compiling"
    make -j 8 &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ${OPENFAST_DIR}
}

compileOpenFAST() {
    echo "Compiling OpenFAST"
    echo -n "   Setting up config"
    [ -d build ] || mkdir build
    cd build/
    make clean &> /dev/null
    # if [ -f CMakeCache.txt] ; then
    #     rm CMakeCache.txt
    # fi
    cp ${OPENFAST_DIR}/share/doConfigOpenFAST_cpp_nospack .
    passFail $?
    echo -n "   Configuring"
    ./doConfigOpenFAST_cpp_nospack
    passFail $?
    echo -n "   Compiling"
    make &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ${OPENFAST_DIR}
}

prepPhiEnv() {
    echo -n "Prepping phi.env"
    echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > phi.env
    echo "export PATH=$PATH" >> phi.env
    echo "export FAST=FAST_ProgC_glin64" >> phi.env
    passFail $?
}

prepInstall
#compileYAMLcpp
#compileHDF5
compileOpenFAST
#rm .prepInstall


