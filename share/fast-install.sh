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
    openfast_dir=$(sed 's:/:\\/:g' /tmp/fastDir)
    echo "export openfast_dir=${openfast_dir}" > .prepInstall
    source .prepInstall
    passFail $?
}

compileLapack() {
#Registry
    echo -n "Compiling Lapack"
    [ -d modules/lapack ] || mkdir modules/lapack
    cd modules/lapack
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
    cd ${openfast_dir}
}

compileYAMLcpp() {
#yaml-cpp
    echo "Compiling yaml-cpp"
    echo -n "   Setting up build directory"
    [ -d modules/yaml-cpp ] || mkdir modules/yaml-cpp
    cd modules/yaml-cpp
    git clone https://github.com/jbeder/yaml-cpp.git &> /dev/null
    [ -d build ] || mkdir build
    cd build
    if [ -f CMakeCache.txt ] ; then
        rm -rf CMakeCache.txt CMakeFiles/
    fi
    passFail $?
    echo -n "   Configuring"
    cmake ../yaml-cpp/ -DCMAKE_INSTALL_PREFIX=${openfast_dir}/install &> log.cmake
    passFail $? 
    echo -n "   Compiling"
    make -j 8 &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ${openfast_dir}
}

compileHDF5() {
    echo "Compiling hdf5"
    echo -n "   Getting source"
    [ -d modules/hdf5 ] || mkdir modules/hdf5
    cd modules/hdf5
    wget --no-check-certificate https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.bz2 &> log.wget
    passFail $?
    echo -n "   Setting up build directory"
    tar -jxf hdf5-1.10.1.tar.bz2 &> log.untar
    cd hdf5-1.10.1
    passFail $?
    echo -n "   Configuring"
    ./configure CC=mpicc FC=mpif90 CXX=mpicxx --enable-parallel --prefix=${openfast_dir}/install &> log.config
    passFail $? 
    echo -n "   Compiling"
    make -j 8 &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ${openfast_dir}
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
    cp ${openfast_dir}/share/fast-build-cpp.sh .
    passFail $?
    echo -n "   Configuring"
    ./fast-build-cpp.sh
    passFail $?
    echo -n "   Compiling"
    make &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ${openfast_dir}
}

prepInstall
compileYAMLcpp
compileHDF5
compileOpenFAST


