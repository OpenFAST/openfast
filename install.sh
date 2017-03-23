#!/bin/bash

passFail() {
    if [ $1 -eq 0 ]
    then
	echo "... PASSED"
    else
	echo "... FAILED"
    fi
}

prepareSourceMods() {
#Prepare sourceMods.sh
    echo -n "Sourcing modules"
    echo -n $PWD > /tmp/fastDir
    echo "export COMPILER=${COMPILER}" > sourceMods.sh
    echo "export BUILD=${BUILD}" >> sourceMods.sh
    echo "export LAPACK=${LAPACK}" >> sourceMods.sh
    echo "COMPILER=${COMPILER}" > sourceMods.sh
    echo "BUILD=${BUILD}" >> sourceMods.sh
    echo "LAPACK=${LAPACK}" >> sourceMods.sh
    sed -e "s/FASTDIR/$(sed 's:/:\\/:g' /tmp/fastDir)/" .sourceMods.sh >> sourceMods.sh
    source sourceMods.sh
    module list
    passFail $?
    
    echo -n "Setting up lib, include and bin directories"
    [ -d install/lib ] || mkdir install/lib
    [ -d install/bin ] || mkdir install/bin
    [ -d install/include ] || mkdir install/include
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
    cd ../../../
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
    if [ "${COMPILER}" == 'gnu' ] ; then
	CC=/usr/local/Cellar/gcc/6.3.0_1/bin/gcc-6 CXX=/usr/local/Cellar/gcc/6.3.0_1/bin/g++-6 cmake ../yaml-cpp/ -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_INSTALL_PREFIX=../../../install &> log.cmake
    elif [ "${COMPILER}" == 'intel' ] ; then
	cmake ../yaml-cpp/ -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_INSTALL_PREFIX=../../../install &> log.cmake
    elif [ "${COMPILER}" == 'intelPhi' ] ; then
	cmake ../yaml-cpp/ -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS="-std=c++11 -mmic" -DCMAKE_INSTALL_PREFIX=../../../install &> log.cmake
    fi
    passFail $? 
    echo -n "   Compiling"
    make -j 8 &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ../../../
}

compileHDF5() {
    echo "Compiling hdf5"
    echo -n "   Getting source"
    [ -d modules-ext/hdf5 ] || mkdir modules-ext/hdf5
    cd modules-ext/hdf5
    curl -k -o hdf5-1.8.18.tar.bz2 https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.18.tar.bz2 &> log.curl
    passFail $?
    echo -n "   Setting up build directory"
    tar -jxf hdf5-1.8.18.tar.bz2 &> log.untar
    [ -d build ] || mkdir build &> log.mkdirBuild
    cd build
    passFail $?
    echo -n "   Configuring"
    if [ "${COMPILER}" == 'gnu' ] ; then
	cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_INSTALL_PREFIX=../../../install/ -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF ../hdf5-1.8.18/ &> log.cmake
    elif [ "${COMPILER}" == 'intel' ] ; then
	cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc  -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_INSTALL_PREFIX=../../../install/ -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF ../ &> log.cmake
    elif [ "${COMPILER}" == 'intelPhi' ] ; then
	cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc  -DCMAKE_CXX_FLAGS="-std=c++11 -mmic" -DCMAKE_INSTALL_PREFIX=../../../install/ -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=OFF ../ &> log.cmake
    fi
    passFail $? 
    echo -n "   Compiling"
    make -j 8 &> log.make
    passFail $?
    echo -n "   Installing"
    make install &> log.makeInstall
    passFail $?
    cd ../../../
}

compileFAST() {
    #FAST
    cd build/
    cp ../doConfigOpenFAST_cpp .
    ./doConfigOpenFAST_cpp
    make
    make install
}

prepPhiEnv() {
    echo -n "Prepping phi.env"
    echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > phi.env
    echo "export PATH=$PATH" >> phi.env
    echo "export FAST=FAST_ProgC_glin64" >> phi.env
    passFail $?
}


#prepareSourceMods
# if [ "${LAPACK}" == 'lapack' ]; then
#     if [ "${COMPILER}" == 'gnu' ] ; then
# 	compileLapack
#     else
# 	echo  "Can't use lapack with Intel compilers. Please use mkl instead."
# 	exit 1
#     fi
# fi
compileYAMLcpp
compileHDF5
compileFAST
if [ "${COMPILER}" == 'intelPhi' ]; then
    prepPhiEnv
fi

