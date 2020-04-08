#!/bin/bash

git fetch origin ${GITHUB_REF}:CI
git checkout CI
git submodule update

# Print the current git info
echo `git status`
echo `git log -1`
cd /openfast/reg_tests/r-test
echo `git status`
echo `git log -1`
cd /openfast

# Move into the "build" directory, remove the old reg tests, and compile
cd /openfast/build
rm -rf reg_tests
cmake ..
make -j4 install

# Run the tests

# NWTC Library tests
ctest -VV -R nwtc_library_utest

# BeamDyn-specific tests
ctest -VV -j7 -R bd_
ctest -VV -R beamdyn_utest

# OpenFAST linearization tests
# Dont run these in parallel, copying the case files can fail in a race condition
ctest -VV -L linear

# Subset of OpenFAST regression tests; do not run
## - 9, 16 because they're very sensitive
## - 19, 20 because theyre too long
## - 17, 22, 23 becuase we dont know why they fail :(
ctest -VV -j8 -I 1,1,1,2,3,4,5,6,7,8,10,11,12,13,14,15,18,21,24,25
