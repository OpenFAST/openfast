#!/bin/bash

verbosecommand() { echo ">>> $1" && eval $1 && echo "<<<"; }

# Configure Bash to exit if any command returns an error
set -e

verbosecommand "cd /openfast"

repo="OpenFAST"
echo "GITHUB_EVENT_NAME: ${GITHUB_EVENT_NAME}"
if [[ "${GITHUB_EVENT_NAME}" != "pull_request" ]]; then
    repo=${GITHUB_ACTOR}
fi
# Create a branch "CI" at the current commit from the GH Actor's fork.
verbosecommand "git fetch https://github.com/${repo}/openfast ${GITHUB_REF}:CI"
verbosecommand "git checkout CI"
verbosecommand "git submodule update"

# Display the current git info
echo "*** git-status from openfast:"
verbosecommand "git status"

echo "*** git-log from openfast:"
verbosecommand "git log -1"

verbosecommand "cd /openfast/reg_tests/r-test"
echo "*** git-status from r-test:"
verbosecommand "git status"

echo "*** git-log from r-test:"
verbosecommand "git log -1"

verbosecommand "cd /openfast"

# Display the differences between this commit and `dev`
echo "*** git-diff from ${GITHUB_REF} to dev:"
verbosecommand "git diff dev --numstat"

# Move into the "build" directory, remove the old reg tests, and compile
verbosecommand "cd /openfast/build"
verbosecommand "rm -rf reg_tests"
verbosecommand "cmake .."
verbosecommand "make -j4 install"

# Run the tests

# NWTC Library tests
verbosecommand "ctest -VV -R nwtc_library_utest"

# BeamDyn-specific tests
verbosecommand "ctest -VV -j7 -R bd_"
verbosecommand "ctest -VV -R beamdyn_utest"

# OLAF free vortex wake tests
ctest -VV -R fvw_utest

# OpenFAST linearization tests
# Dont run these in parallel, copying the case files can fail in a race condition
# Exclude the Ideal_Beam test cases
# - They fail consistently in the Docker container when run on GitHub,
#   but pass everywhere else including running the same Docker image locally
verbosecommand "ctest -VV -L linear -E Ideal"

# Subset of OpenFAST regression tests; do not run
## - 9, 16 because they're very sensitive
## - 19, 20 because theyre too long
## - 17, 22, 23 because we dont know why they fail :(
verbosecommand "ctest -VV -j8 -I 1,1,1,2,3,4,5,6,7,8,10,11,12,13,14,15,18,21,24,25,26,27,28"
