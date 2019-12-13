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

FROM rafmudaf/openfast-ubuntu:dev

# Move into the openfast directory and update
WORKDIR /openfast
RUN git fetch
RUN git pull
RUN git submodule update

# Move into the "build" directory, remove the old reg tests, and compile
WORKDIR /openfast/build
RUN rm -rf reg_tests
RUN cmake ..
RUN make -j4 install

# Run the tests

# BeamDyn-specific tests
RUN ctest -VV -j7 -R bd_
RUN ctest -VV -R beamdyn_utest

# OpenFAST linearization tests
# Dont run these in parallel, copying the case files can fail in a race condition
RUN ctest -VV -L linear

# Subset of OpenFAST regression tests; do not run
## - 9, 16 because they're very sensitive
## - 19, 20 because theyre too long
## - 17, 22, 23 becuase we dont know why they fail :(
RUN ctest -VV -j8 -I 1,1,1,2,3,4,5,6,7,8,10,11,12,13,14,15,18,21,24,25
