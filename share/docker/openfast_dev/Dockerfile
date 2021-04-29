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

FROM rafmudaf/openfast-ubuntu:v2.3.0

# Move into the openfast directory and update
WORKDIR /openfast
RUN git fetch && \
    git checkout -b dev origin/dev && \
    git submodule update

# Move into the "build" directory, remove the old reg tests, and compile
WORKDIR /openfast/build
RUN rm -rf reg_tests && \
    cmake .. && \
    make -j4 install
