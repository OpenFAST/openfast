#*********************************************************************************************************************************
# LICENSING
# Copyright (C) 2017-2025  National Renewable Energy Laboratory
# Copyright (C) 2025-2026  National Laboratory of the Rockies
#
#    This file is part of OpenFAST.
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
#*********************************************************************************************************************************

add_subdirectory(openfast)
add_subdirectory(labview)

if(BUILD_OPENFAST_CPP_API)
   add_subdirectory(openfast-cpp)
endif()

if(BUILD_FASTFARM)
   add_subdirectory(fast-farm)
endif()

if(BUILD_OPENFAST_SIMULINK_API)
   add_subdirectory(simulink)
endif()