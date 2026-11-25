#
# Copyright 2026 National Renewable Energy Laboratory
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

# Only used with LLVM Flang (added from modules/nwtc-library/CMakeLists.txt).
#
# FFTPACK must be compiled with default (4-byte) REAL; see the FFTPACK_SOURCES notes in
# modules/nwtc-library/CMakeLists.txt. LLVM Flang has no -fno-default-real-8, so the
# promotion flags are stripped from CMAKE_Fortran_FLAGS in this directory scope instead.
# fftpack_kind.f must be compiled with the same flags as fftpack5.1.f so NWTC_FFTPACK can
# verify the real kind at run time.

string(REPLACE "-fdefault-real-8"   "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
string(REPLACE "-fdefault-double-8" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

add_library(nwtc_fftpack_r4 OBJECT
   fftpack5.1.f
   fftpack_kind.f
)
set_target_properties(nwtc_fftpack_r4 PROPERTIES POSITION_INDEPENDENT_CODE ON)
