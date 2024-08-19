#
# Copyright 2017 National Renewable Energy Laboratory
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

# -----------------------------------------------------------
# -- OpenFAST Unit Testing
# -----------------------------------------------------------

include(CTest)

# Test-Drive library
add_library(testdrivelib test-drive/testdrive.F90 test-drive/testdrive_version.f90)

# AeroDyn Unit Tests
add_executable(aerodyn_utest 
  ${PROJECT_SOURCE_DIR}/modules/aerodyn/tests/aerodyn_utest.F90
  ${PROJECT_SOURCE_DIR}/modules/aerodyn/tests/test_AD_FVW.F90
)
target_link_libraries(aerodyn_utest aerodynlib versioninfolib testdrivelib)
add_test(NAME aerodyn_utest COMMAND aerodyn_utest)

# BeamDyn Unit Tests
add_executable(beamdyn_utest 
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/beamdyn_utest.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_Crv.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_Misc.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_diffmtc.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_MemberEta.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_QuadraturePointData.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_ShapeFuncs.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_TrapezoidalPointWeight.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_BD_InitializeNodalLocations.F90
  ${PROJECT_SOURCE_DIR}/modules/beamdyn/tests/test_tools.F90
)
target_link_libraries(beamdyn_utest beamdynlib versioninfolib testdrivelib)
add_test(NAME beamdyn_utest COMMAND beamdyn_utest)

# InflowWind Unit Tests
add_executable(inflowwind_utest 
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/inflowwind_utest.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_bladed_wind.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_hawc_wind.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_outputs.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_steady_wind.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_turbsim_wind.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_uniform_wind.F90
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/ifw_test_tools.F90
)
target_link_libraries(inflowwind_utest ifwlib versioninfolib testdrivelib)
add_test(NAME inflowwind_utest COMMAND inflowwind_utest)

# NWTC Library Unit Tests
add_executable(nwtc_library_utest 
  ${PROJECT_SOURCE_DIR}/modules/nwtc-library/tests/nwtc_library_utest.F90
  ${PROJECT_SOURCE_DIR}/modules/nwtc-library/tests/NWTC_Library_test_tools.F90
  ${PROJECT_SOURCE_DIR}/modules/nwtc-library/tests/test_NWTC_IO_FileInfo.F90
  ${PROJECT_SOURCE_DIR}/modules/nwtc-library/tests/test_NWTC_RandomNumber.F90
)
target_link_libraries(nwtc_library_utest nwtclibs testdrivelib)

add_test(NAME nwtc_library_utest COMMAND nwtc_library_utest)

# Version Info Unit Tests
add_executable(versioninfo_utest 
  ${PROJECT_SOURCE_DIR}/modules/version/tests/versioninfo_utest.F90
  ${PROJECT_SOURCE_DIR}/modules/version/tests/VersionInfo_test_tools.F90
  ${PROJECT_SOURCE_DIR}/modules/version/tests/test_VersionInfo_CheckArgs.F90
)
target_link_libraries(versioninfo_utest versioninfolib testdrivelib)
add_test(NAME versioninfo_utest COMMAND versioninfo_utest)

# Custom target combining all the unit tests
add_custom_target(
  unit_tests
  DEPENDS 
    nwtc_library_utest 
    versioninfo_utest
    beamdyn_utest 
    aerodyn_utest 
    inflowwind_utest 
)