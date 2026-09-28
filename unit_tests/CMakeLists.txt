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

add_executable(aerodyn_inflow_utest ${PROJECT_SOURCE_DIR}/modules/aerodyn/tests/test_interpolated_inflow.F90)
target_link_libraries(aerodyn_inflow_utest aerodyn_driver_subs testdrivelib)
set(inflow_test_data_dir "${PROJECT_BINARY_DIR}/reg_tests/modules/aerodyn")
file(COPY "${PROJECT_SOURCE_DIR}/modules/aerodyn/tests/interpolated_inflow/time"
     "${PROJECT_SOURCE_DIR}/modules/aerodyn/tests/interpolated_inflow/steady"
     DESTINATION "${inflow_test_data_dir}")
foreach(inflow_case time time_inputs steady null)
  if (inflow_case STREQUAL "steady")
    set(inflow_fixture "steady")
  else()
    set(inflow_fixture "time")
  endif()
  add_test(NAME aerodyn_utest_inflow_${inflow_case}
           COMMAND ${CMAKE_COMMAND} -E env AD_INFLOW_CASE=${inflow_case}
                   $<TARGET_FILE:aerodyn_inflow_utest>
                   "${inflow_test_data_dir}/${inflow_fixture}/ad_driver.dvr")
  set_tests_properties(aerodyn_utest_inflow_${inflow_case} PROPERTIES RUN_SERIAL TRUE)
endforeach()

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
  ${PROJECT_SOURCE_DIR}/modules/inflowwind/tests/test_grid3d_field.F90
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
  ${PROJECT_SOURCE_DIR}/modules/nwtc-library/tests/test_NWTC_C_Binding.F90
  ${PROJECT_SOURCE_DIR}/modules/nwtc-library/tests/test_NWTC_FFTPACK.F90
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

add_dependencies(unit_tests aerodyn_inflow_utest)

if (AMREX_READER AND BUILD_FASTFARM)
  # Need input files
  add_custom_target(copy_amrex_reader_utest_files ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/modules/awae/tests/data ${CMAKE_BINARY_DIR}/unit_tests/data
  )
  # FAST.Farm AMReX reader test
  add_executable(amrex_reader_utest
    ${CMAKE_SOURCE_DIR}/modules/awae/tests/AMReX_reader_utest.F90
    ${CMAKE_SOURCE_DIR}/modules/awae/tests/test_AMReX_reader.f90
  )
  target_link_libraries(amrex_reader_utest awaelib nwtclibs testdrivelib)
  add_test(NAME amrex_reader_utest COMMAND amrex_reader_utest)
  add_dependencies(amrex_reader_utest copy_amrex_reader_utest_files)     # copy input files needed for this test
  add_dependencies(unit_tests amrex_reader_utest)     # add to the unit_tests target
endif()
