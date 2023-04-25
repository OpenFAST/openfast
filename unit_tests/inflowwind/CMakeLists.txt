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

set_source_files_properties(${pfunit_directory}/include/driver.F90 PROPERTIES GENERATED 1)

set(module_name "inflowwind")
set(module_directory "inflowwind")
set(module_library "ifwlib")

file(MAKE_DIRECTORY ${build_testdirectory}/${module_directory})
file(WRITE ${build_testdirectory}/${module_directory}/testSuites.inc "")

include_directories(${build_testdirectory}/${module_directory})

set(testlist
    ifw_test_tools
    test_steady_wind
    test_turbsim_wind
    test_bladed_wind
    test_hawc_wind
    test_outputs
    test_uniform_wind
)
foreach(test ${testlist})
    set(test_dependency pfunit ${source_modulesdirectory}/${module_directory}/tests/${test}.F90)
    add_custom_command(
        OUTPUT ${build_testdirectory}/${module_directory}/${test}.F90
        COMMAND ${Python_EXECUTABLE} ${pfunit_directory}/bin/pFUnitParser.py ${source_modulesdirectory}/${module_directory}/tests/${test}.F90 ${build_testdirectory}/${module_directory}/${test}.F90
        DEPENDS ${test_dependency}
    )
    set(test_sources ${test_sources} ${build_testdirectory}/${module_directory}/${test}.F90)
    file(APPEND ${build_testdirectory}/${module_directory}/testSuites.inc "ADD_TEST_SUITE(${test}_suite)\n")
endforeach()
  
add_executable(
    ${module_name}_utest
    ${pfunit_directory}/include/driver.F90
    ${test_sources}
)

target_link_libraries(
    ${module_name}_utest
    pfunit_lib
    ${module_library}
)

add_test(
    ${module_name}_utest
    ${PROJECT_BINARY_DIR}/${module_directory}/${module_name}_utest
)
