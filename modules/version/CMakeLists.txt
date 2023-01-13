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

# Allow the git-version information to be set externally for the
# package manager distributions (homebrew, conda, apt) that download
# the tarball source code rather than clone with git.
if( DEFINED GIT_DESCRIBE )
  message( WARNING
    "Version information has been set as a CMake flag. This should only used when the git-version cannot be set automatically."
  )
else()
  include(GetGitRevisionDescription)
  git_describe(GIT_DESCRIBE)
endif()
add_definitions(-DGIT_VERSION_INFO="${GIT_DESCRIBE}")

add_library(versioninfolib_obj OBJECT src/VersionInfo.f90)
target_link_libraries(versioninfolib_obj nwtclibs_obj)

add_library(versioninfolib $<TARGET_OBJECTS:versioninfolib_obj>)
target_link_libraries(versioninfolib nwtclibs)

install(TARGETS versioninfolib versioninfolib_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
