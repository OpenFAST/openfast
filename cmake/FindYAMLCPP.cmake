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

find_path(YAML_INCLUDES
  yaml-cpp/yaml.h
  HINTS ${YAML_ROOT} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES include)

find_library(YAML_LIBRARIES
  NAMES yaml-cpp libyaml-cpp.a
  HINTS ${YAML_ROOT} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAMLCPP DEFAULT_MSG YAML_INCLUDES YAML_LIBRARIES)
mark_as_advanced(YAML_INCLUDES YAML_LIBRARIES)
