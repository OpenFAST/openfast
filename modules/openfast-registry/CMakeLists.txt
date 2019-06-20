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

set(REGISTRY_SOURCES
  src/data.c
  src/gen_c_types.c
  src/gen_module_files.c
  src/misc.c
  src/my_strtok.c
  src/reg_parse.c
  src/registry.c
  src/sym.c
  src/symtab_gen.c
  src/type.c
  )

add_executable(openfast_registry ${REGISTRY_SOURCES})

set_target_properties(openfast_registry PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/modules/openfast-registry
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/modules/openfast-registry
)
