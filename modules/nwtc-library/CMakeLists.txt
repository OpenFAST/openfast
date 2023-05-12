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

#-------------------------------------------------------------------------------
# NWTC System File
#-------------------------------------------------------------------------------

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  if (WIN32)
    set(NWTC_SYS_FILE src/SysGnuWin.f90)
    set(NWTC_SYS_FILE_MATLAB ${CMAKE_CURRENT_SOURCE_DIR}/src/SysMatlabWindows.f90)
  elseif (APPLE OR UNIX OR CYGWIN)
    set(NWTC_SYS_FILE src/SysGnuLinux.f90)
    set(NWTC_SYS_FILE_MATLAB ${CMAKE_CURRENT_SOURCE_DIR}/src/SysMatlabLinuxGnu.f90)
  endif ()
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
  if (APPLE OR UNIX)
    set(NWTC_SYS_FILE src/SysIFL.f90)
    set(NWTC_SYS_FILE_MATLAB ${CMAKE_CURRENT_SOURCE_DIR}/src/SysMatlabLinuxIntel.f90)
  elseif (WIN32)
    set(NWTC_SYS_FILE src/SysIVF.f90)
    set(NWTC_SYS_FILE_MATLAB ${CMAKE_CURRENT_SOURCE_DIR}/src/SysMatlabWindows.f90)
  endif (APPLE OR UNIX)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
  if(APPLE OR UNIX)
    set(NWTC_SYS_FILE src/SysFlangLinux.f90)
  endif()
endif ()

# System file
if (NWTC_SYS_FILE)
  message("-- Setting system file as: ${NWTC_SYS_FILE}")
else ()
  message(FATAL_ERROR "Cannot determine system file used with NWTC_Library")
endif ()

#-------------------------------------------------------------------------------
# NWTC Library
#-------------------------------------------------------------------------------

set(NWTCLIBS_SOURCES

  src/NWTC_Base.f90
  src/SingPrec.f90

  src/ModMesh.f90
  src/ModMesh_Mapping.f90
  src/ModMesh_Types.f90

  src/NWTC_IO.f90
  src/NWTC_Library.f90
  src/NWTC_Num.f90
  src/NWTC_RandomNumber.f90
  
  src/NWTC_Library_Types.f90
  src/VTK.f90
  src/YAML.f90
  src/JSON.f90

  # RanLux sources
  src/ranlux/RANLUX.f90

  # Public Domain Aeronautical Software (PDAS) Polynomial Root Finder
  src/Polynomial/quartic.f90
  
  # NetLib sources
  src/NetLib/fftpack/fftpack4.1.f
  src/NetLib/scalapack/dlasrt2.f
  src/NetLib/scalapack/slasrt2.f
  src/NetLib/fftpack/NWTC_FFTPACK.f90
  src/NetLib/lapack/NWTC_LAPACK.f90
  src/NetLib/scalapack/NWTC_ScaLAPACK.f90

  # NetLib SLATEC sources
  src/NetLib/slatec/NWTC_SLATEC.f90
  src/NetLib/slatec/dqk61.f
  src/NetLib/slatec/qk61.f
  src/NetLib/slatec/d1mach.f
  src/NetLib/slatec/r1mach.f
  src/NetLib/slatec/xercnt.f
  src/NetLib/slatec/xerhlt.f
  src/NetLib/slatec/xerprn.f
  src/NetLib/slatec/xersve.f
  src/NetLib/slatec/fdump.f
  src/NetLib/slatec/i1mach.f
  src/NetLib/slatec/j4save.f
  src/NetLib/slatec/xgetua.f
  src/NetLib/slatec/xermsg.f  
)

get_filename_component(FCNAME ${CMAKE_Fortran_COMPILER} NAME)

# Recursive use of routine in qk61/dqk61 will trigger errors in debug
if (CMAKE_BUILD_TYPE MATCHES Debug)
   if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
      set_source_files_properties(src/NetLib/slatec/NWTC_SLATEC.f90 PROPERTIES COMPILE_FLAGS -frecursive)
      set_source_files_properties(src/NetLib/slatec/dqk61.f         PROPERTIES COMPILE_FLAGS -frecursive)
      set_source_files_properties(src/NetLib/slatec/qk61.f          PROPERTIES COMPILE_FLAGS -frecursive)
   elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
      if("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "19")
         if (WIN32)
            set_source_files_properties(src/NetLib/slatec/NWTC_SLATEC.f90 PROPERTIES COMPILE_FLAGS "/recursive")
            set_source_files_properties(src/NetLib/slatec/dqk61.f         PROPERTIES COMPILE_FLAGS "/recursive")
            set_source_files_properties(src/NetLib/slatec/qk61.f          PROPERTIES COMPILE_FLAGS "/recursive")
         else()
            set_source_files_properties(src/NetLib/slatec/NWTC_SLATEC.f90 PROPERTIES COMPILE_FLAGS "-recursive")
            set_source_files_properties(src/NetLib/slatec/dqk61.f         PROPERTIES COMPILE_FLAGS "-recursive")
            set_source_files_properties(src/NetLib/slatec/qk61.f          PROPERTIES COMPILE_FLAGS "-recursive")
         endif()
      else()   # The above syntax will be deprecated after version 20.  This is supported from version 19 onwards
               # see: https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/floating-point-options/recursive.html
         if (WIN32)
            set_source_files_properties(src/NetLib/slatec/NWTC_SLATEC.f90 PROPERTIES COMPILE_FLAGS "/assume:recursion")
            set_source_files_properties(src/NetLib/slatec/dqk61.f         PROPERTIES COMPILE_FLAGS "/assume:recursion")
            set_source_files_properties(src/NetLib/slatec/qk61.f          PROPERTIES COMPILE_FLAGS "/assume:recursion")
         else()
            set_source_files_properties(src/NetLib/slatec/NWTC_SLATEC.f90 PROPERTIES COMPILE_FLAGS "-assume recursion")
            set_source_files_properties(src/NetLib/slatec/dqk61.f         PROPERTIES COMPILE_FLAGS "-assume recursion")
            set_source_files_properties(src/NetLib/slatec/qk61.f          PROPERTIES COMPILE_FLAGS "-assume recursion")
         endif()
      endif()
   endif()
endif()

# Create NWTC Library
add_library(nwtclibs ${NWTC_SYS_FILE} ${NWTCLIBS_SOURCES})
target_link_libraries(nwtclibs PUBLIC
  ${LAPACK_LIBRARIES} 
  ${CMAKE_DL_LIBS} 
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)
if (USE_DLL_INTERFACE)
  target_compile_definitions(nwtclibs PRIVATE USE_DLL_INTERFACE)
endif (USE_DLL_INTERFACE)
if (USE_LOCAL_STATIC_LAPACK)
  add_dependencies(nwtclibs lapack)
endif()

# Install libraries
install(TARGETS nwtclibs
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

#-------------------------------------------------------------------------------
# NWTC Library for MATLAB
#-------------------------------------------------------------------------------

if (BUILD_OPENFAST_SIMULINK_API)
  
  # Check that matlab system file was found
  if (NWTC_SYS_FILE_MATLAB)
    message("-- Setting Matlab system file as: ${NWTC_SYS_FILE_MATLAB}")
    message(STATUS "In NWTC Matlab_MEX_LIBRARY: ${Matlab_MEX_LIBRARY}")
  else ()
    message(FATAL_ERROR "Cannot determine Matlab system file used with NWTC_Library")
  endif ()

  # NWTC Library with MATLAB system file
  add_library(nwtclibs-matlab ${NWTC_SYS_FILE_MATLAB} ${NWTCLIBS_SOURCES})
  set_target_properties(nwtclibs-matlab PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/matlab)
  if (USE_DLL_INTERFACE)
    target_compile_definitions(nwtclibs-matlab PRIVATE USE_DLL_INTERFACE)
  endif (USE_DLL_INTERFACE)
  target_link_libraries(nwtclibs-matlab PUBLIC
    ${LAPACK_LIBRARIES} 
    ${CMAKE_DL_LIBS} 
    ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
  )
  if (USE_LOCAL_STATIC_LAPACK)
    add_dependencies(nwtclibs lapack)
  endif()

  # Install matlab nwtc libraries
  install(TARGETS nwtclibs-matlab
    EXPORT "${CMAKE_PROJECT_NAME}Libraries"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endif ()

