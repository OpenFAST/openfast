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

set(NWTCLIBS_SOURCES
  src/ModMesh.f90
  src/ModMesh_Mapping.f90
  src/ModMesh_Types.f90
  src/NWTC_Base.f90
  src/NWTC_IO.f90
  src/NWTC_Library.f90
  src/NWTC_Num.f90
  src/NWTC_RandomNumber.f90
  src/SingPrec.f90
  src/NWTC_Library_Types.f90

  # RanLux sources
  src/ranlux/RANLUX.f90
  
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
  src/NetLib/slatec/xermsg.f  )

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


if(BUILD_OPENFAST_SIMULINK_API)

  if (APPLE OR UNIX)
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
      set(NWTC_SYS_FILE src/SysMatlabLinuxGnu.f90)
    elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
      set(NWTC_SYS_FILE src/SysMatlabLinuxIntel.f90)
    endif()
  elseif (WIN32)
    set(NWTC_SYS_FILE src/SysMatlabWindows.f90)
  endif()

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  if (WIN32)
    set(NWTC_SYS_FILE src/SysGnuWin.f90)
  elseif (APPLE OR UNIX OR CYGWIN)
    set(NWTC_SYS_FILE src/SysGnuLinux.f90)
  endif ()
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
  if (APPLE OR UNIX)
    set(NWTC_SYS_FILE src/SysIFL.f90)
  elseif (WIN32)
    set(NWTC_SYS_FILE src/SysIVF.f90)
  endif (APPLE OR UNIX)
endif ()

if (NWTC_SYS_FILE)
  message("-- Setting system file as: ${NWTC_SYS_FILE}")
  list(APPEND NWTCLIBS_SOURCES ${NWTC_SYS_FILE})
else (NWTC_SYS_FILE)
  message(FATAL_ERROR "Cannot determine system file used with NWTC_Library")
endif (NWTC_SYS_FILE)

add_library(nwtclibs ${NWTCLIBS_SOURCES})

if(BUILD_OPENFAST_SIMULINK_API)
  find_package(Matlab REQUIRED)
  target_link_libraries(nwtclibs versioninfolib ${LAPACK_LIBRARIES} ${CMAKE_DL_LIBS} ${Matlab_MEX_LIBRARY})
else()
  target_link_libraries(nwtclibs versioninfolib ${LAPACK_LIBRARIES} ${CMAKE_DL_LIBS})
endif()

install(TARGETS nwtclibs
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
