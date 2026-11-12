#
# Copyright 2026 National Laboratory of the Rockies
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

# check_openblas_threading()
#
# OpenMP builds call LAPACK from several threads at once (TurbSim does this for
# the Cholesky factorizations). Whether that works with OpenBLAS depends on how
# the OpenBLAS library was built, so when the BLAS/LAPACK libraries found by
# CMake are OpenBLAS, ask the library and warn about the configurations that
# give wrong results or poor performance.
#
# Call this after find_package(LAPACK), and only for OpenMP builds.
#
# The check looks at the library found at configure time. If a different
# OpenBLAS is loaded at run time (for example through LD_LIBRARY_PATH or the
# Debian/Ubuntu alternatives system), this check does not see it.
function(check_openblas_threading)

  if(CMAKE_CROSSCOMPILING)
    return()
  endif()

  # openblas_get_parallel() exists only in OpenBLAS and returns how the library
  # was built: 0 = sequential, 1 = pthreads, 2 = OpenMP.
  # The result is cached together with the libraries it was obtained for, so the
  # test program is built and run again only when the BLAS/LAPACK selection changes.
  if(NOT "${OPENBLAS_PARALLEL_CHECKED_LIBS}" STREQUAL "${LAPACK_LIBRARIES}")
    unset(OPENBLAS_PARALLEL_COMPILES CACHE)
    unset(OPENBLAS_PARALLEL_RUN CACHE)
  endif()
  if(NOT DEFINED OPENBLAS_PARALLEL_COMPILES)
    set(_src "${CMAKE_BINARY_DIR}/CMakeFiles/CheckOpenBLASThreading/openblas_parallel.f90")
    file(WRITE "${_src}" "
program openblas_parallel
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   interface
      function openblas_get_parallel() bind(C, name='openblas_get_parallel')
         import :: c_int
         integer(c_int) :: openblas_get_parallel
      end function openblas_get_parallel
      subroutine c_exit(status) bind(C, name='exit')
         import :: c_int
         integer(c_int), value :: status
      end subroutine c_exit
   end interface
   call c_exit(openblas_get_parallel())
end program openblas_parallel
")
    try_run(OPENBLAS_PARALLEL_RUN OPENBLAS_PARALLEL_COMPILES
      "${CMAKE_BINARY_DIR}/CMakeFiles/CheckOpenBLASThreading"
      "${_src}"
      LINK_LIBRARIES ${LAPACK_LIBRARIES}
    )
    set(OPENBLAS_PARALLEL_CHECKED_LIBS "${LAPACK_LIBRARIES}" CACHE INTERNAL "BLAS/LAPACK libraries the OpenBLAS threading check was run for")
  endif()

  if(NOT OPENBLAS_PARALLEL_COMPILES)
    return()   # BLAS/LAPACK is not OpenBLAS; nothing to check
  endif()

  if("${OPENBLAS_PARALLEL_RUN}" STREQUAL "0")
    message(WARNING
      "OpenMP is enabled and the OpenBLAS library found is a single-threaded (sequential) build:\n"
      "  ${LAPACK_LIBRARIES}\n"
      "A single-threaded OpenBLAS is safe to call from several OpenMP threads only if it was "
      "built with USE_LOCKING=1, and CMake cannot tell whether it was. Without it, programs that "
      "call LAPACK from an OpenMP loop (TurbSim) can give wrong results or fail at random with "
      "errors such as 'LAPACK_SPOTRF: Leading minor of order N of A is not positive definite'.\n"
      "Use a pthreads or OpenMP build of OpenBLAS instead, or another BLAS/LAPACK library "
      "(see BLA_VENDOR).")
  elseif("${OPENBLAS_PARALLEL_RUN}" STREQUAL "1")
    message(WARNING
      "OpenMP is enabled and the OpenBLAS library found is a pthreads build:\n"
      "  ${LAPACK_LIBRARIES}\n"
      "Set the environment variable OPENBLAS_NUM_THREADS=1 when running the OpenMP executables. "
      "Otherwise OpenBLAS starts its own threads inside each OpenMP thread, which oversubscribes "
      "the processor (often much slower) and makes results differ slightly with the number of "
      "threads.")
  elseif("${OPENBLAS_PARALLEL_RUN}" STREQUAL "2")
    message(STATUS "OpenBLAS is an OpenMP build: safe to call from OpenMP threads")
  else()
    message(STATUS "Could not determine how OpenBLAS was built (openblas_get_parallel check returned '${OPENBLAS_PARALLEL_RUN}')")
  endif()

endfunction()
