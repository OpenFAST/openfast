# The simple testing framework

[![License](https://img.shields.io/badge/license-MIT%7CApache%202.0-blue)](LICENSE-Apache)
[![Latest Version](https://img.shields.io/github/v/release/fortran-lang/test-drive)](https://github.com/fortran-lang/test-drive/releases/latest)
[![CI](https://github.com/fortran-lang/test-drive/workflows/CI/badge.svg)](https://github.com/fortran-lang/test-drive/actions)
[![codecov](https://codecov.io/gh/fortran-lang/test-drive/branch/main/graph/badge.svg)](https://codecov.io/gh/fortran-lang/test-drive)

This project offers a lightweight, procedural unit testing framework based on nothing but standard Fortran.
Integration with [meson](https://mesonbuild.com), [cmake](https://cmake.org) and [Fortran package manager (fpm)](https://github.com/fortran-lang/fpm) is available.
Alternatively, the [``testdrive.F90``](src/testdrive.F90) source file can be redistributed in the project's testsuite as well.


## Usage

Testsuites are defined by a ``collect_interface`` returning a set of ``unittest_type`` objects.
To create a new test use the ``new_unittest`` constructor, which requires a test identifier and a procedure with a ``test_interface`` compatible signature.
The error status is communicated by the allocation status of an ``error_type``.

The necessary boilerplate code to setup the test entry point is just

```fortran
program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use test_suite1, only : collect_suite1
  use test_suite2, only : collect_suite2
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("suite1", collect_suite1), &
    new_testsuite("suite2", collect_suite2) &
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester
```

Every test is defined in a separate module using a ``collect`` function, which is exported and added to the ``testsuites`` array in the test runner.
All tests have a simple interface with just an allocatable ``error_type`` as output to provide the test results.

```fortran
module test_suite1
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_suite1

contains

!> Collect all exported unit tests
subroutine collect_suite1(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("valid", test_valid), &
    new_unittest("invalid", test_invalid, should_fail=.true.) &
    ]

end subroutine collect_suite1

subroutine test_valid(error)
  type(error_type), allocatable, intent(out) :: error
  ! ...
end subroutine test_valid

subroutine test_invalid(error)
  type(error_type), allocatable, intent(out) :: error
  ! ...
end subroutine test_invalid

end module test_suite1
```


### Checking test conditions

The procedures defining the tests can contain any Fortran code required for checking the correctness of the project.
An easy way to do so is provided by the generic ``check`` function.

```f90
subroutine test_valid(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, 1 + 2 == 3)
  if (allocated(error)) return

  ! equivalent to the above
  call check(error, 1 + 2, 3)
  if (allocated(error)) return
end subroutine test_valid
```

After each check, the status of the error should be checked.
Uncaught errors will not be silently dropped, instead the error will be caught, its message displayed and the run aborted.
Possible ways to use check are listed below

| available checks     | arguments                                                      |
| -------------------- | -------------------------------------------------------------- |
| logical check        | *error*, *logical*, ...                                        |
| status check         | *error*, *integer*, ...                                        |
| logical comparison   | *error*, *logical*, *logical*, ...                             |
| integer comparison   | *error*, *integer*, ...                                        |
| character comparison | *error*, *character*, *character*, ...                         |
| real comparison      | *error*, *real*, *real*, ..., thr=*real*, rel=*logical*        |
| real NaN check       | *error*, *real*, ...                                           |
| complex comparison   | *error*, *complex*, *complex*, ..., thr=*real*, rel=*logical*  |
| complex NaN check    | *error*, *complex*, ...                                        |

Each check will generate a meaningful error message based on the available arguments, but can also be provided with a custom error message instead.

To generate custom checks the ``test_failed`` procedure is available to generate error messages

```f90
subroutine test_custom(error)
  type(error_type), allocatable, intent(out) :: error

  ! ...

  if (.not.cond) then
    call test_failed(error, "Custom check failed")
    return
  end if

  ! ...

  if (.not.cond) then
    call test_failed(error, "Custom check failed", "Additional context")
    return
  end if

end subroutine test_custom
```

To conditionally skip a test use the ``skip_test`` procedure.
It uses the same signature as ``test_failed``, but will mark the respective test as skipped, this is useful to disable tests based on conditional compilation, *e.g.* by using a preprocessor or a different submodule.
An uncaught skipped test will fail regardless, therefore make sure to not run any other checks afterwards.


### Integration in build systems

Finally, for usage with *fpm* it is beneficial to have a single test driver which can run all tests.
While this brings the disadvantage of always having to run the complete testsuite, the main driver can provide the flexibility to select the suite and also the unit test using the boilerplate code shown here:

```f90
!> Driver for unit testing
program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
    & select_suite, run_selected, get_argument
  use test_suite1, only : collect_suite1
  use test_suite2, only : collect_suite2
  implicit none
  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("suite1", collect_suite1), &
    new_testsuite("suite2", collect_suite2) &
    ]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
    is = select_suite(testsuites, suite_name)
    if (is > 0 .and. is <= size(testsuites)) then
      if (allocated(test_name)) then
        write(error_unit, fmt) "Suite:", testsuites(is)%name
        call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
        if (stat < 0) then
          error stop 1
        end if
      else
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end if
    else
      write(error_unit, fmt) "Available testsuites"
      do is = 1, size(testsuites)
        write(error_unit, fmt) "-", testsuites(is)%name
      end do
      error stop 1
    end if
  else
    do is = 1, size(testsuites)
      write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  end if

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

end program tester
```

From *fpm* this allows to run all tests using just the *fpm test* command, but also to debug an individual test in a debugger.
For example to run *broken-test* in *large-suite* with ``gdb`` use

```
fpm test --runner gdb -- large-suite broken-test
```

To make this approach feasible for meson the tests can be created as individual suites.
A usual layout of the test directory like

```
test
├── main.f90
├── meson.build
├── test_suite1.f90
├── test_suite2.f90
└── ...
```

Can use the following snippet to automatically create individual tests running complete suites inside the driver.
Resolution to the unit tests is possible but usually not desired, because the individual runtime of the tests will be short compared to the overhead to start the actual test.

```meson
testdrive_dep = dependency('test-drive', fallback: ['test-drive', 'testdrive_dep'])

tests = [
  'suite1',
  'suite2',
  # ...
]

test_srcs = files(
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [proj_dep, testdrive_dep],
)

test('all tests', tester)

foreach t : tests
  test(t, tester, args: t)
endforeach
```

Similar for a CMake based build the tests can be generated automatically for the layout shown below.

```
test
├── CMakeLists.txt
├── main.f90
├── test_suite1.f90
├── test_suite2.f90
└── ...
```

The CMake file in the test directory should look similar to the one shown here

```cmake
if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

# Unit testing
set(
  tests
  "suite1"
  "suite2"
)
set(
  test-srcs
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.f90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "test-drive::test-drive"
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()
```

<details>
<summary>CMake module to find testing framework</summary>

The following module allows to find or fetch an installation of this project in CMake

```cmake
#[[.rst:
Find test-drive
---------------

Makes the test-drive project available.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``test-drive::test-drive``
  The test-drive library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variables:

``TEST_DRIVE_FOUND``
  True if the test-drive library is available

``TEST_DRIVE_SOURCE_DIR``
  Path to the source directory of the test-drive project,
  only set if the project is included as source.

``TEST_DRIVE_BINARY_DIR``
  Path to the binary directory of the test-drive project,
  only set if the project is included as source.

Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``TEST_DRIVE_FIND_METHOD``
  Methods to find or make the project available. Available methods are
  - ``cmake``: Try to find via CMake config file
  - ``pkgconf``: Try to find via pkg-config file
  - ``subproject``: Use source in subprojects directory
  - ``fetch``: Fetch the source from upstream

``TEST_DRIVE_DIR``
  Used for searching the CMake config file

``TEST_DRIVE_SUBPROJECT``
  Directory to find the test-drive subproject, relative to the project root

#]]

set(_lib "test-drive")
set(_pkg "TEST_DRIVE")
set(_url "https://github.com/fortran-lang/test-drive")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

foreach(method ${${_pkg}_FIND_METHOD})
  if(TARGET "${_lib}::${_lib}")
    break()
  endif()

  if("${method}" STREQUAL "cmake")
    message(STATUS "${_lib}: Find installed package")
    if(DEFINED "${_pkg}_DIR")
      set("_${_pkg}_DIR")
      set("${_lib}_DIR" "${_pkg}_DIR")
    endif()
    find_package("${_lib}" CONFIG)
    if("${_lib}_FOUND")
      message(STATUS "${_lib}: Found installed package")
      break()
    endif()
  endif()

  if("${method}" STREQUAL "pkgconf")
    find_package(PkgConfig QUIET)
    pkg_check_modules("${_pkg}" QUIET "${_lib}")
    if("${_pkg}_FOUND")
      message(STATUS "Found ${_lib} via pkg-config")

      add_library("${_lib}::${_lib}" INTERFACE IMPORTED)
      target_link_libraries(
        "${_lib}::${_lib}"
        INTERFACE
        "${${_pkg}_LINK_LIBRARIES}"
        )
      target_include_directories(
        "${_lib}::${_lib}"
        INTERFACE
        "${${_pkg}_INCLUDE_DIRS}"
        )

      break()
    endif()
  endif()

  if("${method}" STREQUAL "subproject")
    if(NOT DEFINED "${_pkg}_SUBPROJECT")
      set("_${_pkg}_SUBPROJECT")
      set("${_pkg}_SUBPROJECT" "subprojects/${_lib}")
    endif()
    set("${_pkg}_SOURCE_DIR" "${PROJECT_SOURCE_DIR}/${${_pkg}_SUBPROJECT}")
    set("${_pkg}_BINARY_DIR" "${PROJECT_BINARY_DIR}/${${_pkg}_SUBPROJECT}")
    if(EXISTS "${${_pkg}_SOURCE_DIR}/CMakeLists.txt")
      message(STATUS "Include ${_lib} from ${${_pkg}_SUBPROJECT}")
      add_subdirectory(
        "${${_pkg}_SOURCE_DIR}"
        "${${_pkg}_BINARY_DIR}"
      )

      add_library("${_lib}::${_lib}" INTERFACE IMPORTED)
      target_link_libraries("${_lib}::${_lib}" INTERFACE "${_lib}")

      # We need the module directory in the subproject before we finish the configure stage
      if(NOT EXISTS "${${_pkg}_BINARY_DIR}/include")
        make_directory("${${_pkg}_BINARY_DIR}/include")
      endif()

      break()
    endif()
  endif()

  if("${method}" STREQUAL "fetch")
    message(STATUS "Retrieving ${_lib} from ${_url}")
    include(FetchContent)
    FetchContent_Declare(
      "${_lib}"
      GIT_REPOSITORY "${_url}"
      GIT_TAG "HEAD"
      )
    FetchContent_MakeAvailable("${_lib}")

    add_library("${_lib}::${_lib}" INTERFACE IMPORTED)
    target_link_libraries("${_lib}::${_lib}" INTERFACE "${_lib}")

    # We need the module directory in the subproject before we finish the configure stage
    FetchContent_GetProperties("${_lib}" SOURCE_DIR "${_pkg}_SOURCE_DIR")
    FetchContent_GetProperties("${_lib}" BINARY_DIR "${_pkg}_BINARY_DIR")
    if(NOT EXISTS "${${_pkg}_BINARY_DIR}/include")
      make_directory("${${_pkg}_BINARY_DIR}/include")
    endif()

    break()
  endif()

endforeach()

if(TARGET "${_lib}::${_lib}")
  set("${_pkg}_FOUND" TRUE)
else()
  set("${_pkg}_FOUND" FALSE)
endif()

if(DEFINED "_${_pkg}_SUBPROJECT")
  unset("${_pkg}_SUBPROJECT")
  unset("_${_pkg}_SUBPROJECT")
endif()
if(DEFINED "_${_pkg}_DIR")
  unset("${_lib}_DIR")
  unset("_${_pkg}_DIR")
endif()
if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
```
</details>


## License

This project is free software: you can redistribute it and/or modify it under the terms of the [Apache License, Version 2.0](LICENSE-Apache) or [MIT license](LICENSE-MIT) at your opinion.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an _as is_ basis, without warranties or conditions of any kind, either express or implied. See the License for the specific language governing permissions and limitations under the License.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this project by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
