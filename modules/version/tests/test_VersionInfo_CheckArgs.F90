module test_VersionInfo_CheckArgs

use testdrive, only: new_unittest, unittest_type, error_type, check
use VersionInfo
use versioninfo_test_tools

implicit none

private
public :: test_VersionInfo_CheckArgs_suite

contains

!> Collect all exported unit tests
subroutine test_VersionInfo_CheckArgs_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [new_unittest("test_input_file_user_specified", test_input_file_user_specified), &
                new_unittest("test_input_file_default", test_input_file_default), &
                new_unittest("test_input_file_default_user_specified", test_input_file_default_user_specified), &
                new_unittest("test_restart_flag1", test_restart_flag1), &
                new_unittest("test_restart_flag2", test_restart_flag2), &
                new_unittest("test_restart_flag3", test_restart_flag3), &
                new_unittest("test_second_argument", test_second_argument), &
                new_unittest("test_help1", test_help1), &
                new_unittest("test_help2", test_help2), &
                new_unittest("test_version1", test_version1), &
                new_unittest("test_version2", test_version2), &
                new_unittest("test_no_args_no_default", test_no_args_no_default), &
                new_unittest("test_unsupported_flag", test_unsupported_flag), &
                new_unittest("test_restart_bad_syntax", test_restart_bad_syntax) &
                ]
end subroutine

! PASSING CASES

! ************************************************************************
! To read the input file, a default may be provided and a user-specified
! input file may be used instead.

subroutine test_input_file_user_specified(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe FileName

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(1))
   argument_array = ["first_arg.txt   "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

subroutine test_input_file_default(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = "default.txt"
   allocate (argument_array(0))
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "default.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

subroutine test_input_file_default_user_specified(error)
   type(error_type), allocatable, intent(out) :: error

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = "default.txt"
   allocate (argument_array(1))
   argument_array = ["first_arg.txt   "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! ************************************************************************
! Given a restart flag, the flag should be parsed along with the input
! file and second argument regardless of the position of the flag.
! The first argument is optional in this case.

subroutine test_restart_flag1(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe -Flag FileName Arg2

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(3))
   argument_array = ["-restart        ", &
                     "first_arg.txt   ", &
                     "second_arg      "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "second_arg", second_argument); if (allocated(error)) return
   call check(error, "RESTART", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

subroutine test_restart_flag2(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe FileName -Flag Arg2

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(3))
   argument_array = ["first_arg.txt   ", "-restart        ", "second_arg      "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "second_arg", second_argument); if (allocated(error)) return
   call check(error, "RESTART", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

subroutine test_restart_flag3(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe -restart Arg1

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["-restart        ", "first_arg.txt   "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "first_arg.txt", second_argument); if (allocated(error)) return
   call check(error, "RESTART", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! ************************************************************************
! The second argument should be returned when provided even outside
! of the restart functionality.

subroutine test_second_argument(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe FileName Arg2

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["first_arg.txt   ", "second_arg      "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "second_arg", second_argument); if (allocated(error)) return
   call check(error, "", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! ************************************************************************
! The help flag in any position should show the help prompt and exit
! normally.

subroutine test_help1(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe -Flag FileName

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["-h              ", "first_arg.txt   "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "H", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

subroutine test_help2(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe FileName -Flag

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["first_arg.txt   ", "-h              "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "H", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! ************************************************************************
! The version flag in any position should show the version info and exit
! normally.

subroutine test_version1(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe -v FileName

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["-v              ", "first_arg.txt   "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "V", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

subroutine test_version2(error)
   type(error_type), allocatable, intent(out) :: error

   ! executable.exe FileName -VERSION

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["first_arg.txt   ", "-VERSION        "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 0, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "VERSION", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! FAILING CASES

! ************************************************************************
! No arguments and no default input file should exit with an error

subroutine test_no_args_no_default(error)
   type(error_type), allocatable, intent(out) :: error

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(0))
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "", filename); if (allocated(error)) return
   call check(error, 4, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! ************************************************************************
! An unsupported flag should exit with an error

subroutine test_unsupported_flag(error)
   type(error_type), allocatable, intent(out) :: error

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(2))
   argument_array = ["first_arg.txt   ", "-flag           "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "first_arg.txt", filename); if (allocated(error)) return
   call check(error, 4, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "FLAG", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

! ************************************************************************
! The restart flag requires at least one additional argument

subroutine test_restart_bad_syntax(error)
   type(error_type), allocatable, intent(out) :: error

   character(1024) :: filename, second_argument, flag
   integer(IntKi) :: error_status
   character(16), dimension(:), allocatable :: argument_array

   filename = ""
   allocate (argument_array(1))
   argument_array = ["-restart        "]
   call hide_terminal_output()
   call CheckArgs(filename, error_status, second_argument, flag, argument_array)
   call show_terminal_output()
   call check(error, "", filename); if (allocated(error)) return
   call check(error, 4, error_status); if (allocated(error)) return
   call check(error, "", second_argument); if (allocated(error)) return
   call check(error, "RESTART", flag); if (allocated(error)) return
   deallocate (argument_array)
end subroutine

end module
