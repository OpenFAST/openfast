module test_NWTC_IO_FileInfo

use iso_c_binding, only: c_char, c_null_char
use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_IO

implicit none
private
public :: test_NWTC_IO_FileInfo_suite

contains

!> Collect all exported unit tests
subroutine test_NWTC_IO_FileInfo_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_initfileinfo", test_initfileinfo), &
               new_unittest("test_initfileinfo2", test_initfileinfo2) &
               ]
end subroutine

subroutine test_initfileinfo(error)
   type(error_type), allocatable, intent(out) :: error

   ! This case should result in error status 0.
   ! It's a normal initialization of FileInfoType.

   integer, parameter :: num_lines = 5
   integer, parameter :: num_files = 1

   integer(IntKi) :: error_status
   character(ErrMsgLen) :: error_message
   character(MaxFileInfoLineLen) :: input_strings(num_lines)
   type(FileInfoType) :: file_info_type
   integer :: i

   input_strings = (/ &
                   "line 0", &
                   "line 1", &
                   "line 2", &
                   "line 3", &
                   "line 4" &
                   /)
   call InitFileInfo(input_strings, file_info_type, error_status, error_message)

   call check(error, num_lines, file_info_type%NumLines)
   call check(error, num_files, file_info_type%NumFiles)
   do i = 1, num_lines
      call check(error, i, file_info_type%FileLine(i))
   end do
   do i = 1, num_files
      call check(error, i, file_info_type%FileIndx(i))
   end do

   ! TODO: test FileList when we actually use it

   do i = 1, num_lines
      call check(error, trim(input_strings(i)), trim(file_info_type%Lines(i)))
   end do

   call check(error, 0, error_status)
   call check(error, "", error_message)
end subroutine

subroutine test_initfileinfo2(error)
   type(error_type), allocatable, intent(out) :: error

   ! This case should result in a non-zero error status.
   ! It attempts to initialize FileInfoType without having
   ! properly initializing the input strings array.

   integer, parameter :: num_lines = 5
   integer, parameter :: num_files = 1

   integer(IntKi) :: error_status
   character(ErrMsgLen) :: error_message
   character(MaxFileInfoLineLen*2) :: input_strings(num_lines)
   character(MaxFileInfoLineLen*2) :: tmpstring
   type(FileInfoType) :: file_info_type
   integer :: i

   input_strings = (/ &
                   "line 0", &
                   "line 1", &
                   "line 2", &
                   "line 3", &
                   "line 4" &
                   /)
   ! make the last character a + so a trim does not reduce this string length
   tmpstring = input_strings(5)
   tmpstring(MaxFileInfoLineLen + 1:MaxFileInfoLineLen + 1) = 'a'
   input_strings(5) = tmpstring
   call InitFileInfo(input_strings, file_info_type, error_status, error_message)
   call check(error, num_lines, file_info_type%NumLines)
   call check(error, num_files, file_info_type%NumFiles)
   call check(error, 4, error_status)

end subroutine

subroutine test_initfileinfo3(error)
   type(error_type), allocatable, intent(out) :: error

   ! This case should result in zero error status.
   ! It attempts to initialize FileInfoType with a C_NULL_CHAR delimited string and compare with the equivalent array parsing

   integer, parameter :: num_lines = 7
   integer, parameter :: num_files = 1

   integer(IntKi) :: error_status
   character(ErrMsgLen) :: error_message
   character(kind=C_CHAR, len=MaxFileInfoLineLen*2) :: input_string
   character(MaxFileInfoLineLen) :: input_string_array(num_lines)
   type(FileInfoType) :: file_info_type
   type(FileInfoType) :: file_info_type_array
   integer :: i

   ! Fortron string array pass
   input_string_array = (/ &
                        "line  0", &
                        "line  1", &
                        "line  2", &
                        "line  3", &
                        "line  4", &
                        "line  5", &
                        "line  6" &
                        /)
   call InitFileInfo(input_string_array, file_info_type_array, error_status, error_message)
   call check(error, num_lines, file_info_type_array%NumLines)
   call check(error, num_files, file_info_type_array%NumFiles)
   call check(error, 0, error_status)

   ! Equivalent C_CHAR string to pass
   ! Note: the rest of the input string is empty.  This won't pose an issue since the remainder of the empty string is ignored.
   input_string = "line  0"//C_NULL_CHAR// &
                  "line  1"//C_NULL_CHAR// &
                  "line  2"//C_NULL_CHAR// &
                  "line  3"//C_NULL_CHAR// &
                  "line  4"//C_NULL_CHAR// &
                  "line  5"//C_NULL_CHAR// &
                  "line  6"//C_NULL_CHAR
   call InitFileInfo(input_string, file_info_type, error_status, error_message)
   call check(error, num_lines, file_info_type%NumLines)
   call check(error, num_files, file_info_type%NumFiles)
   call check(error, 0, error_status)

   ! Now check that the entries are identical
   do i = 1, num_lines
      call check(error, trim(file_info_type_array%Lines(i)), trim(file_info_type%Lines(i)))
   end do

end subroutine

end module
