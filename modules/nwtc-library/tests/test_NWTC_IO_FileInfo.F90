module test_NWTC_IO_FileInfo

    use pFUnit_mod
    use NWTC_IO
    
    implicit none

contains

@test
subroutine test_initfileinfo()

    ! This case should result in error status 0.
    ! It's a normal initialization of FileInfoType.

    integer, parameter :: num_lines = 5
    integer, parameter :: num_files = 1

    integer(IntKi) :: error_status
    character(ErrMsgLen) :: error_message
    character(1024) :: input_strings(num_lines)
    type(FileInfoType) :: file_info_type
    integer :: i

    input_strings = (/ &
        "line 0",  &
        "line 1",  &
        "line 2",  &
        "line 3",  &
        "line 4"   &
    /)
    call InitFileInfo( input_strings, file_info_type, error_status, error_message )

    @assertEqual(num_lines, file_info_type%NumLines)
    @assertEqual(num_files, file_info_type%NumFiles)
    do i = 1, num_lines
        @assertEqual(i, file_info_type%FileLine(i))
    end do
    do i = 1, num_files
        @assertEqual(i, file_info_type%FileIndx(i))
    end do

    ! TODO: test FileList when we actually use it

    do i = 1, num_lines
        @assertEqual( trim(input_strings(i)), trim(file_info_type%Lines(i)) )
    end do

    @assertEqual( 0, error_status )
    @assertEqual( "", error_message )
end subroutine

@test
subroutine test_initfileinfo2()
    
    ! This case should result in a non-zero error status.
    ! It attempts to initialize FileInfoType without having
    ! properly initializing the input strings array.

    integer, parameter :: num_lines = 5
    integer, parameter :: num_files = 1

    integer(IntKi) :: error_status
    character(ErrMsgLen) :: error_message
    character(1025) :: input_strings(num_lines)
    type(FileInfoType) :: file_info_type
    integer :: i

    input_strings = (/ &
        "line 0",  &
        "line 1",  &
        "line 2",  &
        "line 3",  &
        "line 4"   &
    /)
    call InitFileInfo( input_strings, file_info_type, error_status, error_message )

    @assertEqual(num_lines, file_info_type%NumLines)
    @assertEqual(num_files, file_info_type%NumFiles)
    @assertEqual( 4, error_status )
    
end subroutine

end module