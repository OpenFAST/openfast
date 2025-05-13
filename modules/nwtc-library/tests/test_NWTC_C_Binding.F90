module test_NWTC_C_Binding

    use iso_c_binding
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NWTC_C_Binding
    
    implicit none
    private
    public :: test_NWTC_C_Binding_suite
    
    ! C string for testing
    integer(c_int), parameter :: test_c_string_len = 10 + 1 ! +1 for null terminator
    character(kind=c_char) :: test_c_string(test_c_string_len) = (/ "1", "2", "3", "4", "5", "6", "7", "8", "9", "0", C_NULL_CHAR /)

    contains
    
    !> Collect all exported unit tests
    subroutine test_NWTC_C_Binding_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("test_SetErrStat_F2C", test_SetErrStat_F2C),               &
            new_unittest("test_SetErrStat_C", test_SetErrStat_C),                   &
            new_unittest("test_StringConvert_F2C", test_StringConvert_F2C),         &
            new_unittest("test_StringConvert_C2F", test_StringConvert_C2F),         &
            new_unittest("test_RemoveCStringNullChar", test_RemoveCStringNullChar), &
            new_unittest("test_FileNameFromCString", test_FileNameFromCString)      & !, &
        ]
    end subroutine

    subroutine test_SetErrStat_F2C(error)
        type(error_type), allocatable, intent(out) :: error
    
        ! This case should result in error status 0.
        ! It tests that the Fortran ErrStat and ErrMsg are appropriately converted to C-based types

        ! Fortran-based inputs
        integer :: error_status_f = 11
        character(ErrMsgLen) :: error_message_f = "Error message"

        ! C-based outputs to test against inputs
        integer(c_int) :: error_status_c
        character(kind=c_char) :: error_message_c(ErrMsgLen_C)

        call SetErrStat_F2C(error_status_f, error_message_f, error_status_c, error_message_c)

        call check(error, error_status_f, error_status_c)
        ! call check(error, "Error message", error_message_c)


        !!!

    end subroutine

    subroutine test_SetErrStat_C(error)
        type(error_type), allocatable, intent(out) :: error
    
        ! This case should result in error status 0.
        ! It tests that the local error status and error message are incorporated into
        ! the global error status and error message variables.
        ! All data types are C-based.

        integer(c_int) :: error_status_local
        integer(c_int) :: error_status_global
        character(kind=c_char) :: error_message_local(ErrMsgLen_C)
        character(kind=c_char) :: error_message_global(ErrMsgLen_C)
        character(kind=c_char) :: routine_name(IntfStrLen)

        character(len=ErrMsgLen) :: error_message_f
        integer :: loc_global, loc_local, loc_routine

        error_message_local = ""//C_NULL_CHAR
        error_message_global = ""//C_NULL_CHAR
        CALL StringConvert_F2C("test_SetErrStat_C", routine_name)

        ! Check ErrStat handling

        ! If the local error status is 0, the global variables should not change
        error_status_local = 0
        error_status_global = 1
        call SetErrStat_C(error_status_local, error_message_local, error_status_global, error_message_global, routine_name)
        call check(error, 1, error_status_global)

        ! If the local error status is larger than the global, the global should be updated to the local
        error_status_global = 0
        error_status_local = 1
        call SetErrStat_C(error_status_local, error_message_local, error_status_global, error_message_global, routine_name)
        call check(error, 1, error_status_global)

        ! If the local error status is smaller than the global, the global should not change
        error_status_global = 1
        error_status_local = 0
        call SetErrStat_C(error_status_local, error_message_local, error_status_global, error_message_global, routine_name)
        call check(error, 1, error_status_global)

        ! Check ErrMsg handling

        error_status_local = 1
        error_status_global = 2
        CALL StringConvert_F2C("Global error message", error_message_global)        ! Store a Fortran string into the C string test variable

        ! The local error message should be appended to the global error message

        CALL StringConvert_F2C("Local error message", error_message_local)
        CALL SetErrStat_C(error_status_local, error_message_local, error_status_global, error_message_global, routine_name)
        CALL StringConvert_C2F(error_message_global, error_message_f)               ! Convert the C string back to Fortran for checking

        loc_global = INDEX(error_message_f, "Global error message")
        loc_local = INDEX(error_message_f, "Local error message")
        loc_routine = INDEX(error_message_f, "test_SetErrStat_C")
        call check(error, .TRUE., loc_global > 0)
        call check(error, .TRUE., loc_local > 0)
        call check(error, .TRUE., loc_routine > 0)
        call check(error, 1, loc_global)                        ! The global error message should be first
        call check(error, .TRUE., loc_routine < loc_local)      ! Then the routine name should be appended

    end subroutine

    subroutine test_StringConvert_F2C(error)
        type(error_type), allocatable, intent(out) :: error
    
        ! This case should result in error status 0.
        ! It tests that the Fortran string is converted to a C string correctly.

        character(len=test_c_string_len - 1) :: f_string
        character(kind=c_char) :: c_string(test_c_string_len)
        ! character(len=test_c_string_len + 3) :: f_string_longer
        
        integer :: i
        integer :: index_array(11)

        ! Create an array to test for where the C null character should be in the C string
        ! Note that INDEX of a C string returns the index of the substring for each character in the string
        ! It essentially looks like a binary result for each element of the string
        index_array = 0
        index_array(11) = 1

        f_string = "1234567890"
        call StringConvert_F2C(f_string, c_string)
        do i=1, 11
            call check(error, index_array(i), index(c_string(i), C_NULL_CHAR))  ! Check that the C null character is added
            if (i<11) then
                call check(error, f_string(i:i), c_string(i))                   ! Check that the C string matches the original Fortran string
            end if
        end do

        ! Fortran string is longer than C string
        ! f_string_longer = "12345678901234"
        ! call StringConvert_F2C(f_string_longer, c_string)
        ! Note this case is not handled in the subroutine. It MUST be handled in the calling code.
        ! As is, the subroutine will write past the c string array and there is no error.

    end subroutine

    subroutine test_StringConvert_C2F(error)
        type(error_type), allocatable, intent(out) :: error

        ! This case should result in error status 0.
        ! It tests that the C string is converted to Fortran correctly.

        character(len=test_c_string_len - 1) :: f_string
        character(len=test_c_string_len - 2) :: f_string_shorter

        ! Convert the C string to Fortran
        call StringConvert_C2F(test_c_string, f_string)
        call check(error, 0, index(f_string, C_NULL_CHAR))  ! Check that the C null character is removed
        call check(error, "1234567890", f_string)           ! Check that the Fortran string matches the original C string
        
        ! Verify that the C string is truncated when it is longer than the Fortran string
        call StringConvert_C2F(test_c_string, f_string_shorter)
        call check(error, "123456789", f_string_shorter)           ! Check that the Fortran string matches the original C string

    end subroutine

    subroutine test_RemoveCStringNullChar(error)
        type(error_type), allocatable, intent(out) :: error
    
        ! This case should result in error status 0.
        ! It tests that RemoveCStringNullChar indeed removes the C null character.

        ! Make this longer than test_c_string so that we don't accidentally chop off
        ! the null terminator without checking that the function removed it
        integer, parameter :: result_length = test_c_string_len + 10
        character(len=result_length) :: result_string

        result_string = RemoveCStringNullChar(test_c_string, test_c_string_len)
        call check(error, 0, INDEX(result_string, C_NULL_CHAR))

    end subroutine

    subroutine test_FileNameFromCString(error)
        type(error_type), allocatable, intent(out) :: error
    
        ! This case should result in error status 0.
        ! It tests that FileNameFromString extracts the filename from a C string correctly.

        integer, parameter :: input_string_len = IntfStrLen + 8  ! Make the input string slightly longer than IntfStrLen to check that it truncates correctly
        character(kind=c_char, len=input_string_len) :: input_string
        character(len=IntfStrLen) :: file_name

        integer :: i

        ! Fill the test string with 'A's
        do i=1, input_string_len
            input_string(i:i) = "A"
        end do

        ! Assign the character at IntfStrLen to 'B'
        input_string(IntfStrLen:IntfStrLen) = "B"

        ! Terminate with C null character
        input_string(input_string_len:input_string_len) = C_NULL_CHAR

        ! Get the FileName from the C string
        file_name = FileNameFromCString(input_string, IntfStrLen+8)

        ! Check that the extracted filename matches the input up to the IntfStrLen
        call check(error, input_string(1:IntfStrLen), file_name)

        ! Check that the test string is truncated at the IntfStrLen correctly
        call check(error, "B", input_string(IntfStrLen:IntfStrLen))

        ! Check the C_NULL_CHAR is removed
        call check(error, 0, index(file_name, C_NULL_CHAR))

    end subroutine

end module
