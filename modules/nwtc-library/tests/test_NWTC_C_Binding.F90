module test_NWTC_C_Binding

    use iso_c_binding
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use NWTC_C_Binding
    
    implicit none
    private
    public :: test_NWTC_C_Binding_suite
    
    contains
    
    !> Collect all exported unit tests
    subroutine test_NWTC_C_Binding_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("test_FileNameFromCString", test_FileNameFromCString) & !, &
            ! new_unittest("test_initfileinfo2", test_initfileinfo2) &
        ]
    end subroutine
    
    subroutine test_FileNameFromCString(error)
        type(error_type), allocatable, intent(out) :: error
    
        ! This case should result in error status 0.
        ! It tests that FileNameFromString extracts the filename from a string correctly.
        
        integer(IntKi) :: error_status
        character(ErrMsgLen) :: error_message
        character(len=IntfStrLen+8) :: input_string     ! Make the input string slightly longer than IntfStrLen to check that it truncates correctly
        character(len=IntfStrLen) :: file_name

        ! Each long line below has 256 characters
        input_string = &
            "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf"&
            "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf"&
            "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf"&
            "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf"&
            "asdfasdf"

        file_name = FileNameFromCString(input_string, IntfStrLen+8)
    
        call check(error, file_name, input_string(1:IntfStrLen))

    end subroutine


end module
