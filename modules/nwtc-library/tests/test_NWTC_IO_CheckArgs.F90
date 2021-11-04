module test_NWTC_IO_CheckArgs

    use pFUnit_mod
    use NWTC_IO
    use nwtc_library_test_tools
    
    implicit none

contains

    ! PASSING CASES

    ! ************************************************************************
    ! To read the input file, a default may be provided and a user-specified
    ! input file may be used instead.

    @test
    subroutine test_input_file_user_specified()

        ! executable.exe FileName

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(1))
        argument_array = (/      &
            "first_arg.txt   "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "", flag )
        deallocate(argument_array)
    end subroutine

    @test
    subroutine test_input_file_default()

        ! executable.exe

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = "default.txt"
        allocate(argument_array(0))
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "default.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "", flag )
        deallocate(argument_array)
    end subroutine

    @test
    subroutine test_input_file_default_user_specified()

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = "default.txt"
        allocate(argument_array(1))
        argument_array = (/      &
            "first_arg.txt   "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "", flag )
        deallocate(argument_array)
    end subroutine

    ! ************************************************************************
    ! Given a restart flag, the flag should be parsed along with the input
    ! file and second argument regardless of the position of the flag.
    ! The first argument is optional in this case.

    @test
    subroutine test_restart_flag1()

        ! executable.exe -Flag FileName Arg2

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(3))
        argument_array = (/      &
            "-restart        ",  &
            "first_arg.txt   ",  &
            "second_arg      "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "second_arg", second_argument )
        @assertEqual( "RESTART", flag )
        deallocate(argument_array)
    end subroutine

    @test
    subroutine test_restart_flag2()

        ! executable.exe FileName -Flag Arg2

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(3))
        argument_array = (/      &
            "first_arg.txt   ",  &
            "-restart        ",  &
            "second_arg      "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "second_arg", second_argument )
        @assertEqual( "RESTART", flag )
        deallocate(argument_array)
    end subroutine

    @test
    subroutine test_restart_flag3()

        ! executable.exe -restart Arg1

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "-restart        ",  &
            "first_arg.txt   "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "first_arg.txt", second_argument )
        @assertEqual( "RESTART", flag )
        deallocate(argument_array)
    end subroutine

    ! ************************************************************************
    ! The second argument should be returned when provided even outside
    ! of the restart functionality.

    @test
    subroutine test_second_argument()

        ! executable.exe FileName Arg2

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "first_arg.txt   ",  &
            "second_arg      "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "second_arg", second_argument )
        @assertEqual( "", flag )
        deallocate(argument_array)
    end subroutine

    ! ************************************************************************
    ! The help flag in any position should show the help prompt and exit
    ! normally.

    @test
    subroutine test_help1()

        ! executable.exe -Flag FileName

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "-h              ",  &
            "first_arg.txt   "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "H", flag )
        deallocate(argument_array)
    end subroutine

    @test
    subroutine test_help2()

        ! executable.exe FileName -Flag 

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "first_arg.txt   ",  &
            "-h              "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "H", flag )
        deallocate(argument_array)
    end subroutine

    ! ************************************************************************
    ! The version flag in any position should show the version info and exit
    ! normally.

    @test
    subroutine test_version1()

        ! executable.exe -v FileName

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "-v              ",  &
            "first_arg.txt   "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "V", flag )
        deallocate(argument_array)
    end subroutine

    @test
    subroutine test_version2()

        ! executable.exe FileName -VERSION 

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "first_arg.txt   ",  &
            "-VERSION        "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 0, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "VERSION", flag )
        deallocate(argument_array)
    end subroutine

    ! FAILING CASES

    ! ************************************************************************
    ! No arguments and no default input file should exit with an error

    @test
    subroutine test_no_args_no_default()

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(0))
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "", filename )
        @assertEqual( 4, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "", flag )
        deallocate(argument_array)
    end subroutine

    ! ************************************************************************
    ! An unsupported flag should exit with an error

    @test
    subroutine test_unsupported_flag()

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(2))
        argument_array = (/      &
            "first_arg.txt   ",  &
            "-flag           "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "first_arg.txt", filename )
        @assertEqual( 4, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "FLAG", flag )
        deallocate(argument_array)
    end subroutine

    ! ************************************************************************
    ! The restart flag requires at least one additional argument

    @test
    subroutine test_restart_bad_syntax()

        character(1024) :: filename, second_argument, flag
        integer(IntKi) :: error_status
        character(16), dimension(:), allocatable :: argument_array

        filename = ""
        allocate(argument_array(1))
        argument_array = (/      &
            "-restart        "   &
        /)
        call hide_terminal_output()
        call CheckArgs( filename, error_status, second_argument, flag, argument_array )
        call show_terminal_output()
        @assertEqual( "", filename )
        @assertEqual( 4, error_status )
        @assertEqual( "", second_argument )
        @assertEqual( "RESTART", flag )
        deallocate(argument_array)
    end subroutine

end module
