module nwtc_library_test_tools

use NWTC_IO

implicit none

#ifdef _WIN32
    character(9), parameter :: nullfile="NUL"
    character(11), parameter :: terminal="CON"
#else
    character(9), parameter :: nullfile="/dev/null"
    character(11), parameter :: terminal="/dev/stdout"
#endif

integer, parameter :: stdout=CU

contains

subroutine hide_terminal_output()
    open(unit=stdout, file=trim(nullfile))
end subroutine

subroutine show_terminal_output()
    open(unit=stdout, file=terminal, status="old")
end subroutine

end module
