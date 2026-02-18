module versioninfo_test_tools

use NWTC_IO

implicit none

#ifdef _WIN32
    character(9), parameter :: nullfile="NUL"
    character(11), parameter :: terminal="CON"
#else
    character(9), parameter :: nullfile="/dev/null"
    character(11), parameter :: terminal="/dev/stdout"
#endif

contains

subroutine hide_terminal_output()
    open(unit=CU, file=trim(nullfile))
end subroutine

subroutine show_terminal_output()
    open(unit=CU, file=terminal, status="old")
end subroutine
    
end module
