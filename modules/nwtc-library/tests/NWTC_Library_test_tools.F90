!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020-2025  National Renewable Energy Laboratory
! Copyright (C) 2025-2026  National Laboratory of the Rockies
!
!    This file is part of the NWTC Subroutine Library.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!*********************************************************************************************************************************
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

contains

subroutine hide_terminal_output()
    open(unit=CU, file=trim(nullfile))
end subroutine

subroutine show_terminal_output()
    open(unit=CU, file=terminal, status="old")
end subroutine

end module
