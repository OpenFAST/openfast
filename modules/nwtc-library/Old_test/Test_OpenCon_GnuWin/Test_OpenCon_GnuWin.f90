!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2025  National Renewable Energy Laboratory
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
PROGRAM Test_OpenCon_GnuWin


! This program test OpenCon() and WrOvr() in SysGnuWin.

! This was created because the console could not be redirected in v1.06.00c-mlb
! and plus signs showed up on screen when WrOvr() was used.

! Programmer: Marshall Buhl
!             National Wind Technology Center
!             National Renewable Energy Laboratory


USE                                NWTC_Library

IMPLICIT                           NONE


   ! Local declarations.

TYPE(ProgDesc), PARAMETER     :: ProgInfo = ProgDesc( 'NWTC Test_OpenCon_GnuWin', '', '')       ! The name, version, and date of this program.


   ! Initialize the NWTC Library, which will initialize Pi-based constants.

CALL NWTC_Init


   ! Print out program and library name, version, and date.

CALL DispNVD ( ProgInfo )


   ! Test overwriting of output.

CALL WrScr  ( '' )
CALL WrNR   ( 'This line should be overwritten.  What happens if it is longer than the second line?' )
CALL WrOver ( 'This line should overwrite the first.' )
CALL WrOver ( 'This line should overwrite the second.' )


STOP
END

