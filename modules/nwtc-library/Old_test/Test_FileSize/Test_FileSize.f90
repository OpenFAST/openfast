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
PROGRAM Test_FileSize

   ! This program tests the FileSize() function in the library.

USE NWTC_Library

IMPLICIT NONE

INTEGER                         :: Unit      = 11

CALL NWTC_Init ( 'Test_FileSize', 'v1.00.01, 01-Mar-2013' )
CALL WrScr     ( ' Running Test_FileSize (v1.00.01, 01-Mar-2013)' )

OPEN ( Unit, FILE='Test_FileSize.f90', STATUS='OLD' )

PRINT '(/,A,/)', '  Size of Test_FileSize.f90 = '//TRIM( Num2LStr( Int( FileSize( Unit ) ) ) )//' bytes.'

CALL ProgPause

STOP
END
