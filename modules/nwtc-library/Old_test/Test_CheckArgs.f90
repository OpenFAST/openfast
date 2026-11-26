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
   PROGRAM Test_CheckArgs


      ! This program is used to test the use of CheckArgs() in NWTC_IO.f90.


   USE                     :: NWTC_Library

   IMPLICIT                   NONE


   CHARACTER(1024)         :: InpFile = ''

   ProgName = 'Test_CheckArgs'
   ProgVer  = '(v1.00.00, 12-Dec-2012)'


   ! Initialize the NWTC Library, which will initialize Pi-based constants.

CALL NWTC_Init ( ProgName, ProgVer )


   ! Print out program and library name, version, and date.

CALL DispNVD


      ! Get the name of the primary input file from the command line.

   CALL CheckArgs ( InpFile )

   IF ( LEN_TRIM( InpFile ) == 0 )  THEN
      CALL WrScr1    ( ' Syntax is:' )
      CALL WrScr     ( '    '//TRIM( ProgName )//' ['//SwChar//'h] <infile>' )
      CALL WrScr     ( ' where:' )
      CALL WrScr     ( '    '//SwChar//'h generates this help message.' )
!      CALL ProgAbort ( '    <infile> is the name of the required primary input file.', TimeWait=4.0, ErrLevel=1 )
!      CALL ProgAbort ( '    <infile> is the name of the required primary input file.', TimeWait=0.0, ErrLevel=1 )
      CALL ProgAbort ( '    <infile> is the name of the required primary input file.', TimeWait=-1.0, ErrLevel=1 )
   ENDIF


   STOP
   END PROGRAM Test_CheckArgs
