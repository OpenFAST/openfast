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
PROGRAM Test_ChkRealFmtStr


   ! This program tests NWTC_IO.f90\ChkRealFmtStr().


USE NWTC_Library

IMPLICIT                NONE

INTEGER(IntKi)      :: ErrStat                               ! An error status to be returned by ChkRealFmtStr.
INTEGER(IntKi)      :: FmtWidth                              ! The number of characters that will result from writes using RealFmt.

CHARACTER(200)      :: ErrMsg                                ! An error message to be returned by ChkRealFmtStr.



CALL NWTC_Init ( 'Test_ChkRealFmtStr', 'v1.00.00, 08-Jan-2013' )
CALL WrScr     ( ' Running Test_ChkRealFmtStr (v1.00.00, 07-Jan-2013)' )

CALL ChkRealFmtStr ( '1ES12.4', 'RealFmt', FmtWidth, ErrStat, ErrMsg )
CALL WrScr1 ( ' Message = "'//TRIM(ErrMsg)//'"' )
CALL WrScr  ( ' Format = "1ES12.4", Field width = '//TRIM( Num2LStr( FmtWidth ) ) )

CALL ChkRealFmtStr ( 'F11.4', 'RealFmt', FmtWidth, ErrStat, ErrMsg )
CALL WrScr1 ( ' Message = "'//TRIM(ErrMsg)//'"' )
CALL WrScr  ( ' Format = "F11.4", Field width = '//TRIM( Num2LStr( FmtWidth ) ) )

CALL ChkRealFmtStr ( 'xyzzy', 'RealFmt', FmtWidth, ErrStat, ErrMsg )
CALL WrScr1 ( ' Message = "'//TRIM(ErrMsg)//'"' )
CALL WrScr  ( ' Format = "xyzzy", Field width = '//TRIM( Num2LStr( FmtWidth ) ) )

CALL ChkRealFmtStr ( 'F1.0', 'RealFmt', FmtWidth, ErrStat, ErrMsg )
CALL WrScr1 ( ' Message = "'//TRIM(ErrMsg)//'"' )
CALL WrScr  ( ' Format = "F1.0", Field width = '//TRIM( Num2LStr( FmtWidth ) ) )

CALL ChkRealFmtStr ( '1ES2.4', 'RealFmt', FmtWidth, ErrStat, ErrMsg )
CALL WrScr1 ( ' Message = "'//TRIM(ErrMsg)//'"' )
CALL WrScr  ( ' Format = "1ES2.4", Field width = '//TRIM( Num2LStr( FmtWidth ) ) )

CALL ChkRealFmtStr ( '1ES2.1', 'RealFmt', FmtWidth, ErrStat, ErrMsg )
CALL WrScr1 ( ' Message = "'//TRIM(ErrMsg)//'"' )
CALL WrScr  ( ' Format = "1ES2.1", Field width = '//TRIM( Num2LStr( FmtWidth ) ) )

CALL WrScr ( '' )


STOP
END