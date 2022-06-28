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