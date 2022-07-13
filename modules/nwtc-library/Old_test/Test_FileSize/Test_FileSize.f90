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
