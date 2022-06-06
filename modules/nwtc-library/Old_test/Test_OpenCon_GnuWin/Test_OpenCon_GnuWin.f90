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

