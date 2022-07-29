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
