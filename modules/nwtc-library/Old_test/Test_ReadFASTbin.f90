PROGRAM Test_ReadFASTbin


   ! This program is used to test the TestReadFASTbin_Mod module.


   USE                                 :: NWTC_Library

   IMPLICIT                               NONE

      ! Local declarations.

   INTEGER(IntKi)                      :: ErrLev                     ! The error level returned by ReadFASTbin.
   INTEGER(IntKi)                      :: IO_Unit   = 1              ! The IO unit for the FAST binary file.
   INTEGER(IntKi)                      :: IRec                       ! The record index used for DO loops.

   LOGICAL                             :: Error                      ! A flag returned by Get_Arg().

   CHARACTER(200)                      :: ErrMsg                     ! A possible error message returned by ReadFASTbin.
   CHARACTER( 21)                      :: Fmt                        ! The format specifier for the channel headings.
   CHARACTER(  4)                      :: OutExt    = 'outc'         ! The extension of the converted output file.
   CHARACTER(256)                      :: OutFile                    ! The extension of the converted output file.
   CHARACTER(200)                      :: RootName                   ! The root name of the FAST data files.

   TYPE (FASTdataType)                 :: FASTdata                   ! The derived type for holding FAST output data.



      ! Initialize the NWTC Library.

   CALL NWTC_Init ( 'Test_ReadFASTbin', 'v1.01.00a-mlb, 10-Jan-2003' )


      ! Get the name of the FAST binary data file.

   CALL GET_COMMAND_ARGUMENT( 1, FASTdata%File, STATUS=ErrLev )

   IF ( ErrLev > 0 )  THEN
      CALL ProgAbort ( ' Syntax: Test_ReadFASTbin <input_file>' )
   ENDIF

   CALL GetRoot( FASTdata%File, RootName )

   OutFile = TRIM( RootName )//'.'//OutExt


      ! Read the FAST binary file.

   CALL ReadFASTbin ( IO_Unit, .FALSE., FASTdata, ErrLev, ErrMsg )

   IF ( ErrLev /= 0 )  THEN
      CALL ProgAbort ( TRIM( ErrMsg ) )
   END IF ! ( ErrLev /= 0 )


      ! Open the formatted output file.

   CALL OpenFOutFile  ( IO_Unit, TRIM( OutFile ), ErrLev )

   IF ( ErrLev /= 0 )  THEN
      CALL ProgAbort ( 'Unable to open the formatted output file, "'//TRIM( OutFile )//'".' )
   END IF ! ( ErrLev /= 0 )


      ! Write out the contents of the formatted output file (description, channel names and units, and data).

   CALL WrScr1 ( ' Writing formatted data to '//TRIM( OutFile )//'".' )

   WRITE (IO_Unit,'(/,A,/)')  FASTdata%Descr

   Fmt = '(     (1X,A10))'                                           ! Format specified for channel names and units.

   WRITE (Fmt(2:6),'(I5)')  FASTdata%NumChans+1

   WRITE (IO_Unit,Fmt)  FASTdata%ChanNames
   WRITE (IO_Unit,Fmt)  FASTdata%ChanUnits

   Fmt = '(F10.3,     (ES11.3))'                                  ! Format specified for data.  Time is in fixed-point form.

   WRITE (Fmt(8:12),'(I5)')  FASTdata%NumChans+1

   DO IRec=1,FASTdata%NumRecs
      WRITE (IO_Unit,Fmt)  FASTdata%Data(IRec,:)
   END DO


      ! Close the output file.

   CLOSE ( IO_Unit )


      ! Deallocate the arrays.

   IF ( ALLOCATED( FASTdata%ChanNames ) ) DEALLOCATE( FASTdata%ChanNames )
   IF ( ALLOCATED( FASTdata%ChanUnits ) ) DEALLOCATE( FASTdata%ChanUnits )
   IF ( ALLOCATED( FASTdata%Data      ) ) DEALLOCATE( FASTdata%Data      )


      ! We be done.

   CALL WrScr1 ( ' Test_ReadFASTbin processing complete.' )


   STOP

END
