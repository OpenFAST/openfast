!**********************************************************************************************************************************
!
!  MODULE: InflowWind_Driver_Subs  - This module contains subroutines used by the InflowWind Driver program
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
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
!
!**********************************************************************************************************************************
! File last committed: $Date: 2014-07-29 13:30:04 -0600 (Tue, 29 Jul 2014) $
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE InflowWind_Driver_Subs

   USE NWTC_Library
   USE InflowWind_Driver_Types
   IMPLICIT NONE


CONTAINS
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
SUBROUTINE DispHelpText( ErrStat, ErrMsg )
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
      ! Print out help information  !
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

   USE NWTC_Library

   IMPLICIT NONE

      ! Error Handling
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

   ErrStat  =   ErrID_None
   ErrMsg   =   ''


      !  Statement about usage
   CALL WrScr("")
   CALL WrScr("  Syntax:  InlowWind_Driver <filename> [options]")
   CALL WrScr("")
   CALL WrScr("       where:     <filename>     -- Name of driver input file to use")
   CALL WrScr("     options:     "//SWChar//"ifw           -- treat <filename> as name of InflowWind input file [N/A]")
   CALL WrScr("")
   CALL WrScr("              The following options will overwrite values in the driver input file:")
   CALL WrScr("                  "//SwChar//"DT[#]         -- timestep                                        [N/A]")
   CALL WrScr("                  "//SwChar//"TStart[#]     -- start time                                      [N/A]")
   CALL WrScr("                  "//SwChar//"TSteps[#]     -- number of timesteps                             [N/A]")
   CALL WrScr("                  "//SwChar//"xrange[#:#]   -- range of x (#'s are reals)                      [N/A]")
   CALL WrScr("                  "//SwChar//"yrange[#:#]   -- range of y                                      [N/A]")
   CALL WrScr("                  "//SwChar//"zrange[#:#]   -- range in z (ground = 0.0)                       [N/A]")
   CALL WrScr("                  "//SwChar//"Dx[#]         -- spacing in x                                    [N/A]")
   CALL WrScr("                  "//SwChar//"Dy[#]         -- spacing in y                                    [N/A]")
   CALL WrScr("                  "//SwChar//"Dz[#]         -- spacing in z                                    [N/A]")
   CALL WrScr("                  "//SwChar//"sum           -- summarize wind file info                        [N/A]")
   CALL WrScr("                  "//SwChar//"FFT[X,Y,Z]    -- an fft over all t at X,Y,Z (outputs .fft file)  [N/A]")
   CALL WrScr("                  "//SwChar//"points[FILE]  -- calculates at x,y,z coordinates specified in a  [N/A]")
   CALL WrScr("                                    white space delimited FILE")
   CALL WrScr("                  "//SwChar//"help          -- print this help menu and exit")
   CALL WrScr("")
   CALL WrScr("   Notes:")
   CALL WrScr("   -- Unspecified ranges and resolutions default to what is in the file.")
   CALL WrScr("   -- Features marked [N/A] have not been implimented in this version.")
   CALL WrScr("   -- Options are not case sensitive.")
   CALL WrScr("   -- If no XRange is specified, assumed to be only at X=0")
   CALL WrScr("")


END SUBROUTINE DispHelpText


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine retrieves the command line arguments and passes them to the
!! ::ParseArgs routine for processing.
SUBROUTINE RetrieveArgs( CLSettings, CLSettingsFlags, ErrStat, ErrMsg )

   USE NWTC_Library
   USE InflowWind_Driver_Types

   IMPLICIT NONE

      ! Storing the arguments
   TYPE( IfWDriver_Flags ),            INTENT(  OUT)  :: CLSettingsFlags      !< Flags indicating which command line arguments were specified
   TYPE( IfWDriver_Settings ),         INTENT(  OUT)  :: CLSettings           !< Command line arguments passed in

      ! Error Handling
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

      ! Local variables
   INTEGER(IntKi)                                     :: i                    !< Generic counter
   CHARACTER(1024)                                    :: Arg                  !< argument given
   CHARACTER(1024)                                    :: ArgUC                !< Upper case argument to check
   INTEGER(IntKi)                                     :: NumInputArgs         !< Number of argements passed in from command line
   LOGICAL                                            :: ifwFlag              !< The -ifw flag was set
   CHARACTER(1024)                                    :: FileName             !< Filename from the command line.
   LOGICAL                                            :: FileNameGiven        !< Flag indicating if a filename was given.

   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status (for calls)
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error message (for calls)


      ! initialize some things
   CLSettingsFlags%DvrIptFile =  .FALSE.
   ErrStat                    =  ErrID_None
   ErrStatTmp                 =  ErrID_None
   ErrMsg                     =  ""
   ErrMsgTmp                  =  ""
   ifwFlag                    =  .FALSE.
   FileNameGiven              =  .FALSE.
   FileName                   =  ""


      ! Check how many arguments are passed in
   NumInputArgs = COMMAND_ARGUMENT_COUNT()

      ! exit if we don't have enough
   IF (NumInputArgs == 0) THEN
      CALL SetErrStat(ErrID_Fatal," Insufficient Arguments. Use option "//SwChar//"help for help menu.",  &
         ErrStat,ErrMsg,'RetrieveArgs')
      RETURN
   ENDIF


      ! Loop through all the arguments, and store them
   DO i=1,NumInputArgs
         ! get the ith argument
      CALL get_command_argument(i, Arg)
      ArgUC =  Arg

         ! convert to uppercase
      CALL Conv2UC( ArgUC )


         ! Check to see if it is a control parameter or the filename
      IF ( INDEX( SwChar, ArgUC(1:1) ) > 0 ) THEN

            ! check to see if we asked for help
         IF ( ArgUC(2:5) == "HELP" ) THEN
            CALL DispHelpText( ErrStat, ErrMsg )
            CALL ProgExit(0)
         ENDIF


            ! Check the argument and put it where it belongs
            ! chop the SwChar off before passing the argument
         CALL ParseArg( CLSettings, CLSettingsFlags, ArgUC(2:), Arg(2:), ifwFlag, ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'RetrieveArgs')
         IF (ErrStat>AbortErrLev) RETURN

      ELSE

            ! since there is no switch character, assume it is the filename, unless we already set one
         IF ( FileNameGiven ) THEN
            CALL SetErrStat(ErrID_Fatal," Multiple driver input filenames given: "//TRIM(FileName)//", "//TRIM(Arg),   &
               ErrStat,ErrMsg,'RetrieveArgs')
            RETURN
         ELSE
            FileName       = TRIM(Arg)
            FileNameGiven  = .TRUE.
         ENDIF

      ENDIF
   END DO


      ! Was a filename given?
   IF ( .NOT. FileNameGiven ) THEN
      CALL SetErrStat( ErrID_Fatal, " No filename given.", ErrStat, ErrMsg, 'RetrieveArgs' )
      RETURN
   ENDIF

      ! Was the -ifw flag set?  If so, the filename is the InflowWind input file.  Otherwise
      ! it is the driver input file.
   IF ( ifwFlag ) THEN
      CLSettings%IfWIptFileName  =  TRIM(FileName)
      CLSettingsFlags%IfWIptFile =  .TRUE.
   ELSE
      CLSettings%DvrIptFileName  =  TRIM(FileName)
      CLSettingsFlags%DvrIptFile =  .TRUE.
   ENDIF



   !-------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
      CONTAINS


   !-------------------------------------------------------------------------------
   FUNCTION StringToReal( StringIn, ErrStat )
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
         ! Convert a string to a real number !
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

      IMPLICIT NONE

         ! Error Handling
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat

         ! Input
      CHARACTER(*),                       INTENT(IN   )  :: StringIn

         ! Returned value
      REAL(ReKi)                                         :: StringToReal

         ! Local Variables
      INTEGER(IntKi)                                     :: ErrStatTmp         ! Temporary variable to hold the error status

         read( StringIn, *, iostat=ErrStatTmp) StringToReal

            ! If that isn't a number, only warn since we can continue by skipping this value
         IF ( ErrStatTmp .ne. 0 ) ErrStat  = ErrID_Warn

   END FUNCTION StringToReal



   !-------------------------------------------------------------------------------
   SUBROUTINE ParseArg( CLSettings, CLSettingsFlags, ThisArgUC, ThisArg, ifwFlagSet, ErrStat, ErrMsg )
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
         ! Parse and store the input argument  !
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

      USE NWTC_Library
      USE InflowWind_Driver_Types
      USE InflowWind_Types

      IMPLICIT NONE

         ! Storing the arguments
      TYPE( IfWDriver_Flags ),            INTENT(INOUT)  :: CLSettingsFlags      ! Flags indicating which arguments were specified
      TYPE( IfWDriver_Settings ),         INTENT(INOUT)  :: CLSettings           ! Arguments passed in

      CHARACTER(*),                       INTENT(IN   )  :: ThisArgUC            ! The current argument (upper case for testing)
      CHARACTER(*),                       INTENT(IN   )  :: ThisArg              ! The current argument (as passed in for error messages)
      LOGICAL,                            INTENT(INOUT)  :: ifwFlagSet           ! Was the -ifw flag given?

         ! Error Handling
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


         ! local variables
      INTEGER(IntKi)                                     :: Delim1               ! where the [ is
      INTEGER(IntKi)                                     :: Delim2               ! where the ] is
      INTEGER(IntKi)                                     :: DelimSep             ! where the : is
      INTEGER(IntKi)                                     :: DelimSep2            ! where the : is
      INTEGER(IntKi)                                     :: DelimSep3            ! where the : is
      REAL(ReKi)                                         :: TempReal             ! temp variable to hold a real

      INTEGER(IntKi)                                     :: ErrStatTmp           ! Temporary error status for calls



         ! Initialize some things
      ErrStat     = ErrID_None
      ErrStatTmp  =  ErrID_None
      ErrMsg      = ""

         ! Get the delimiters -- returns 0 if there isn't one
      Delim1   = INDEX(ThisArgUC,'[')
      Delim2   = INDEX(ThisArgUC,']')
      DelimSep = INDEX(ThisArgUC,':')


         ! check that if there is an opening bracket, then there is a closing one
      IF ( (Delim1 > 0_IntKi ) .and. (Delim2 < Delim1) ) THEN
         CALL SetErrStat(ErrID_Warn," Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
            ErrStat,ErrMsg,'ParseArg')
         RETURN
      ENDIF

         ! check that if there is a colon, then there are brackets
      IF ( (DelimSep > 0_IntKi) .and. (Delim1 == 0_IntKi) ) THEN
         CALL SetErrStat(ErrID_Warn," Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
            ErrStat,ErrMsg,'ParseArg')
         RETURN
      ENDIF


         ! If no delimeters were given, than this option is simply a flag
      IF ( Delim1 == 0_IntKi ) THEN
            ! check to see if the filename is the name of the IfW input file
         IF       ( ThisArgUC(1:3) == "IFW" )   THEN
            ifwFlagSet              = .TRUE.             ! More logic in the routine that calls this one to set things.
            RETURN
         ELSEIF   ( ThisArgUC(1:3) == "SUM" )   THEN
            CLSettingsFlags%Summary = .TRUE.
            RETURN
         ELSE
            CALL SetErrStat( ErrID_Warn," Unrecognized option '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options.",  &
               ErrStat,ErrMsg,'ParseArg')
         ENDIF

      ENDIF



          ! "FFT[X,Y,Z]"
      IF     ( ThisArgUC(1:Delim1) == "FFT["  ) THEN
         DelimSep = INDEX(ThisArgUC,',')
         DelimSep2= INDEX(ThisArgUC(DelimSep+1:),',') + DelimSep
         IF ( DelimSep2 <= DelimSep ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,'ParseArg')
            RETURN
         ENDIF
            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%FFTcalc    = .TRUE.
            CLSettings%FFTcoord(1)      = TempReal
         ELSE
            CLSettingsFlags%FFTcalc    = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, 'ParseArg')
            ELSE
               CALL SetErrStat(ErrID_FATAL," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:DelimSep2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%FFTcalc    = .TRUE.
            CLSettings%FFTcoord(2)      = TempReal
         ELSE
            CLSettingsFlags%FFTcalc    = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, 'ParseArg')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Third Value
         TempReal = StringToReal( ThisArgUC(DelimSep2+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%FFTcalc    = .TRUE.
            CLSettings%FFTcoord(3)      = TempReal
         ELSE
            CLSettingsFlags%FFTcalc    = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

          ! "XRANGE[#:#]"
      ELSEIF ( ThisArgUC(1:Delim1) == "XRANGE["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF     ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%XRange = .TRUE.
            CLSettings%XRange(1)   = TempReal
         ELSE
            CLSettingsFlags%XRange = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%XRange = .TRUE.
            CLSettings%XRange(2)   = TempReal
         ELSE
            CLSettingsFlags%XRange = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Check the order of values
         IF ( CLSettings%XRange(1) > CLSettings%Xrange(2) ) THEN
            CLSettings%XRange(1)   = 0.0
            CLSettings%XRange(2)   = 0.0
            CLSettingsFlags%XRange = .FALSE.
            CALL SetErrStat(ErrID_Warn," Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
               ErrStat, ErrMsg, 'ParseArgs')
         ENDIF



         ! "YRANGE[#:#]"
      ELSEIF ( ThisArgUC(1:Delim1) == "YRANGE["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF     ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%YRange = .TRUE.
            CLSettings%YRange(1)   = TempReal
         ELSE
            CLSettingsFlags%YRange = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%YRange = .TRUE.
            CLSettings%YRange(2)   = TempReal
         ELSE
            CLSettingsFlags%YRange = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Check the order of values
         IF ( CLSettings%YRange(1) > CLSettings%Yrange(2) ) THEN
            CLSettings%YRange(1)   = 0.0
            CLSettings%YRange(2)   = 0.0
            CLSettingsFlags%YRange = .FALSE.
            CALL SetErrStat(ErrID_Warn," Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring.",   &
               ErrStat, ErrMsg, 'ParseArgs')
         ENDIF



         ! "ZRANGE[#:#]"
      ELSEIF ( ThisArgUC(1:Delim1) == "ZRANGE["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF     ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%ZRange = .TRUE.
            CLSettings%ZRange(1)   = TempReal
         ELSE
            CLSettingsFlags%ZRange = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLSettingsFlags%ZRange = .TRUE.
            CLSettings%ZRange(2)   = TempReal
         ELSE
            CLSettingsFlags%ZRange = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF

            ! Check the order of values
         IF ( CLSettings%ZRange(1) > CLSettings%Zrange(2) ) THEN
            CLSettings%ZRange(1)   = 0.0
            CLSettings%ZRange(2)   = 0.0
            CLSettingsFlags%ZRange = .FALSE.
            CALL SetErrStat(ErrID_Warn," Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring.",   &
               ErrStat, ErrMsg, 'ParseArgs')
         ENDIF


         ! "DX[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "DX["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLSettingsFlags%Dx      = .TRUE.
            CLSettings%GridDelta(1) = abs(TempReal)
         ELSE
            CLSettingsFlags%Dx      = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


         ! "DY[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "DY["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLSettingsFlags%Dy      = .TRUE.
            CLSettings%GridDelta(2) = abs(TempReal)
         ELSE
            CLSettingsFlags%Dy      = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


         ! "DZ[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "DZ["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLSettingsFlags%Dz      = .TRUE.
            CLSettings%GridDelta(3) = abs(TempReal)
         ELSE
            CLSettingsFlags%Dz      = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


         ! "DT[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "DT["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLSettingsFlags%Dt      = .TRUE.
            CLSettings%DT           = abs(TempReal)
         ELSE
            CLSettingsFlags%Dt      = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


         ! "TSTEPS[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "TSTEPS["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLSettingsFlags%NumTimeSteps  = .TRUE.
            CLSettings%NumTimeSteps       = nint(TempReal)
         ELSE
            CLSettingsFlags%NumTimeSteps  = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF



         ! "TSTART[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "TSTART["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLSettingsFlags%TStart  = .TRUE.
            CLSettings%TStart       = abs(TempReal)
         ELSE
            CLSettingsFlags%TStart  = .FALSE.
            IF ( ErrStat == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
                  ErrStat,ErrMsg,'ParseArgs')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF



         ! "POINTS[FILE]"
      ELSEIF( ThisArgUC(1:Delim1) == "POINTS["    ) THEN
         CLSettingsFlags%PointsFile= .TRUE.
         CLSettings%PointsFileName = ThisArgUC(Delim1+1:Delim2-1)
      ELSE
         ErrMsg  = " Unrecognized option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options."
         ErrStat = ErrID_Warn
      ENDIF

   END SUBROUTINE ParseArg
   !-------------------------------------------------------------------------------



END SUBROUTINE RetrieveArgs


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine reads the driver input file and sets up the flags and settings
!! for the driver code.  Any settings from the command line options will override
!! this.
SUBROUTINE ReadDvrIptFile( DvrFileName, DvrFlags, DvrSettings, ProgInfo, ErrStat, ErrMsg )

   CHARACTER(1024),                    INTENT(IN   )  :: DvrFileName
   TYPE(IfWDriver_Flags),              INTENT(INOUT)  :: DvrFlags
   TYPE(IfWDriver_Settings),           INTENT(INOUT)  :: DvrSettings
   TYPE(ProgDesc),                     INTENT(IN   )  :: ProgInfo
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)                                     :: UnIn                 ! Unit number for the driver input file
   CHARACTER(1024)                                    :: FileName             ! Name of InflowWind driver input file

      ! Input file echoing
   LOGICAL                                            :: EchoFileContents     ! Do we echo the driver file out or not?
   INTEGER(IntKi)                                     :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                    :: EchoFileName         ! Name of InflowWind driver echo file

      ! Time steps
   CHARACTER(1024)                                    :: NumTimeStepsChr      ! Character string for number of timesteps (to handle DEFAULT value)
   CHARACTER(1024)                                    :: DTChr                ! Character string for timesteps size (to handle DEFAULT value)

      ! Gridded data
   INTEGER(IntKi)                                     :: TmpIntAr2(2)         ! Temporary array for reading in a pair of integer values from the input file
   REAL(ReKi)                                         :: TmpRealAr2(2)        ! Temporary array for reading in a pair of real values from the input file
   REAL(ReKi)                                         :: GridCtrCoord(3)      ! Center coordinate of the grid read in


      ! Local error handling
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error messages for calls


      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1

   FileName = TRIM(DvrFileName)

   CALL GetNewUnit( UnIn )
   CALL OpenFInpFile( UnIn, FileName, ErrStatTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,' Failed to open InflowWind Driver input file: '//FileName,   &
         ErrStat,ErrMsg,'ReadDvrIptFile')
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL WrScr( 'Opening InflowWind Driver input file:  '//FileName )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName,' InflowWind Driver input file header line 1', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL ReadCom( UnIn, FileName, 'InflowWind Driver input file header line 2', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CLOSE( UnIn )
      RETURN
   ENDIF


     ! Echo Input Files.
   CALL ReadVar ( UnIn, FileName, EchoFileContents, 'Echo', 'Echo Input', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.

   IF ( EchoFileContents ) THEN

      EchoFileName = TRIM(FileName)//'.ech'
      CALL GetNewUnit( UnEchoLocal )
      CALL OpenEcho ( UnEchoLocal, EchoFileName, ErrStatTmp, ErrMsgTmp, ProgInfo )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CLOSE( UnIn )
         RETURN
      ENDIF


      REWIND(UnIn)


         ! Reread and echo
      CALL ReadCom( UnIn, FileName,' InflowWind Driver input file header line 1', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF


      CALL ReadCom( UnIn, FileName, 'InflowWind Driver input file header line 2', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF


        ! Echo Input Files.
      CALL ReadVar ( UnIn, FileName, EchoFileContents, 'Echo', 'Echo Input', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF



   ENDIF


   !-------------------------------------------------------------------------------------------------
   !  Driver setup section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' Driver setup section, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! InflowWind input file
   CALL ReadVar( UnIn, FileName,DvrSettings%IfWIptFileName,'IfWIptFileName',' InflowWind input filename',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%IfWIptFile  =  .TRUE.
   ENDIF


      ! Number of timesteps
   CALL ReadVar( UnIn, FileName,NumTimeStepsChr,'NumTimeStepsChr',' Character string for number of timesteps to read.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( NumTimeStepsChr )
   IF ( TRIM(NumTimeStepsChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%NumTimeSteps         =  .TRUE.
      DvrFlags%NumTimeStepsDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.         
      READ (NumTimeStepsChr,*,IOSTAT=ErrStatTmp)   DvrSettings%NumTimeSteps
      IF ( ErrStatTmp /= ErrID_None )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( ErrStatTmp, "", 'NumTimeSteps',NumType, .TRUE., ErrMsgTmp )
         RETURN
      ELSE     ! Was ok, so set the flags
         DvrFlags%NumTimeSteps         =  .TRUE.
         DvrFlags%NumTimeStepsDefault  =  .FALSE.
      ENDIF
   ENDIF
   

      ! TStart    -- start time
   CALL ReadVar( UnIn, FileName,DvrSettings%TStart,'TStart',' Time in wind file to start parsing.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%TStart   =  .TRUE.
   ENDIF


      ! Summarize the extents in the windfile
   CALL ReadVar( UnIn, FileName,DvrFlags%Summary,'Summary',' Summarize data extents in the windfile', &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%Summary  =  .TRUE.
   ENDIF


      ! Summarize everything in a summary file/
   CALL ReadVar( UnIn, FileName,DvrFlags%SummaryFile,'SummaryFile',' Summarize the results in a .sum file', &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%SummaryFile =  .TRUE.
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !  InflowWind setup section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' InflowWind setup section, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! DT    -- Timestep size for the driver to take (or DEFAULT for what the file contains)
   CALL ReadVar( UnIn, FileName,DTChr,'DTChr',' Character string for Timestep size for the driver to take (or DEFAULT for what the file contains).',  &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( DTChr )
   IF ( TRIM(DTChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%DT         =  .TRUE.
      DvrFlags%DTDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.         
      READ (DTChr,*,IOSTAT=ErrStatTmp)   DvrSettings%DT
      IF ( ErrStatTmp /= ErrID_None )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( ErrStatTmp, "", 'DT',NumType, .TRUE., ErrMsgTmp )
         RETURN
      ELSE     ! Was ok, so set the flags
         DvrFlags%DT         =  .TRUE.
         DvrFlags%DTDefault  =  .FALSE.
      ENDIF
   ENDIF
 

      ! TurbineHeight    -- Height of the hub of the turbine.  Used for error checking (m).
   CALL ReadVar( UnIn, FileName,DvrSettings%TurbineHeight,'TurbineHeight',' Height of the turbine hub (for checking only).',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%TurbineHeight   =  .TRUE.
   ENDIF



       ! Width    -- Width of the windfield needed.
   CALL ReadVar( UnIn, FileName,DvrSettings%Width,'Width',' Width of windfield needed (estimate).',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%Width   =  .TRUE.
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !  FFT calculations
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' FFT calculations, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


       ! FFTcalc    -- FFTcalc of the windfield needed.
   CALL ReadVar( UnIn, FileName,DvrFlags%FFTcalc,'FFTcalc',' Perform an FFT?',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! Read the coordinate for the FFT if the flag is set, otherwise skip the line
   IF ( DvrFlags%FFTcalc   ) THEN
         ! FFTcoord     -- The coordinates to perform the FFT at
      CALL ReadAry ( UnIn, FileName, DvrSettings%FFTcoord, 3, 'FFTcoord', &
         'FFT coordinate', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   ELSE
      CALL ReadCom( UnIn, FileName,' Skipping the FFT coordinate since not doint an FFT.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   ENDIF
 
      

   !-------------------------------------------------------------------------------------------------
   !  gridded data output
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' Gridded data output, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


       ! WindGrid    -- Gridded data output
   CALL ReadVar( UnIn, FileName,DvrFlags%WindGrid,'WindGrid',' Output a grid of data?',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! Read the coordinate for the FFT if the flag is set, otherwise skip the line
   IF ( DvrFlags%WindGrid   ) THEN

         ! GridCtrCoord     -- The coordinates to center the gridded data at
      CALL ReadAry ( UnIn, FileName, GridCtrCoord, 3, 'GridCtrCoord', &
         'Coordinate of the center of the gridded data', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF

         ! Read the DY and DZ stepsize
      CALL ReadAry ( UnIn, FileName, TmpRealAr2, 2, 'GridDY, GridDZ', &
         'GridDY, GridDZ', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF

         ! Save the DY and DZ values
      DvrSettings%GridDelta(2)   =  TmpRealAr2(1)     ! Y direction
      DvrSettings%GridDelta(3)   =  TmpRealAr2(2)     ! Z direction
      DvrFlags%Dx                =  .FALSE.           ! no x direction gridding
      DvrFlags%Dy                =  .TRUE.            ! read in value for the Y direction gridding
      DvrFlags%Dz                =  .TRUE.            ! read in value for the Z direction gridding


         ! Read the number of points in the Y and Z directions
      CALL ReadAry ( UnIn, FileName, TmpIntAr2, 2, 'GridNY, GridNZ', &
         'GridNY, GridNZ', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF

         ! Save the GridNY and GridNZ values
      DvrSettings%GridN(2)   =  TmpIntAr2(1)          ! Y direction
      DvrSettings%GridN(3)   =  TmpIntAr2(2)          ! Z direction

   ELSE
         ! Skip the next three entries of the gridded data section.
      CALL ReadCom( UnIn, FileName,' Skipping the gridded data section since not calculating it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the gridded data section since not calculating it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the gridded data section since not calculating it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadDvrIptFile')
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   ENDIF
 

      ! Now need to set the XRange, YRange, and ZRange values based on what we read in
      ! For the XRange, we don't have any depth (not allowed in the input file)
   DvrSettings%XRange      =  GridCtrCoord(1)

      ! For YRange, check that we have an actual value for the number of points
   IF ( DvrSettings%GridN(2) <= 0 )  THEN
      DvrSettings%YRange   =  abs(GridCtrCoord(2))       ! shouldn't have a negative value anyhow
      DvrFlags%Dy          =  .FALSE.
      DvrSettings%GridN(2) =  0_IntKi

      IF ( DvrSettings%GridN(2) < 0 )  THEN
         CALL SetErrStat(ErrID_Warn,' Negative number for number of grid points along Y direction.  Ignoring.',ErrStat,ErrMsg, 'ReadDvrIptFile')
      ELSE
         CALL SetErrStat(ErrID_Warn,' No points along Y direction.  Ignoring.',ErrStat,ErrMsg, 'ReadDvrIptFile')
      ENDIF

   ELSE
         ! Set the YRange values
      DvrSettings%YRange(1)   =  GridCtrCoord(2) - (REAL(DvrSettings%GridN(2) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(2)
      DvrSettings%YRange(2)   =  GridCtrCoord(2) + (REAL(DvrSettings%GridN(2) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(2)
      DvrFlags%YRange         =  .TRUE.
   ENDIF

      ! For ZRange, check that we have an actual value for the number of points
   IF ( DvrSettings%GridN(3) <= 0 )  THEN
      DvrSettings%ZRange   =  abs(GridCtrCoord(3))       ! shouldn't have a negative value anyhow
      DvrFlags%Dz          =  .FALSE.
      DvrSettings%GridN(3) =  0_IntKi

      IF ( DvrSettings%GridN(3) < 0 )  THEN
         CALL SetErrStat(ErrID_Warn,' Negative number for number of grid points along Z direction.  Ignoring.',ErrStat,ErrMsg, 'ReadDvrIptFile')
      ELSE
         CALL SetErrStat(ErrID_Warn,' No points along Z direction.  Ignoring.',ErrStat,ErrMsg, 'ReadDvrIptFile')
      ENDIF

   ELSE
         ! Set the ZRange values
      DvrSettings%ZRange(1)   =  GridCtrCoord(3) - (REAL(DvrSettings%GridN(3) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(3)
      DvrSettings%ZRange(2)   =  GridCtrCoord(3) + (REAL(DvrSettings%GridN(3) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(3)
      DvrFlags%ZRange         =  .TRUE.
   ENDIF



      ! Close the echo and input file
   CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
   CLOSE( UnIn )


CONTAINS

   !----------------------------------------------------------------------------------------------------
   !> The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to
   !! any existing echo information
   SUBROUTINE CleanupEchoFile( EchoFlag, UnEcho)
      LOGICAL,                      INTENT(IN   )  :: EchoFlag          ! local version of echo flag
      INTEGER(IntKi),               INTENT(IN   )  :: UnEcho            !  echo unit number

         ! Close this module's echo file
      IF ( EchoFlag ) THEN
         CLOSE(UnEcho)
      ENDIF
   END SUBROUTINE CleanupEchoFile



END SUBROUTINE ReadDvrIptFile


!> This subroutine copies an command line (CL) settings over to the program settings.  Warnings are
!! issued if anything is changed from what the driver input file requested.
SUBROUTINE UpdateSettingsWithCL( DvrFlags, DvrSettings, CLFlags, CLSettings, DVRIPT, ErrStat, ErrMsg )

   TYPE(IfWDriver_Flags),              INTENT(INOUT)  :: DvrFlags
   TYPE(IfWDriver_Settings),           INTENT(INOUT)  :: DvrSettings
   TYPE(IfWDriver_Flags),              INTENT(IN   )  :: CLFlags
   TYPE(IfWDriver_Settings),           INTENT(IN   )  :: CLSettings
   LOGICAL,                            INTENT(IN   )  :: DVRIPT
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


      ! Local variables
   INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp         !< Temporary error status for calls
   LOGICAL                                            :: WindGridModify    !< Did we modify any of the WindGrid related settings?

      ! Initialization
   WindGridModify =  .FALSE.

      ! Initialize the error handling
   ErrStat     =  ErrID_None
   ErrMsg      =  ""
   ErrStatTmp  =  ErrID_None
   ErrMsgTmp   =  ""


      ! Due to the complexity, we are handling overwriting driver input file settings with
      ! command line settings and the instance where no driver input file is read separately.
   IF ( DVRIPT .AND. DvrFlags%WindGrid ) THEN

         !--------------------------------------------
         !  Did we change any time information?
         !--------------------------------------------

         !  Check TStart
      IF ( CLFlags%TStart ) THEN
         IF ( DvrFlags%TStart .AND. ( .NOT. EqualRealNos(DvrSettings%TStart, CLSettings%TStart) ) )   THEN
            CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for TStart with '//TRIM(Num2LStr(CLSettings%TStart))//'.', &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ELSE
            DvrFlags%TStart   =  .TRUE.
         ENDIF
         DvrSettings%TStart   =  CLSettings%TStart
      ENDIF

         ! Check DT
      IF ( CLFlags%DT ) THEN
         IF ( DvrFlags%DT .AND. ( .NOT. EqualRealNos(DvrSettings%DT, CLSettings%DT) ) )   THEN
            CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for DT with '//TRIM(Num2LStr(CLSettings%DT))//'.',  &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ELSE
            DvrFlags%DT   =  .TRUE.
         ENDIF
         DvrSettings%DT   =  CLSettings%DT
      ENDIF

         ! Check NumTimeSteps
      IF ( CLFlags%NumTimeSteps ) THEN
         IF ( DvrFlags%NumTimeSteps .AND. ( DvrSettings%NumTimeSteps /= CLSettings%NumTimeSteps ) )   THEN
            CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for NumTimeSteps with '//  &
               TRIM(Num2LStr(CLSettings%NumTimeSteps))//'.',&
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ELSE
            DvrFlags%NumTimeSteps   =  .TRUE.
         ENDIF
         DvrSettings%NumTimeSteps   =  CLSettings%NumTimeSteps
      ENDIF


         !--------------------------------------------
         ! Did we change the ranges for the WindGrid?
         !--------------------------------------------

         ! XRange -- The xrange calculations are a little different than the other two ranges.  The Xrange can only be specified
         !           from the command line.  In order for it to work, we need both an XRange and Dx.  Otherwise we will ignore
         !           this setting.
         !
         !     XRange      Dx       What to change
         !        y        n        Nothing.  Cannot calculate since no Dx or Nx values from driver input
         !        n        y        Nothing.  Cannot calculate since no XRange or Nx values from driver input
         !        y        y        Calculate Nx from XRange and Dx, reset XRange(2) after calculating (extend
         !                          a bit if necessary).  Also set WindGrid to true if it is false.
         !        n        n        Nothing.
         !

      IF       ( CLFlags%XRange .AND. ( .NOT. CLFlags%Dx ) )   THEN

            ! Only the XRange was specfied, so throw warning and don't change anything.
         CALL SetErrStat( ErrID_Warn, ' XRange argument given, but no Dx.  Ignoring since no points to calculate.',  &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
         DvrFlags%XRange   =  .FALSE.     ! Should already be false, but just in case
         DvrFlags%Dx       =  .FALSE.     ! Should already be false, but just in case

      ELSEIF   ( ( .NOT. CLFlags%XRange ) .AND. CLFlags%Dx )   THEN

            ! Only Dx was specfied, so throw warning and don't change anything.
         CALL SetErrStat( ErrID_Warn, ' Dx argument given, but no XRange.  Ignoring since no points to calculate.',  &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
         DvrFlags%XRange   =  .FALSE.     ! Should already be false, but just in case
         DvrFlags%Dx       =  .FALSE.     ! Should already be false, but just in case

      ELSEIF   ( CLFlags%XRange .AND. CLFlags%Dx ) THEN

         DvrFlags%Dx       =  .TRUE.
         DvrFlags%XRange   =  .TRUE.
         WindGridModify    =  .TRUE.

            ! Copy over the range and stepsize
         DvrSettings%XRange(1)      =  minval(CLSettings%XRange)     ! just in case the order is funky
         DvrSettings%XRange(2)      =  maxval(CLSettings%XRange)
         DvrSettings%GridDelta(1)   =  abs(CLSettings%GridDelta(1))  ! just in case there was a negative number

            ! Set number of points in X direction
         DvrSettings%GridN(1) =  CEILING( abs( ( DvrSettings%XRange(2) - DvrSettings%XRange(1) ) / DvrSettings%GridDelta(1) ) + 1_IntKi )

            ! Now adjust the upper limit of the range slightly for the integer number of steps
         DvrSettings%XRange(2)=  DvrSettings%XRange(1)   + ( DvrSettings%GridN(1) - 1_IntKi ) * DvrSettings%GridDelta(1)

            ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
         IF ( .NOT. EqualRealNos(DvrSettings%XRange(2),maxval(CLSettings%XRange)) )  THEN
            CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of XRange to '//TRIM(Num2LStr(DvrSettings%XRange(2)))//   &
               ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(1)))//' points with a requested spacing '//       &
               TRIM(Num2LStr(DvrSettings%GridDelta(1)))//'.',                                                           &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ENDIF

      ELSE
         ! Nothing was specified, so we don't do anything
      ENDIF



         ! YRange -- There are three bits of information required for the full specification of the Y
         !           data: Dy, YRange, and Ny (the number of steps).  In the driver input file, Dy and Ny are
         !           given, and YRange is calculated using the specified center point.  From the command line,
         !           Dy, and YRange are specified.
         !           In order for the command line values to be used, we must know which were specified or not
         !           as it will change the values we use for calculating.  
         !
         !  First, check if WindGrid is set.  If not, set it.
         !
         !     YRange      Dy       What to change
         !        y        n        Assume Dy good.  Calculate new Ny.  Increase YRange(2) if necessary and warn.
         !        n        y        Assume YRange good (if Ny == 0, issue warning about no points, set to none)
         !        y        y        Calculate new Ny.  Increase Yrange(2) if necessary and warn.
         !        n        n        Nothing to do. leave alone.
         !

      IF       ( CLFlags%YRange .AND. ( .NOT. CLFlags%Dy ) )   THEN

         IF ( .NOT. DvrFlags%Dy )         DvrFlags%Dy       =  .TRUE.
         IF ( .NOT. DvrFlags%YRange )     DvrFlags%YRange   =  .TRUE.
         IF ( .NOT. WindGridModify )      WindGridModify    =  .TRUE.

            ! Calculate new Ny.
         DvrSettings%YRange(1)      =  minval(CLSettings%YRange)     ! just in case the order is funky
         DvrSettings%YRange(2)      =  maxval(CLSettings%YRange)
         
            ! Set number of points in Y direction
         DvrSettings%GridN(2) =  CEILING( abs( ( DvrSettings%YRange(2) - DvrSettings%YRange(1) ) / DvrSettings%GridDelta(2) ) + 1_IntKi )

            ! Now adjust the upper limit of the range slightly for the integer number of steps
         DvrSettings%YRange(2)=  DvrSettings%YRange(1)   + ( DvrSettings%GridN(2) - 1_IntKi ) * DvrSettings%GridDelta(2)

            ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
         IF ( .NOT. EqualRealNos(DvrSettings%YRange(2),maxval(CLSettings%YRange)) )  THEN
            CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of YRange to '//TRIM(Num2LStr(DvrSettings%YRange(2)))//   &
               ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(2)))//' points with a requested spacing '//       &
               TRIM(Num2LStr(DvrSettings%GridDelta(2)))//'.',                                                           &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ENDIF

      ELSEIF   ( ( .NOT. CLFlags%YRange ) .AND. CLFlags%Dy )   THEN

            ! Make sure we have points to calculate
         IF ( DvrSettings%GridN(2)  == 0 )   THEN
               ! We only issue this warning.
            CALL SetErrStat( ErrID_Warn,' No points specified in driver input file for Y direction.  Ignoring Dy input.',  &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
            DvrFlags%Dy       =  .FALSE.
            DvrFlags%YRange   =  .FALSE.
         ELSE

            IF ( .NOT. DvrFlags%Dy )         DvrFlags%Dy       =  .TRUE.
            IF ( .NOT. DvrFlags%YRange )     DvrFlags%YRange   =  .TRUE.
            IF ( .NOT. WindGridModify )      WindGridModify    =  .TRUE.

               ! Save the Dy value
            DvrSettings%GridDelta(2)   =  abs(CLSettings%GridDelta(2))  ! just in case there was a negative number

               ! Calculate the new value of Ny
            DvrSettings%GridN(2) =  CEILING( abs( ( DvrSettings%YRange(2) - DvrSettings%YRange(1) ) / DvrSettings%GridDelta(2) ) + 1_IntKi )

               ! Now adjust the upper limit of the range slightly for the integer number of steps
            DvrSettings%YRange(2)=  DvrSettings%YRange(1)   + ( DvrSettings%GridN(2) - 1_IntKi ) * DvrSettings%GridDelta(2)

               ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
            IF ( .NOT. EqualRealNos(DvrSettings%YRange(2),maxval(CLSettings%YRange)) )  THEN
               CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of YRange to '//TRIM(Num2LStr(DvrSettings%YRange(2)))//   &
                  ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(2)))//' points with a requested spacing '//       &
                  TRIM(Num2LStr(DvrSettings%GridDelta(2)))//'.',                                                           &
                  ErrStat,ErrMsg,'UpdateSettingsWithCL')
            ENDIF

         ENDIF


      ELSEIF   ( CLFlags%YRange .AND. CLFlags%Dy ) THEN

            ! Set flags if not already set
         IF ( .NOT. DvrFlags%Dy )         DvrFlags%Dy       =  .TRUE.
         IF ( .NOT. DvrFlags%YRange )     DvrFlags%YRange   =  .TRUE.
         IF ( .NOT. WindGridModify )      WindGridModify    =  .TRUE.

            ! Copy over the range and stepsize
         DvrSettings%YRange(1)      =  minval(CLSettings%YRange)     ! just in case the order is funky
         DvrSettings%YRange(2)      =  maxval(CLSettings%YRange)
         DvrSettings%GridDelta(2)   =  abs(CLSettings%GridDelta(2))  ! just in case there was a negative number

            ! Set number of points in Y direction
         DvrSettings%GridN(2) =  CEILING( abs( ( DvrSettings%YRange(2) - DvrSettings%YRange(1) ) / DvrSettings%GridDelta(2) ) + 1_IntKi )

            ! Now adjust the upper limit of the range slightly for the integer number of steps
         DvrSettings%YRange(2)=  DvrSettings%YRange(1)   + ( DvrSettings%GridN(2) - 1_IntKi ) * DvrSettings%GridDelta(2)

            ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
         IF ( .NOT. EqualRealNos(DvrSettings%YRange(2),maxval(CLSettings%YRange)) )  THEN
            CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of YRange to '//TRIM(Num2LStr(DvrSettings%YRange(2)))//   &
               ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(2)))//' points with a requested spacing '//       &
               TRIM(Num2LStr(DvrSettings%GridDelta(2)))//'.',                                                           &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ENDIF

      ELSE
         ! Nothing was specified, so we don't do anything
      ENDIF




         ! ZRange / Dz  -- This is processed exactly the same as for YRange.




         !--------------------------------------------
         ! Did we change the FFT info?
         !--------------------------------------------

         !--------------------------------------------
         ! Did we request a summary file?
         !--------------------------------------------

         !--------------------------------------------
         ! Did we request a different Points file?
         !--------------------------------------------




         ! Did we modify any of the WindGrid information?
      IF ( WindGridModify ) THEN
         !FIXME: do something here
      ENDIF



   ELSE  ! No driver input file was read, or WindGrid was not set, so we have to populate a few extra things.

         ! Time settings


!FIXME: impose these restrictions either here or at the driver level.
         ! XRange -- if not specified on input, will need to assume only at x=0 (otherwise becomes too huge a dataset)

         ! YRange

         ! ZRange

         ! FFT info

         ! Sum file

         ! Points file         


   ENDIF



END SUBROUTINE UpdateSettingsWithCL


!SUBROUTINE WaveElevGrid_Output (drvrInitInp, HDynInitInp, HDynInitOut, HDyn_p, ErrStat, ErrMsg)
!
!   TYPE(HD_drvr_InitInput),       INTENT( IN    )   :: drvrInitInp
!   TYPE(HydroDyn_InitInputType),  INTENT( IN    )   :: HDynInitInp
!   TYPE(HydroDyn_InitOutputType), INTENT( IN    )   :: HDynInitOut          ! Output data from initialization
!   TYPE(HydroDyn_ParameterType),  INTENT( IN    )   :: HDyn_p               ! Output data from initialization
!   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
!   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
!
!         ! Temporary local variables
!   INTEGER(IntKi)                                   :: ErrStatTmp           !< Temporary variable for the status of error message
!   CHARACTER(*)                                     :: ErrMsgTmp            !< Temporary variable for the error message
!
!   INTEGER(IntKi)                                   :: WaveElevFileUn       !< Number for the output file for the wave elevation series
!   CHARACTER(1024)                                  :: WaveElevFileName     !< Name for the output file for the wave elevation series
!   CHARACTER(128)                                   :: WaveElevFmt          !< Format specifier for the output file for wave elevation series
!
!
!   WaveElevFmt = "(F14.7,3x,F14.7,3x,F14.7)"
!
!   ErrMsg      = ""
!   ErrStat     = ErrID_None
!   ErrMsgTmp   = ""
!   ErrStatTmp  = ErrID_None
!
!
!      ! If we calculated the wave elevation at a set of coordinates for use with making movies, put it into an output file
!   WaveElevFileName  =  TRIM(drvrInitInp%OutRootName)//".WaveElev.out"
!   CALL GetNewUnit( WaveElevFileUn )
!
!   CALL OpenFOutFile( WaveElevFileUn, WaveElevFileName, ErrStat, ErrMsg )
!   IF ( ErrStat /= ErrID_None) THEN
!      IF ( ErrStat >= AbortErrLev ) RETURN
!   ENDIF
!
!      ! Write some useful header information
!!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file was generated by '//TRIM(GetNVD(HDyn_Drv_ProgDesc))// &
!!         ' on '//CurDate()//' at '//CurTime()//'.'
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file was generated on '//CurDate()//' at '//CurTime()//'.'
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file contains the wave elevations at a series of points '// &
!         'through the entire timeseries.'
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## It is arranged as blocks of X,Y,Elevation at each timestep'
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## Each block is separated by two blank lines for use in gnuplot'
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# '
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# WaveTMax    =  '//TRIM(Num2LStr(HDyn_p%WaveTime(HDyn_P%NStepWave)))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# NStepWave   =  '//TRIM(Num2LStr(HDyn_p%NStepWave))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridXPoints =  '//TRIM(Num2LStr(drvrInitInp%WaveElevNX))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridYPoints =  '//TRIM(Num2LStr(drvrInitInp%WaveElevNY))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridDX      =  '//TRIM(Num2LStr(drvrInitInp%WaveElevDX))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridDY      =  '//TRIM(Num2LStr(drvrInitInp%WaveElevDY))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# MaxWaveElev =  '//TRIM(Num2LStr(MAXVAL(HDynInitOut%WaveElevSeries)))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# MinWaveElev =  '//TRIM(Num2LStr(MINVAL(HDynInitOut%WaveElevSeries)))
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# '
!
!      ! Timestep looping
!   DO I = 0,HDyn_p%NStepWave
!      WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp ) NewLine
!      WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp ) '# Time: '//TRIM(Num2LStr(HDyn_p%WaveTime(I)))
!         ! Now output the X,Y, Elev info for this timestep
!      DO J=1,SIZE(HDynInitInp%WaveElevXY,DIM=2)
!         WRITE (WaveElevFileUn,WaveElevFmt, IOSTAT=ErrStatTmp ) HDynInitInp%WaveElevXY(1,J),&
!                  HDynInitInp%WaveElevXY(2,J),HDynInitOut%WaveElevSeries(I,J)
!      ENDDO
!
!   ENDDO
!
!      ! Done.  Close the file
!   CLOSE (WaveElevFileUn)
!
!END SUBROUTINE WaveElevGrid_Output




!> This routine exists only to support the development of the module.  It will not be kept after the module is complete.
SUBROUTINE  printSettings( DvrFlags, DvrSettings )
      ! The arguments
   TYPE( IfWDriver_Flags ),            INTENT(IN   )  :: DvrFlags           !< Flags indicating which settings were set
   TYPE( IfWDriver_Settings ),         INTENT(IN   )  :: DvrSettings        !< Stored settings
   print*,' DvrIptFile:          ',FLAG(DvrFlags%DvrIptFile),        '      ',TRIM(DvrSettings%DvrIptFileName)
   print*,' IfWIptFile:          ',FLAG(DvrFlags%IfWIptFile),        '      ',TRIM(DvrSettings%IfwIptFileName)
   print*,' Summary:             ',FLAG(DvrFlags%Summary)
   print*,' SummaryFile:         ',FLAG(DvrFlags%SummaryFile),       '      ',TRIM(DvrSettings%SummaryFileName)
   print*,' TStart:              ',FLAG(DvrFlags%TStart),            '      ',DvrSettings%TStart
   print*,' DT:                  ',FLAG(DvrFlags%DT),                '      ',DvrSettings%DT
   print*,' NumTimeSteps:        ',FLAG(DvrFlags%NumTimeSteps),      '      ',DvrSettings%NumTimeSteps
   print*,' NumTimeStepsDefault: ',FLAG(DvrFlags%NumTimeStepsDefault)
   print*,' XRange:              ',FLAG(DvrFlags%XRange),            '      ',DvrSettings%XRange
   print*,' YRange:              ',FLAG(DvrFlags%YRange),            '      ',DvrSettings%YRange
   print*,' ZRange:              ',FLAG(DvrFlags%ZRange),            '      ',DvrSettings%ZRange
   print*,' Dx:                  ',FLAG(DvrFlags%Dx),                '      ',DvrSettings%GridDelta(1)
   print*,' Dy:                  ',FLAG(DvrFlags%Dy),                '      ',DvrSettings%GridDelta(2)
   print*,' Dz:                  ',FLAG(DvrFlags%Dz),                '      ',DvrSettings%GridDelta(3)
   print*,' FFTcalc:             ',FLAG(DvrFlags%FFTcalc),           '     [',DvrSettings%FFTcoord(1),',',DvrSettings%FFTcoord(2),',',DvrSettings%FFTcoord(3),']'
   print*,' WindGrid:            ',FLAG(DvrFlags%WindGrid)
   print*,' GridN:               ',DvrSettings%GridN
   RETURN
END SUBROUTINE printSettings

!> This routine exists only to support the development of the module.  It will not be kept after the module is complete.
!! This routine takes a flag setting (LOGICAL) and exports either 'T' or '-' for T/F (respectively)
FUNCTION FLAG(flagval)
   LOGICAL,       INTENT(IN   )  :: flagval        !< Value of the flag
   CHARACTER(1)                  :: FLAG           !< character interpretation (for prettiness when printing)
   IF ( flagval ) THEN
      FLAG  =  'T'
   ELSE
      FLAG  =  '-'
   ENDIF
   RETURN
END FUNCTION FLAG

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
END MODULE InflowWind_Driver_Subs
