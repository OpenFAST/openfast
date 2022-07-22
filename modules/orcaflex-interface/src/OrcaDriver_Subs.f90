!**********************************************************************************************************************************
!
!  MODULE: OrcaDriver_Subs  - This module contains subroutines used by the OrcaFlexInterface Driver program
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of OrcaFlexInterface.
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
! (File) Revision #: $Rev: 169 $
! URL: $HeadURL: https://wind-dev.nrel.gov/svn/OrcaFlexInterface/Trunk/Source/Driver/OrcaDriver_Subs.f90 $
!**********************************************************************************************************************************
MODULE OrcaDriver_Subs

   USE NWTC_Library
   USE OrcaDriver_Types
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
   CALL WrScr("")
   CALL WrScr("              The following options will overwrite values in the driver input file:")
   CALL WrScr("                  "//SwChar//"DT[#]         -- timestep                                         ")
   CALL WrScr("                  "//SwChar//"degrees       -- input angles specified in degrees                ")
   CALL WrScr("                  "//SwChar//"pointsdegrees -- input angles in points filespecified in degrees  ")
   CALL WrScr("                  "//SwChar//"Coord[X,Y,Z,R1,R2,R3]                                             ")
   CALL WrScr("                                 -- platform origin centered at [X,Y,Z]                         ")
   CALL WrScr("                                    with Roll / Pitch / Yaw of [R1,R2,R3]                       ")
   CALL WrScr("                  "//SwChar//"points[FILE]  -- calculates at a given position specified in a    ")
   CALL WrScr("                                    comma delimited FILE.")
   CALL WrScr("                  "//SwChar//"v             -- increase verbose level to 7 ")
   CALL WrScr("                  "//SwChar//"vv            -- increase verbose level to 10 ")
   CALL WrScr("                  "//SwChar//"help          -- print this help menu and exit")
   CALL WrScr("")
   CALL WrScr("   Notes:")
   CALL WrScr("   -- Options are not case sensitive.")
   CALL WrScr("   -- If no coordinates are specified, assumed to be at (0,0,0) with no rotation")
   CALL WrScr("")


END SUBROUTINE DispHelpText


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine retrieves the command line arguments and passes them to the
!! ::ParseArg routine for processing.
SUBROUTINE RetrieveArgs( CLSettings, CLFlags, ErrStat, ErrMsg )

   USE NWTC_Library
   USE OrcaDriver_Types

   IMPLICIT NONE

      ! Storing the arguments
   TYPE( OrcaDriver_Flags ),            INTENT(  OUT)  :: CLFlags      !< Flags indicating which command line arguments were specified
   TYPE( OrcaDriver_Settings ),         INTENT(  OUT)  :: CLSettings           !< Command line arguments passed in

      ! Error Handling
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

      ! Local variable
   INTEGER(IntKi)                                     :: i                    !< Generic counter
   CHARACTER(1024)                                    :: Arg                  !< argument given
   CHARACTER(1024)                                    :: ArgUC                !< Upper case argument to check
   INTEGER(IntKi)                                     :: NumInputArgs         !< Number of argements passed in from command line
   LOGICAL                                            :: ifwFlag              !< The -ifw flag was set
   CHARACTER(1024)                                    :: FileName             !< Filename from the command line.
   LOGICAL                                            :: FileNameGiven        !< Flag indicating if a filename was given.

   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status (for calls)
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error message (for calls)
   CHARACTER(*),              PARAMETER               :: RoutineName = 'RetrieveArgs'


      ! initialize some things
   ErrStat                    =  ErrID_None
   ErrStatTmp                 =  ErrID_None
   ErrMsg                     =  ''
   ErrMsgTmp                  =  ''
   FileNameGiven              =  .FALSE.
   FileName                   =  ''


      ! Check how many arguments are passed in
   NumInputArgs = COMMAND_ARGUMENT_COUNT()

      ! exit if we don't have enough
   IF (NumInputArgs == 0) THEN
      CALL SetErrStat(ErrID_Fatal," Insufficient Arguments. Use option "//SwChar//"help for help menu.",  &
         ErrStat,ErrMsg,RoutineName)
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
         CALL ParseArg( CLSettings, CLFlags, ArgUC(2:), Arg(2:), ifwFlag, ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat>AbortErrLev) RETURN

      ELSE

            ! since there is no switch character, assume it is the filename, unless we already set one
         IF ( FileNameGiven ) THEN
            CALL SetErrStat(ErrID_Fatal," Multiple driver input filenames given: "//TRIM(FileName)//", "//TRIM(Arg),   &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ELSE
            FileName       = TRIM(Arg)
            FileNameGiven  = .TRUE.
         ENDIF

      ENDIF
   END DO


      ! Was a filename given?
   IF ( .NOT. FileNameGiven ) THEN
      CALL SetErrStat( ErrID_Fatal, " No filename given.", ErrStat, ErrMsg, RoutineName )
      RETURN
   ELSE
      CLSettings%DvrIptFileName  =  TRIM(FileName)
      CLFlags%DvrIptFile =  .TRUE.
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
   SUBROUTINE ParseArg( CLSettings, CLFlags, ThisArgUC, ThisArg, ifwFlagSet, ErrStat, ErrMsg )
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
         ! Parse and store the input argument  !
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

      USE NWTC_Library
      USE OrcaDriver_Types
      USE OrcaFlexInterface_Types

      IMPLICIT NONE

         ! Storing the arguments
      TYPE( OrcaDriver_Flags ),            INTENT(INOUT)  :: CLFlags      ! Flags indicating which arguments were specified
      TYPE( OrcaDriver_Settings ),         INTENT(INOUT)  :: CLSettings           ! Arguments passed in

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
      INTEGER(IntKi)                                     :: DelimSep4            ! where the : is
      INTEGER(IntKi)                                     :: DelimSep5            ! where the : is
      REAL(ReKi)                                         :: TempReal             ! temp variable to hold a real

      INTEGER(IntKi)                                     :: ErrStatTmp           ! Temporary error status for calls
      CHARACTER(*),              PARAMETER               :: RoutineName = 'ParseArg'



         ! Initialize some things
      ErrStat     = ErrID_None
      ErrStatTmp  =  ErrID_None
      ErrMsg      = ''

         ! Get the delimiters -- returns 0 if there isn't one
      Delim1   = INDEX(ThisArgUC,'[')
      Delim2   = INDEX(ThisArgUC,']')
      DelimSep = INDEX(ThisArgUC,':')


         ! check that if there is an opening bracket, then there is a closing one
      IF ( (Delim1 > 0_IntKi ) .and. (Delim2 < Delim1) ) THEN
         CALL SetErrStat(ErrID_Warn," Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
            ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF

         ! check that if there is a colon, then there are brackets
      IF ( (DelimSep > 0_IntKi) .and. (Delim1 == 0_IntKi) ) THEN
         CALL SetErrStat(ErrID_Warn," Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
            ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF


         ! If no delimeters were given, than this option is simply a flag
      IF ( Delim1 == 0_IntKi ) THEN
            ! check to see if the filename is the name of the Orca input file
         IF       ( ThisArgUC(1:13)== "POINTSDEGREES" )  THEN
            CLFlags%PointsDegrees   = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:13)== "AddedMassFile" )  THEN
            CLFlags%AddedMassFile   = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:9) == "AddedMass"     )  THEN
            CLFlags%AddedMass       = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:7) == "DEGREES"       )  THEN
            CLFlags%Degrees         = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:2) == "VV"            )  THEN
            CLFlags%VVerbose        = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:1) == "V"             )  THEN
            CLFlags%Verbose         = .TRUE.
            RETURN
         ELSE
            CALL SetErrStat( ErrID_Warn," Unrecognized option '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options.",  &
               ErrStat,ErrMsg,RoutineName)
         ENDIF

      ENDIF



          ! "Veloc[X,Y,Z,R1,R2,R3]"
      IF     ( ThisArgUC(1:Delim1) == "Veloc["  ) THEN
         DelimSep = INDEX(ThisArgUC,',')
         DelimSep2= INDEX(ThisArgUC(DelimSep+1:),',') + DelimSep
         IF ( DelimSep2 <= DelimSep ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep3= INDEX(ThisArgUC(DelimSep2+1:),',') + DelimSep
         IF ( DelimSep3 <= DelimSep2 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep4= INDEX(ThisArgUC(DelimSep3+1:),',') + DelimSep
         IF ( DelimSep4 <= DelimSep3 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep5= INDEX(ThisArgUC(DelimSep4+1:),',') + DelimSep
         IF ( DelimSep5 <= DelimSep4 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep5= INDEX(ThisArgUC(DelimSep5+1:),',') + DelimSep
         IF ( DelimSep5 <= DelimSep5 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmVeloc            = .TRUE.
            CLSettings%PtfmVeloc(1)     = TempReal
         ELSE
            CLFlags%PtfmVeloc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat(ErrID_FATAL," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:DelimSep2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmVeloc            = .TRUE.
            CLSettings%PtfmVeloc(2)     = TempReal
         ELSE
            CLFlags%PtfmVeloc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Third Value
         TempReal = StringToReal( ThisArgUC(DelimSep2+1:DelimSep3-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmVeloc            = .TRUE.
            CLSettings%PtfmVeloc(3)     = TempReal
         ELSE
            CLFlags%PtfmVeloc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Fourth Value
         TempReal = StringToReal( ThisArgUC(DelimSep3+1:DelimSep4-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmVeloc            = .TRUE.
            CLSettings%PtfmVeloc(4)     = TempReal
         ELSE
            CLFlags%PtfmVeloc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Fifth Value
         TempReal = StringToReal( ThisArgUC(DelimSep4+1:DelimSep5-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmVeloc            = .TRUE.
            CLSettings%PtfmVeloc(5)     = TempReal
         ELSE
            CLFlags%PtfmVeloc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Sixth Value
         TempReal = StringToReal( ThisArgUC(DelimSep5+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmVeloc            = .TRUE.
            CLSettings%PtfmVeloc(6)     = TempReal
         ELSE
            CLFlags%PtfmVeloc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF



          ! "Accel[X,Y,Z,R1,R2,R3]"
      ELSEIF ( ThisArgUC(1:Delim1) == "Accel["  ) THEN
         DelimSep = INDEX(ThisArgUC,',')
         DelimSep2= INDEX(ThisArgUC(DelimSep+1:),',') + DelimSep
         IF ( DelimSep2 <= DelimSep ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep3= INDEX(ThisArgUC(DelimSep2+1:),',') + DelimSep
         IF ( DelimSep3 <= DelimSep2 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep4= INDEX(ThisArgUC(DelimSep3+1:),',') + DelimSep
         IF ( DelimSep4 <= DelimSep3 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep5= INDEX(ThisArgUC(DelimSep4+1:),',') + DelimSep
         IF ( DelimSep5 <= DelimSep4 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep5= INDEX(ThisArgUC(DelimSep5+1:),',') + DelimSep
         IF ( DelimSep5 <= DelimSep5 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmAccel            = .TRUE.
            CLSettings%PtfmAccel(1)     = TempReal
         ELSE
            CLFlags%PtfmAccel            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat(ErrID_FATAL," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:DelimSep2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmAccel            = .TRUE.
            CLSettings%PtfmAccel(2)     = TempReal
         ELSE
            CLFlags%PtfmAccel            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Third Value
         TempReal = StringToReal( ThisArgUC(DelimSep2+1:DelimSep3-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmAccel            = .TRUE.
            CLSettings%PtfmAccel(3)     = TempReal
         ELSE
            CLFlags%PtfmAccel            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Fourth Value
         TempReal = StringToReal( ThisArgUC(DelimSep3+1:DelimSep4-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmAccel            = .TRUE.
            CLSettings%PtfmAccel(4)     = TempReal
         ELSE
            CLFlags%PtfmAccel            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Fifth Value
         TempReal = StringToReal( ThisArgUC(DelimSep4+1:DelimSep5-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmAccel            = .TRUE.
            CLSettings%PtfmAccel(5)     = TempReal
         ELSE
            CLFlags%PtfmAccel            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Sixth Value
         TempReal = StringToReal( ThisArgUC(DelimSep5+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmAccel            = .TRUE.
            CLSettings%PtfmAccel(6)     = TempReal
         ELSE
            CLFlags%PtfmAccel            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

          ! "Coord[X,Y,Z,R1,R2,R3]"
      ELSEIF ( ThisArgUC(1:Delim1) == "Coord["  ) THEN
         DelimSep = INDEX(ThisArgUC,',')
         DelimSep2= INDEX(ThisArgUC(DelimSep+1:),',') + DelimSep
         IF ( DelimSep2 <= DelimSep ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep3= INDEX(ThisArgUC(DelimSep2+1:),',') + DelimSep
         IF ( DelimSep3 <= DelimSep2 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep4= INDEX(ThisArgUC(DelimSep3+1:),',') + DelimSep
         IF ( DelimSep4 <= DelimSep3 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep5= INDEX(ThisArgUC(DelimSep4+1:),',') + DelimSep
         IF ( DelimSep5 <= DelimSep4 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF
         DelimSep5= INDEX(ThisArgUC(DelimSep5+1:),',') + DelimSep
         IF ( DelimSep5 <= DelimSep5 ) THEN
            CALL SetErrStat(ErrID_Warn," Unrecognized coordinate in '"//SwChar//TRIM(ThisArg)//"'.  Ignoring.", &
               ErrStat,ErrMsg,RoutineName)
            RETURN
         ENDIF

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmCoord            = .TRUE.
            CLSettings%PtfmCoord(1)     = TempReal
         ELSE
            CLFlags%PtfmCoord            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat(ErrID_FATAL," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArgUC(DelimSep+1:DelimSep2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmCoord            = .TRUE.
            CLSettings%PtfmCoord(2)     = TempReal
         ELSE
            CLFlags%PtfmCoord            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat(ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",  &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Third Value
         TempReal = StringToReal( ThisArgUC(DelimSep2+1:DelimSep3-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmCoord            = .TRUE.
            CLSettings%PtfmCoord(3)     = TempReal
         ELSE
            CLFlags%PtfmCoord            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Fourth Value
         TempReal = StringToReal( ThisArgUC(DelimSep3+1:DelimSep4-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmCoord            = .TRUE.
            CLSettings%PtfmCoord(4)     = TempReal
         ELSE
            CLFlags%PtfmCoord            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Fifth Value
         TempReal = StringToReal( ThisArgUC(DelimSep4+1:DelimSep5-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmCoord            = .TRUE.
            CLSettings%PtfmCoord(5)     = TempReal
         ELSE
            CLFlags%PtfmCoord            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF

            ! Sixth Value
         TempReal = StringToReal( ThisArgUC(DelimSep5+1:Delim2-1), ErrStatTmp )
         IF ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%PtfmCoord            = .TRUE.
            CLSettings%PtfmCoord(6)     = TempReal
         ELSE
            CLFlags%PtfmCoord            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, RoutineName)
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, RoutineName)
            ENDIF
            RETURN
         ENDIF


         ! "POINTS[FILE]"
      ELSEIF( ThisArgUC(1:Delim1)   == "POINTS["    ) THEN
         CLFlags%PointsFile         = .TRUE.
         CLSettings%PointsFileName  = ThisArg(Delim1+1:Delim2-1)
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
   TYPE(OrcaDriver_Flags),             INTENT(INOUT)  :: DvrFlags
   TYPE(OrcaDriver_Settings),          INTENT(INOUT)  :: DvrSettings
   TYPE(ProgDesc),                     INTENT(IN   )  :: ProgInfo
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)                                     :: UnIn                 ! Unit number for the driver input file
   CHARACTER(1024)                                    :: FileName             ! Name of OrcaFlexInterface driver input file

      ! Input file echoing
   LOGICAL                                            :: EchoFileContents     ! Do we echo the driver file out or not?
   INTEGER(IntKi)                                     :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                    :: EchoFileName         ! Name of OrcaFlexInterface driver echo file

      ! Time steps
   CHARACTER(1024)                                    :: DTChr                ! Character string for timesteps size (to handle DEFAULT value)

      ! Local error handling
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   INTEGER(IntKi)                                     :: ErrStatTmp2          !< Temporary error status for IO checks
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error messages for calls
   CHARACTER(*),              PARAMETER               :: RoutineName = 'ReadDvrIptFile'

      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1

   FileName = TRIM(DvrFileName)

   CALL GetNewUnit( UnIn )
   CALL OpenFInpFile( UnIn, FileName, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL WrScr( 'Opening OrcaFlexInterface Driver input file:  '//FileName )


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName,' OrcaFlexInterface Driver input file header line 1', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL ReadCom( UnIn, FileName, 'OrcaFlexInterface Driver input file header line 2', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CLOSE( UnIn )
      RETURN
   ENDIF


     ! Echo Input Files.
   CALL ReadVar ( UnIn, FileName, EchoFileContents, 'Echo', 'Echo Input', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
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
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CLOSE( UnIn )
         RETURN
      ENDIF


      REWIND(UnIn)


         ! Reread and echo
      CALL ReadCom( UnIn, FileName,' OrcaFlexInterface Driver input file header line 1', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF


      CALL ReadCom( UnIn, FileName, 'OrcaFlexInterface Driver input file header line 2', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF


        ! Echo Input Files.
      CALL ReadVar ( UnIn, FileName, EchoFileContents, 'Echo', 'Echo Input', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF



   ENDIF


   !-------------------------------------------------------------------------------------------------
   !  OrcaFlexInterface setup section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' OrcaFlexInterface setup section, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! DT    -- Timestep size for the driver to take (or DEFAULT for what the file contains)
   CALL ReadVar( UnIn, FileName,DTChr,'DTChr',' Character string for Timestep size for the driver to take (or DEFAULT for what the file contains).',  &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
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
      READ (DTChr,*,IOSTAT=ErrStatTmp2)   DvrSettings%DT
      IF ( ErrStatTmp /= ErrID_None )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( ErrStatTmp2, '', 'DT',NumType, ErrStatTmp, ErrMsgTmp, .TRUE. )
         RETURN
      ELSE     ! Was ok, so set the flags
         DvrFlags%DT         =  .TRUE.
         DvrFlags%DTDefault  =  .FALSE.
      ENDIF
   ENDIF


      ! OrcaFlexInterface input file
   CALL ReadVar( UnIn, FileName,DvrSettings%OrcaIptFileName,'OrcaIptFileName',' OrcaFlexInterface input filename',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%OrcaIptFile  =  .TRUE.
   ENDIF



   !-------------------------------------------------------------------------------------------------
   !  PtfmCoordinates
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' Coordinates, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


       ! PtfmCoord    -- PtfmCoord of the windfield needed.
   CALL ReadVar( UnIn, FileName,DvrFlags%PtfmCoord,'PtfmCoord',' Use a set of coordinates?',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   DvrFlags%PtfmVeloc   =  DvrFlags%PtfmCoord
   DvrFlags%PtfmAccel   =  DvrFlags%PtfmCoord
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! Read the coordinates if the flag is set, otherwise skip the line
   IF ( DvrFlags%PtfmCoord   ) THEN
   
         ! Degrees    -- PtfmCoord of the windfield needed.
      CALL ReadVar( UnIn, FileName,DvrFlags%Degrees,'Degrees',' Angles specified in degrees?',   &
         ErrStatTmp,ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   
         ! PtfmCoord     -- The coordinates to pass to the DLL
      CALL ReadAry ( UnIn, FileName, DvrSettings%PtfmCoord(1:6), 6, 'PtfmCoord(1:6)', &
         'platform coordinate', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   
         ! PtfmVeloc     -- The coordinates to pass to the DLL
      CALL ReadAry ( UnIn, FileName, DvrSettings%PtfmVeloc(1:6), 6, 'PtfmVeloc(1:6)', &
         'platform coordinate', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   
         ! PtfmAccel     -- The coordinates to pass to the DLL
      CALL ReadAry ( UnIn, FileName, DvrSettings%PtfmAccel(1:6), 6, 'PtfmAccel(1:6)', &
         'platform coordinate', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   
         ! Write the added mass matrix to the screen
      CALL ReadVar( UnIn, FileName,DvrFlags%AddedMass,'AddedMass',' Write table of added mass to screen.', &
         ErrStatTmp,ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ELSE
         DvrFlags%AddedMass  =  .TRUE.
      ENDIF
   
         ! Write the added mass matrix to a file
      CALL ReadVar( UnIn, FileName,DvrFlags%AddedMassFile,'AddedMassFile',' Write added Mass matrix to file.', &
         ErrStatTmp,ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ELSE
         DvrFlags%AddedMassFile =  .TRUE.
      ENDIF
   
   ELSE
      CALL ReadCom( UnIn, FileName,' Skipping the degrees flag since not calculating anything.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the platform coordinate since not calculating anything.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the platform velocity since not calculating anything.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the platform acceleration since not calculating anything.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the Added mass matrix output since not calculating anything.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the added mass matrix file output since not calculating anything.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
   ENDIF




   !-------------------------------------------------------------------------------------------------
   !  points file input
   !-------------------------------------------------------------------------------------------------

      ! Header line
   CALL ReadCom( UnIn, FileName,' Points file input, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! PointsFile    -- Read a points file
   CALL ReadVar( UnIn, FileName,DvrFlags%PointsFile,'PointsFile',' Read a points file?',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


   IF ( DvrFlags%PointsFile ) THEN
         ! Points file in degrees
      CALL ReadVar( UnIn, FileName,DvrFlags%PointsDegrees,'PointsDegrees',' Angles in points file given in degrees?',   &
         ErrStatTmp,ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
         ! Points input file
      CALL ReadVar( UnIn, FileName,DvrSettings%PointsFileName,'PointsFileName',' Points file input filename',   &
         ErrStatTmp,ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   ELSE
         ! Skip the next entry points file section.
      CALL ReadCom( UnIn, FileName,' Skipping the degreespoints flag since not using it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
         ! Skip the next entry points file section.
      CALL ReadCom( UnIn, FileName,' Skipping the points filename since not using it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
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

   TYPE(OrcaDriver_Flags),             INTENT(INOUT)  :: DvrFlags
   TYPE(OrcaDriver_Settings),          INTENT(INOUT)  :: DvrSettings
   TYPE(OrcaDriver_Flags),             INTENT(IN   )  :: CLFlags
   TYPE(OrcaDriver_Settings),          INTENT(IN   )  :: CLSettings
   LOGICAL,                            INTENT(IN   )  :: DVRIPT
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


      ! Local variables
   INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp         !< Temporary error status for calls
   CHARACTER(*),              PARAMETER               :: RoutineName = 'UpdateSettingsWithCL'
   LOGICAL                                            :: WindGridModify    !< Did we modify any of the WindGrid related settings?

   INTEGER(IntKi)                                     :: I                 !< local counter

      ! Initialization
   WindGridModify =  .FALSE.

      ! Initialize the error handling
   ErrStat     =  ErrID_None
   ErrMsg      =  ''
   ErrStatTmp  =  ErrID_None
   ErrMsgTmp   =  ''


      !--------------------------------------------
      !  Did we change any time information?
      !--------------------------------------------


      ! Check DT
   IF ( CLFlags%DT ) THEN
      IF ( DvrFlags%DT .AND. ( .NOT. EqualRealNos(DvrSettings%DT, CLSettings%DT) ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for DT with '//TRIM(Num2LStr(CLSettings%DT))//'.',  &
            ErrStat,ErrMsg,RoutineName)
      ELSE
         DvrFlags%DT   =  .TRUE.
      ENDIF
      DvrSettings%DT   =  CLSettings%DT
   ENDIF


      !--------------------------------------------
      ! Did we change the coordinate info?
      !--------------------------------------------

   IF ( CLFlags%PtfmCoord ) THEN
        ! If we are overriding driver input file settings, tell user
      IF ( DvrFlags%PtfmCoord ) THEN
         CALL SetErrStat( ErrID_Warn,' Overriding driver input file settings for platform coordinate.',   &
            ErrStat,ErrMsg,RoutineName )
      ENDIF
      DvrSettings%PtfmCoord =  CLSettings%PtfmCoord
      DvrFlags%PtfmCoord     =  .TRUE.
   ENDIF

   IF ( CLFlags%PtfmVeloc ) THEN
        ! If we are overriding driver input file settings, tell user
      IF ( DvrFlags%PtfmVeloc ) THEN
         CALL SetErrStat( ErrID_Warn,' Overriding driver input file settings for platform velocities.',   &
            ErrStat,ErrMsg,RoutineName )
      ENDIF
      DvrSettings%PtfmVeloc =  CLSettings%PtfmVeloc
      DvrFlags%PtfmVeloc     =  .TRUE.
   ENDIF

   IF ( CLFlags%PtfmAccel ) THEN
        ! If we are overriding driver input file settings, tell user
      IF ( DvrFlags%PtfmAccel ) THEN
         CALL SetErrStat( ErrID_Warn,' Overriding driver input file settings for platform velocities.',   &
            ErrStat,ErrMsg,RoutineName )
      ENDIF
      DvrSettings%PtfmAccel =  CLSettings%PtfmAccel
      DvrFlags%PtfmAccel     =  .TRUE.
   ENDIF

      ! If only one of the PtfmCoord, PtfmVeloc, or PtfmAccel flags is set to true, the other should be set also
   IF ( DvrFlags%PtfmCoord .OR. DvrFlags%PtfmVeloc .OR. DvrFlags%PtfmAccel ) THEN
      DvrFlags%PtfmCoord   =  .TRUE.
      DvrFlags%PtfmVeloc   =  .TRUE.
      DvrFlags%PtfmAccel   =  .TRUE.
   ENDIF


      !--------------------------------------------
      ! Are PtfmCoord angles in degrees?
      !--------------------------------------------

   IF ( CLFlags%Degrees ) THEN
         ! No need to tell the user.  They likely only specified this flag with command line coords.
         ! The logic needed to check if it is a status change otherwise is not worth the effort.
      DvrFlags%Degrees  =  .TRUE.
   ENDIF


      !--------------------------------------------
      ! Did we request Added Mass matrix results?
      !--------------------------------------------

   IF ( CLFlags%AddedMass ) THEN
      IF ( DvrFlags%PtfmCoord .AND. DvrFlags%PtfmVeloc .AND. DvrFlags%PtfmAccel ) THEN
         DvrFlags%AddedMass   =  .TRUE.
      ELSE  ! give a warning and set the flags.  The coordinate is already initialized to (0,0,0,0,0,0)
         CALL SetErrStat( ErrID_Warn,' Added mass matrix requested, but no platform location specified.  Setting location to (0,0,0,0,0,0).',  &
            ErrStat,ErrMsg,RoutineName)
         DvrFlags%AddedMass      =  .FALSE.
         DvrFlags%PtfmCoord      =  .TRUE.
         DvrFlags%PtfmVeloc      =  .TRUE.
         DvrFlags%PtfmAccel      =  .TRUE.
      ENDIF
   ENDIF

   IF ( CLFlags%AddedMassFile ) THEN
      IF ( DvrFlags%PtfmCoord .AND. DvrFlags%PtfmVeloc .AND. DvrFlags%PtfmAccel ) THEN
         DvrFlags%AddedMassFile   =  .TRUE.
      ELSE  ! give a warning and set the flags.  The coordinate is already initialized to (0,0,0,0,0,0)
         CALL SetErrStat( ErrID_Warn,' Added mass matrix file requested, but no platform location specified.  Setting location to (0,0,0,0,0,0).',  &
            ErrStat,ErrMsg,RoutineName)
         DvrFlags%AddedMassFile  =  .TRUE.
         DvrFlags%PtfmCoord      =  .TRUE.
         DvrFlags%PtfmVeloc      =  .TRUE.
         DvrFlags%PtfmAccel      =  .TRUE.
      ENDIF
   ENDIF



      !--------------------------------------------
      ! Did we request a different Points file?
      !--------------------------------------------

   IF ( CLFlags%PointsFile ) THEN
         ! If a name was given in the driver input file, then warn the user.
      IF ( DvrFlags%PointsFile ) THEN
         CALL SetErrStat( ErrID_Warn,' Overriding driver input file settings for Points file.',  &
            ErrStat,ErrMsg,RoutineName )
      ENDIF
      DvrFlags%PointsFile        =  .TRUE.
      DvrSettings%PointsFileName =  CLSettings%PointsFileName
      IF ( CLFlags%PointsDegrees ) THEN
         DvrFLags%PointsDegrees = .TRUE.
      ENDIF
   ELSE
      IF ( CLFlags%PointsDegrees ) THEN
         DvrFlags%PointsDegrees = .TRUE.
         CALL SetErrStat( ErrID_Warn,' Overriding driver input file points file angles in degrees.',   &
            ErrStat,ErrMsg,RoutineName )
      ENDIF
   ENDIF


      ! If no DT value has been set (DEFAULT requested), we need to set a default to pass into Orca
   IF ( .NOT. DvrFlags%DT ) THEN
      DvrSettings%DT =  0.025_DbKi     ! This value gets passed into the Orca_Init routine, so something must be set.
      DvrFlags%DT    =  .TRUE.
   ENDIF


      ! If the angles for PtfmCoord, PtfmVeloc, and PtfmAccel are in degrees, then convert them into radians now
   IF ( DvrFlags%Degrees ) THEN
      DO I=4,6
         DvrSettings%PtfmCoord(I)   =  DvrSettings%PtfmCoord(I)   *  D2R      ! D2R is from the library
         DvrSettings%PtfmVeloc(I)   =  DvrSettings%PtfmVeloc(I)   *  D2R      ! D2R is from the library
         DvrSettings%PtfmAccel(I)   =  DvrSettings%PtfmAccel(I)   *  D2R      ! D2R is from the library
      ENDDO
   ENDIF


END SUBROUTINE UpdateSettingsWithCL


SUBROUTINE ReadPointsFile( PointsFileName, AnglesInDegrees, TimeList, CoordList, VelocList, AccelList, ErrStat, ErrMsg )

   CHARACTER(1024),                    INTENT(IN   )  :: PointsFileName       !< Name of the points file to read
   LOGICAL,                            INTENT(IN   )  :: AnglesInDegrees      !< Are the angles specified in degrees?
   REAL(ReKi), ALLOCATABLE,            INTENT(  OUT)  :: TimeList(:)          !< TimeStamps
   REAL(ReKi), ALLOCATABLE,            INTENT(  OUT)  :: CoordList(:,:)       !< The coordinates we read in
   REAL(ReKi), ALLOCATABLE,            INTENT(  OUT)  :: VelocList(:,:)       !< The velocities we read in
   REAL(ReKi), ALLOCATABLE,            INTENT(  OUT)  :: AccelList(:,:)       !< The accelerations we read in
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< The error status
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< The message for the status

      ! Local variables
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error message for calls
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   CHARACTER(*),              PARAMETER               :: RoutineName = 'ReadPointsFile'
   INTEGER(IntKi)                                     :: FiUnitPoints         !< Unit number for points file to open

   INTEGER(IntKi)                                     :: NumDataColumns       !< Number of data columns
   INTEGER(IntKi)                                     :: NumDataPoints        !< Number of lines of data (one point per line)
   INTEGER(IntKi)                                     :: NumHeaderLines       !< Number of header lines to ignore

   INTEGER(IntKi)                                     :: I                    !< Generic counter

   REAL(ReKi)                                         :: TmpArray(19)         !< Temporary array to hold one line of data from the points file
   REAL(ReKi)                                         :: ConvToRadians        !< Conversion to radians multiplier

      ! Initialization of subroutine
   ErrMsg      =  ''
   ErrMsgTmp   =  ''
   ErrStat     =  ErrID_None
   ErrStatTmp  =  ErrID_None

      ! Set the ConvToRadians multiplier
   IF ( AnglesInDegrees) THEN
      ConvToRadians  =  D2R      ! Set to the library constant
   ELSE
      ConvToRadians  =  1.0_ReKi
   ENDIF



      ! Now open file
   CALL GetNewUnit(    FiUnitPoints )
   CALL OpenFInpFile(   FiUnitPoints,  TRIM(PointsFileName), ErrStatTmp, ErrMsgTmp )   ! Unformatted input file
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF

      ! Find out how long the file is
   CALL GetFileLength( FiUnitPoints, PointsFileName, NumDataColumns, NumDataPoints, NumHeaderLines, ErrMsgTmp, ErrStatTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF
   IF ( NumDataColumns /= 19 ) THEN
      CALL SetErrStat( ErrID_Fatal,' Expecting 19 columns in '//TRIM(PointsFileName)//' corresponding to '//   &
         'timestamp, 3 translation coordinates, 3 rotation angles, 3 translational velocities, 3 rotational velocities, and '// &
         '3 translational velocities, 3 rotational velocities.  '//&
         'Instead found '//TRIM(Num2LStr(NumDataColumns))//' columns.', &
         ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Allocate the storage for the data
   CALL AllocAry( TimeList, NumDataPoints, "Array of timestamp data", ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Allocate the storage for the data
   CALL AllocAry( CoordList, 6, NumDataPoints, "Array of Points and rotation data", ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Read in the headers and throw them away
      ! Allocate the storage for the data
   CALL AllocAry( VelocList, 6, NumDataPoints, "Array of translation and rotation derivative data", ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Read in the headers and throw them away
      ! Allocate the storage for the data
   CALL AllocAry( AccelList, 6, NumDataPoints, "Array of translation and rotation 2nd derivative data", ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Read in the headers and throw them away
   DO I=1,NumHeaderLines
      CALL ReadCom( FiUnitPoints, PointsFileName,' Points file header line', ErrStatTmp, ErrMsgTmp )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CLOSE( FiUnitPoints )
         RETURN
      ENDIF
   ENDDO

      ! Read in the datapoints
   DO I=1,NumDataPoints
      CALL ReadAry ( FiUnitPoints, PointsFileName, TmpArray(:), 19, 'Temporary coordinate', &
         'Coordinate point from Points file', ErrStatTmp, ErrMsgTmp)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CLOSE( FiUnitPoints )
         RETURN
      ENDIF
      TimeList(I)      =  TmpArray(1)
      CoordList(1:3,I) =  TmpArray(2:4)
      CoordList(4:6,I) =  TmpArray(5:7)   * ConvToRadians
      VelocList(1:3,I) =  TmpArray(8:10)
      VelocList(4:6,I) =  TmpArray(11:13) * ConvToRadians
      AccelList(1:3,I) =  TmpArray(14:16)
      AccelList(4:6,I) =  TmpArray(17:19) * ConvToRadians
   ENDDO

   CLOSE( FiUnitPoints )

CONTAINS

  !-------------------------------------------------------------------------------------------------------------------------------
   !>    This subroutine looks at a file that has been opened and finds out how many header lines there are, how many columns there
   !!    are, and    how many lines of data there are in the file.
   !!
   !!    A few things are assumed about the file:
   !!       1. Any header lines are the first thing in the file.
   !!       2. No text appears anyplace other than in first part of the file
   !!       3. The datalines only contain numbers that can be read in as reals.
   !!
   !!    Limitations:
   !!       1. only handles up to 20 words (columns) on a line
   !!       2. empty lines are considered text lines
   !!       3. All data rows must contain the same number of columns
   !!
   !!
   SUBROUTINE GetFileLength(UnitDataFile, DataFileName, NumDataColumns, NumDataLines, NumHeaderLines, ErrMsg, ErrStat)

      IMPLICIT NONE

         ! Passed variables
      INTEGER(IntKi),                     INTENT(IN   )  :: UnitDataFile      !< Unit number of the file we are looking at.
      CHARACTER(*),                       INTENT(IN   )  :: DataFileName      !< The name of the file we are looking at.
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataColumns    !< The number of columns in the data file.
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataLines      !< Number of lines containing data
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumHeaderLines    !< Number of header lines at the start of the file
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< Error Message to return (empty if all good)
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< Status flag if there were any problems (ErrID_None if all good)

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp         !< Temporary message variable.  Used in calls.
      INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status.  Used in calls.
      CHARACTER(*),              PARAMETER               :: RoutineName = 'GetFileLength'
      INTEGER(IntKi)                                     :: LclErrStat        !< Temporary error status.  Used locally to indicate when we have reached the end of the file.
      INTEGER(IntKi)                                     :: TmpIOErrStat      !< Temporary error status for the internal read of the first word to a real number
      LOGICAL                                            :: IsRealNum         !< Flag indicating if the first word on the line was a real number

      CHARACTER(1024)                                    :: TextLine          !< One line of text read from the file
      INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
      CHARACTER(1024)                                    :: StrRead           !< String containing the first word read in
      REAL(ReKi)                                         :: RealRead          !< Returns value of the number (if there was one), or NaN (as set by NWTC_Num) if there wasn't
      CHARACTER(1024)                                    :: VarName           !< Name of the variable we are trying to read from the file
      CHARACTER(24)                                      :: Words(20)         !< Array of words we extract from a line.  We shouldn't have more than 20.
      INTEGER(IntKi)                                     :: i,j,k             !< simple integer counters
      INTEGER(IntKi)                                     :: LineNumber        !< the line I am on
      LOGICAL                                            :: LineHasText       !< Flag indicating if the line I just read has text.  If so, it is a header line.
      LOGICAL                                            :: HaveReadData      !< Flag indicating if I have started reading data.
      INTEGER(IntKi)                                     :: NumWords          !< Number of words on a line
      INTEGER(IntKi)                                     :: FirstDataLineNum  !< Line number of the first row of data in the file


         ! Initialize the error handling
      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None
      LclErrStat  = ErrID_None
      ErrMsg      = ''
      ErrMsgTmp   = ''


         ! Set some of the flags and counters
      HaveReadData   = .FALSE.
      NumDataColumns = 0
      NumHeaderLines = 0
      NumDataLines   = 0
      LineNumber     = 0


         ! Just in case we were handed a file that we are part way through reading (should never be true), rewind to the start

      REWIND( UnitDataFile )


         !------------------------------------
         !> The variable LclErrStat is used to indicate when we have reached the end of the file or had an error from
         !! ReadLine.  Until that occurs, we read each line, and decide if it contained any non-numeric data.  The
         !! first group of lines containing non-numeric data is considered the header.  The first line of all numeric
         !! data is considered the start of the data section.  Any non-numeric containing found within the data section
         !! will be considered as an invalid file format at which point we will return a fatal error from this routine.

      DO WHILE ( LclErrStat == ErrID_None )

            !> Reset the indicator flag for the non-numeric content
         LineHasText = .FALSE.

            !> Read in a single line from the file
         CALL ReadLine( UnitDataFile, '', TextLine, LineLen, LclErrStat )

            !> If there was an error in reading the file, then exit.
            !!    Possible causes: reading beyond end of file in which case we are done so don't process it.
         IF ( LclErrStat /= ErrID_None ) EXIT

            !> Increment the line counter.
         LineNumber  = LineNumber + 1

            !> Read all the words on the line into the array called 'Words'.  Only the first words will be encountered
            !! will be stored.  The others are empty (i.e. only three words on the line, so the remaining 17 are empty).
         CALL GetWords( TextLine, Words, 20 )

            !> Cycle through and count how many are not empty.  Once an empty value is encountered, all the rest should
            !! be empty if GetWords worked correctly.  The index of the last non-empty value is stored.
         DO i=1,20
            IF (TRIM(Words(i)) .ne. '') NumWords=i
         ENDDO


            !> Now cycle through the first 'NumWords' of non-empty values stored in 'Words'.  Words should contain
            !! everything that is one the line.  The subroutine ReadRealNumberFromString will set a flag 'IsRealNum'
            !! when the value in Words(i) can be read as a real(ReKi).  'StrRead' will contain the string equivalent.
         DO i=1,NumWords
            CALL ReadRealNumberFromString( Words(i), RealRead, StrRead, IsRealNum, ErrStatTmp, ErrMsgTmp, TmpIOErrStat )
            IF ( .NOT. IsRealNum) THEN
               LineHasText = .TRUE.
            ENDIF
         ENDDO

            !> If all the words on that line had no text in them, then it must have been a line of data.
            !! If not, then we have either a header line, which is ok, or a line containing text in the middle of the
            !! the data section, which is not good (the flag HaveReadData tells us which case this is).
         IF ( LineHasText ) THEN
            IF ( HaveReadData ) THEN      ! Uh oh, we have already read a line of data before now, so there is a problem
               CALL SetErrStat( ErrID_Fatal, ' Found text on line '//TRIM(Num2LStr(LineNumber))//' of '//TRIM(DataFileName)// &
                           ' when real numbers were expected.  There may be a problem with format of the file: '// &
                           TRIM(DataFileName)//'.', ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  RETURN
               ENDIF
            ELSE
               NumHeaderLines = NumHeaderLines + 1
            ENDIF
         ELSE     ! No text, must be data line
            NumDataLines = NumDataLines + 1
               ! If this is the first row of data, then store the number of words that were on the line
            IF ( .NOT. HaveReadData )  THEN
                  ! If this is the first line of data, keep some relevant info about it and the number of columns in it
               HaveReadData      = .TRUE.
               FirstDataLineNum  = LineNumber         ! Keep the line number of the first row of data (for error reporting)
               NumDataColumns    = NumWords
            ELSE
                  ! Make sure that the number columns on the row matches the number of columnns on the first row of data.
               IF ( NumWords /= NumDataColumns ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error in file: '//TRIM(DataFileName)//'.'// &
                           ' The number of data columns on line '//TRIM(Num2LStr(LineNumber))// &
                           '('//TRIM(Num2LStr(NumWords))//' columns) is different than the number of columns on first row of data '// &
                           ' (line: '//TRIM(Num2LStr(FirstDataLineNum))//', '//TRIM(Num2LStr(NumDataColumns))//' columns).', &
                           ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev ) THEN
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      ENDDO


      REWIND( UnitDataFile )

   END SUBROUTINE GetFileLength

   !-------------------------------------------------------------------------------
   !> This subroutine takes a line of text that is passed in and reads the first
   !! word to see if it is a number.  An internal read is used to do this.  If
   !! it is a number, it is started in ValueRead and returned. The flag IsRealNum
   !! is set to true.  Otherwise, ValueRead is set to NaN (value from the NWTC_Num)
   !! and the flag is set to false.
   !!
   !! The IsRealNum flag is set to indicate if we actually have a real number or
   !! not.  After calling this routine, a simple if statement can be used:
   !!
   !!       @code
   !!    IF (IsRealNum) THEN
   !!       ! do something
   !!    ELSE
   !!       ! do something else
   !!    ENDIF
   !!       @endcode
   !!
   !-------------------------------------------------------------------------------
   SUBROUTINE ReadRealNumberFromString(StringToParse, ValueRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)

      CHARACTER(*),        INTENT(IN   )           :: StringToParse  !< The string we were handed.
      REAL(ReKi),          INTENT(  OUT)           :: ValueRead      !< The variable being read.  Returns as NaN (library defined) if not a Real.
      CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
      LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
      INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
      CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
      INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.



         ! Initialize some things
      ErrStat     = ErrID_None
      ErrMsg      = ''


         ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
      read(StringToParse,*,IOSTAT=IOErrStat)   StrRead
      read(StringToParse,*,IOSTAT=IOErrStat)   ValueRead


         ! If IOErrStat==0, then we have a real number, anything else is a problem.
      if (IOErrStat==0) then
         IsRealNum   = .TRUE.
      else
         IsRealNum   = .FALSE.
         ValueRead   = NaN                ! This is NaN as defined in the NWTC_Num.
         ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
         ErrSTat     = ErrID_Severe
      endif



      RETURN
   END SUBROUTINE ReadRealNumberFromString


   !-------------------------------------------------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
   !> This subroutine works with the ReadNum routine from the library.  ReadNum is
   !! called to read a word from the input file.  An internal read is then done to
   !! convert the string to a number that is stored in VarRead and returned.
   !!
   !! The IsRealNum flag is set to indicate if we actually have a real number or
   !! not.  After calling this routine, a simple if statement can be used:
   !!
   !!       @code
   !!    IF (ISRealNum) THEN
   !!       ! do something
   !!    ELSE
   !!       ! do something else
   !!    ENDIF
   !!       @endcode
   !!
   !-------------------------------------------------------------------------------
   SUBROUTINE ReadRealNumber(UnitNum, FileName, VarName, VarRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)

      INTEGER(IntKi),      INTENT(IN   )           :: UnitNum        !< The unit number of the file being read
      CHARACTER(*),        INTENT(IN   )           :: FileName       !< The name of the file being read.  Used in the ErrMsg from ReadNum (Library routine).
      CHARACTER(*),        INTENT(IN   )           :: VarName        !< The variable we are reading.  Used in the ErrMsg from ReadNum (Library routine)'.
      REAL(ReKi),          INTENT(  OUT)           :: VarRead        !< The variable being read.  Returns as NaN (library defined) if not a Real.
      CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
      LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
      INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
      CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
      INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.

         ! Local vars
      INTEGER(IntKi)                      :: ErrStatTmp
      CHARACTER(2048)                     :: ErrMsgTmp



         ! Initialize some things
      ErrStat     = ErrID_None
      ErrMsg      = ''


         ! Now call the ReadNum routine to get the number
         ! If it is a word that does not start with T or F, then ReadNum won't give any errors.
      CALL ReadNum( UnitNum, FileName, StrRead, VarName, ErrStatTmp, ErrMsgTmp)


         ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
      read(StrRead,*,IOSTAT=IOErrStat)   VarRead


         ! If IOErrStat==0, then we have a real number, anything else is a problem.
      if (IOErrStat==0) then
         IsRealNum   = .TRUE.
      else
         IsRealNum   = .FALSE.
         VarRead     = NaN                ! This is NaN as defined in the NWTC_Num.
         ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
         ErrStat     = ErrStatTmp         ! The ErrStatTmp returned by the ReadNum routine is an ErrID level.
      endif



      RETURN
   END SUBROUTINE ReadRealNumber


END SUBROUTINE ReadPointsFile



SUBROUTINE AddedMassMessage( AM, ToFile, Msg, MsgLen )

   REAL(ReKi),                         INTENT(IN   )  :: AM(6,6)        !< Added mass matrix
   LOGICAL,                            INTENT(IN   )  :: ToFile         !< Prepend comment character
   CHARACTER(2048),                    INTENT(  OUT)  :: Msg
   INTEGER(IntKi),                     INTENT(  OUT)  :: MsgLen

      ! Local Variables
   CHARACTER(15)                                      :: TmpNumString
   INTEGER(IntKi)                                     :: ErrStatTmp
   INTEGER(IntKi)                                     :: I              !< Simple counter

   Msg   =  ''

   IF ( ToFile ) THEN
      Msg='#  Added Mass Values (kg, kg-m, kg-m^2):'//NewLine//NewLine
   ELSE
      Msg="Added Mass values (kg, kg-m, kg-m^2):"//NewLine//NewLine
   ENDIF

   ! Header info:
   Msg   =  TRIM(Msg)
   IF ( ToFile ) Msg=TRIM(Msg)//'#'
   Msg   =           TRIM(Msg)//"  Dim      TDxi        TDyi        TDzi        RDxi        RDyi        RDzi    "//NewLine
   IF ( ToFile ) Msg=TRIM(Msg)//'#'
   Msg   =           TRIM(Msg)//" ------------------------------------------------------------------------------"//NewLine

   MsgLen=  LEN_TRIM(Msg)-1         ! Not sure why an extra count exists here.


   CALL printDirection(" TDxi",1)
   CALL printDirection(" TDyi",2)
   CALL printDirection(" TDzi",3)
   CALL printDirection(" RDxi",4)
   CALL printDirection(" RDyi",5)
   CALL printDirection(" RDzi",6)



   RETURN

   CONTAINS
   SUBROUTINE printDirection( NameIn, IndexNum)

      CHARACTER(*),                    INTENT(IN   )  :: NameIn
      INTEGER(IntKi),                  INTENT(IN   )  :: IndexNum

      IF ( ToFile ) THEN
         Msg=  TRIM(Msg)//"#"//TRIM(NameIn)
      ELSE
         Msg=  TRIM(Msg)//" "//TRIM(NameIn)
      ENDIF
      MsgLen=  MsgLen+8
      Msg   =  Msg(1:MsgLen)//" "
      MsgLen=  MsgLen+1
      DO I=1,6
         WRITE(TmpNumString,'(ES10.3E2)',IOSTAT=ErrStatTmp)           AM(IndexNum,I)
         Msg   =  Msg(1:MsgLen)//TmpNumString(1:10)
         MsgLen=  MsgLen+2+10
      ENDDO
      Msg   =  Msg(1:MsgLen)//NewLine
   
   

      RETURN

   END SUBROUTINE printDirection

END SUBROUTINE AddedMassMessage


!> This subroutine writes the Added Mass matrix to a file
SUBROUTINE AddedMass_OutputWrite (DvrSettings, Initialized, PtfmAM, ErrStat, ErrMsg)

   TYPE( OrcaDriver_Settings ),        INTENT(INOUT)  :: DvrSettings          !< Stored settings
   LOGICAL,                            INTENT(INOUT)  :: Initialized          !< Was this file started before?
   REAL(ReKi),                         INTENT(IN   )  :: PtfmAM(6,6)          !< The added mass matrix
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None

         ! Temporary local variables
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary variable for the status of error message
   CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary variable for the error message
   CHARACTER(*),              PARAMETER               :: RoutineName = 'AddedMass_OutputWrite'
   INTEGER(IntKi)                                     :: LenErrMsgTmp         !< Length of ErrMsgTmp (for getting WindGrid info)

   CHARACTER(25)                                      :: AMfmt                !< Format specifier for the output file for wave elevation series
   INTEGER(IntKi)                                     :: I                    !< generic counter


   AMfmt = "(ES10.3E2,5(3x,ES10.3E2))"

   ErrMsg      = ''
   ErrStat     = ErrID_None
   ErrMsgTmp   = ''
   ErrStatTmp  = ErrID_None


      ! If it hasn't been initially written to, do this then exit. Otherwise set a few things and continue.
   IF ( .NOT. Initialized ) THEN

      CALL GetNewUnit( DvrSettings%AddedMassOutputUnit )
      CALL OpenFOutFile( DvrSettings%AddedMassOutputUnit, TRIM(DvrSettings%AddedMassFileName), ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      Initialized =  .TRUE.

         ! Write header section
      WRITE( DvrSettings%AddedMassOutputUnit,'(A)', IOSTAT=ErrStatTmp )   '## This file was generated by '//TRIM(GetNVD(DvrSettings%ProgInfo))//  &
            ' on '//CurDate()//' at '//CurTime()//'.'
      WRITE( DvrSettings%AddedMassOutputUnit,'(A)', IOSTAT=ErrStatTmp )   '## This file contains the added mass matrix that is returned from OrcaFlex'
      WRITE (DvrSettings%AddedMassOutputUnit,'(A)', IOSTAT=ErrStatTmp  )  '## It is arranged in a 6x6 matrix'
      WRITE (DvrSettings%AddedMassOutputUnit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      CALL AddedMassMessage( PtfmAM, .FALSE., ErrMsgTmp, LenErrMsgTmp )
      WRITE (DvrSettings%AddedMassOutputUnit,'(A)', IOSTAT=ErrStatTmp  )  ErrMsgTmp(1:LenErrMsgTmp)
   ELSE
      ! keep this as a placeholder in case we decide to write out at each timestep.
   ENDIF


END SUBROUTINE AddedMass_OutputWrite


SUBROUTINE PointsForce_OutputWrite(ProgInfo, OutUnit, OutFileName, InputFileName, Initialized, AnglesInDegrees, TotalPoints, &
                     Time, InitOutData, p, u, y, ErrStat, ErrMsg)

   TYPE(ProgDesc),                     INTENT(IN   )  :: ProgInfo             !< Program info
   INTEGER(IntKi),                     INTENT(INOUT)  :: OutUnit              !< Output Unit number
   CHARACTER(1024),                    INTENT(IN   )  :: OutFileName          !< Name of the file to write to
   CHARACTER(1024),                    INTENT(IN   )  :: InputFileName        !< Name of the file the points came from
   LOGICAL,                            INTENT(INOUT)  :: Initialized          !< Is the file initialized
   LOGICAL,                            INTENT(IN   )  :: AnglesInDegrees      !< The angles are in degrees.
   INTEGER(IntKi),                     INTENT(IN   )  :: TotalPoints          !< The total number of points in the points file
   REAL(DbKi),                         INTENT(IN   )  :: Time                 !< Current time
   TYPE(Orca_InitOutputType),          INTENT(IN   )  :: InitOutData          !< InitOutData -- need the header info
   TYPE(Orca_ParameterType),           INTENT(IN   )  :: p                    !< p
   TYPE(Orca_InputType),               INTENT(IN   )  :: u                    !< u
   TYPE(Orca_OutputType),              INTENT(IN   )  :: y                    !< y
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None

         ! Temporary local variables
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary variable for the status of error message
   CHARACTER(*),                       PARAMETER      :: RoutineName = 'PointsForce_OutputWrite'
   CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary variable for the error message
   INTEGER(IntKi)                                     :: LenErrMsgTmp         !< Length of ErrMsgTmp (for getting WindGrid info)
   REAL(ReKi)                                         :: rotdisp(3)           !< Rotational displacement (euler angles)
   INTEGER(IntKi)                                     :: I                    !< Generic counter
   REAL(ReKi)                                         :: outputArray(13)
   
   CHARACTER(47)                                      :: PointsOutputFmt      !< Format specifier for the output file for wave elevation series
   CHARACTER(3)                                       :: AngleUnit            !< Units for the angle

   PointsOutputFmt = "(ES10.3E2,18(3x,ES10.3E2))"

   ErrMsg      = ''
   ErrStat     = ErrID_None
   ErrMsgTmp   = ''
   ErrStatTmp  = ErrID_None


      ! If it hasn't been initially written to, do this then exit. Otherwise set a few things and continue.
   IF ( .NOT. Initialized ) THEN

      IF ( AnglesInDegrees ) THEN
         AngleUnit   =  "deg"
      ELSE
         AngleUnit   =  "rad"
      ENDIF

      CALL GetNewUnit( OutUnit )
      CALL OpenFOutFile( OutUnit, TRIM(OutFileName), ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      Initialized =  .TRUE.

         ! Write header section
      WRITE( OutUnit,'(A)', IOSTAT=ErrStatTmp )   '## This file was generated by '//TRIM(GetNVD(ProgInfo))//  &
            ' on '//CurDate()//' at '//CurTime()//'.'
      IF ( TotalPoints >= 1_IntKi ) THEN
         WRITE( OutUnit,'(A)', IOSTAT=ErrStatTmp )   '## This file contains the resulting forces and moments for the '//   &
                                                      TRIM(Num2LStr(TotalPoints))//' points specified in the '// &
                                                      'file '//TRIM(InputFileName)//'.'
      ENDIF
      WRITE (OutUnit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      CALL WrFileNR( OutUnit,                      '#   Time     '//  &
                                                   '   TDxi         TDyi         TDzi      '       //  &
                                                   '   RDxi         RDyi         RDzi      '       //  &
                                                   '   TVxi         TVyi         TVzi      '       //  &
                                                   '   RVxi         RVyi         RVzi      '          )
      DO I=1,SIZE(InitOutData%WriteOutputHdr)
         CALL WrFileNR ( OutUnit, '   '//InitOutData%WriteOutputHdr(I) )
      ENDDO ! I
      WRITE (OutUnit,'(A)', IOSTAT=ErrStatTmp ) ''



      CALL WrFileNR( OutUnit,                      '#   (s)      '//  &
                                                   '   (m)          (m)          (m)       '       //  &
                                                   '   ('//AngleUnit//')        ('//AngleUnit//')        ('//AngleUnit//')     '//  &
                                                   '   (m/s)        (m/s)        (m/s)     '       //  &
                                                   '   ('//AngleUnit//'/s)      ('//AngleUnit//'/s)      ('//AngleUnit//'/s)   ' )
      DO I=1,SIZE(InitOutData%WriteOutputHdr)
         CALL WrFileNR ( OutUnit, '   '//InitOutData%WriteOutputUnt(I) )
      ENDDO ! I
      WRITE (OutUnit,'(A)', IOSTAT=ErrStatTmp ) ''
   ENDIF

   rotdisp = GetSmllRotAngs ( u%PtfmMesh%Orientation(:,:,1), ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )

   IF ( AnglesInDegrees ) THEN
      outputArray = (/ REAL(Time, ReKi),                                                                                     &
                    REAL(u%PtfmMesh%TranslationDisp(1,1), ReKi),                                                             &
                    REAL(u%PtfmMesh%TranslationDisp(2,1), ReKi),                                                             &
                    REAL(u%PtfmMesh%TranslationDisp(3,1), ReKi),                                                             &
                    rotdisp(1)*R2D,                  rotdisp(2)*R2D,                  rotdisp(3)*R2D,                        &
                    u%PtfmMesh%TranslationVel(1,1),  u%PtfmMesh%TranslationVel(2,1),  u%PtfmMesh%TranslationVel(3,1),        &
                    u%PtfmMesh%RotationVel(1,1)*R2D, u%PtfmMesh%RotationVel(2,1)*R2D, u%PtfmMesh%RotationVel(3,1)*R2D  /)
   ELSE
      outputArray = (/ REAL(Time, ReKi),                                                                                     &
                    REAL(u%PtfmMesh%TranslationDisp(1,1), ReKi),                                                             &
                    REAL(u%PtfmMesh%TranslationDisp(2,1), ReKi),                                                             &
                    REAL(u%PtfmMesh%TranslationDisp(3,1), ReKi),                                                             &
                    rotdisp(1),                      rotdisp(2),                      rotdisp(3),                            &
                    u%PtfmMesh%TranslationVel(1,1),  u%PtfmMesh%TranslationVel(2,1),  u%PtfmMesh%TranslationVel(3,1),        &
                    u%PtfmMesh%RotationVel(1,1),     u%PtfmMesh%RotationVel(2,1),     u%PtfmMesh%RotationVel(3,1)      /)
   ENDIF
   
   CALL WrNumAryFileNR( OutUnit, outputArray, '3x,ES10.3E2', ErrStatTmp, ErrMsgTmp )
   CALL WrNumAryFileNR( OutUnit, y%WriteOutput, '3x,ES10.3E2', ErrStatTmp, ErrMsgTmp )
   WRITE (OutUnit,'(A)', IOSTAT=ErrStatTmp ) ''

END SUBROUTINE PointsForce_OutputWrite




!> This routine exists only to support the development of the module.  It will not be needed after the module is complete.
SUBROUTINE  printSettings( DvrFlags, DvrSettings )
      ! The arguments
   TYPE( OrcaDriver_Flags ),            INTENT(IN   )  :: DvrFlags           !< Flags indicating which settings were set
   TYPE( OrcaDriver_Settings ),         INTENT(IN   )  :: DvrSettings        !< Stored settings

   CALL WrsCr(TRIM(GetNVD(DvrSettings%ProgInfo)))
   CALL WrScr(' DvrIptFile:          '//FLAG(DvrFlags%DvrIptFile)//        '      '//TRIM(DvrSettings%DvrIptFileName))
   CALL WrScr(' OrcaIptFile:         '//FLAG(DvrFlags%OrcaIptFile)//       '      '//TRIM(DvrSettings%OrcaIptFileName))
!   CALL WrScr(' DLLPathFileName:     '//FLAG(DvrFlags%DLLPathFileName)//   '      '//TRIM(DvrSettings%DLLPathFileName))
   CALL WrScr(' PointsFile:          '//FLAG(DvrFlags%PointsFile)//        '      '//TRIM(DvrSettings%PointsFileName))
   CALL WrScr(' AddedMass:           '//FLAG(DvrFlags%AddedMass))
   CALL WrScr(' AddedMassFile:       '//FLAG(DvrFlags%AddedMassFile)//     '      '//TRIM(DvrSettings%AddedMassFileName))
   IF ( DvrFlags%DTDefault) THEN
      CALL WrScr(' DT:                  '//FLAG(DvrFlags%DT)//                '      DEFAULT')
   ELSE
      CALL WrScr(' DT:                  '//FLAG(DvrFlags%DT)//                '      '//TRIM(Num2LStr(DvrSettings%DT)))
   ENDIF
   IF ( DvrFlags%Degrees ) THEN
      CALL WrScr(' PtfmCoord:           '//FLAG(DvrFlags%PtfmCoord)//         '     ['//TRIM(Num2LStr(DvrSettings%PtfmCoord(1)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(2)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(3)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(4)*R2D))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(5)*R2D))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(6)*R2D))//']')
      CALL WrScr(' PtfmVeloc:           '//FLAG(DvrFlags%PtfmVeloc)//         '     ['//TRIM(Num2LStr(DvrSettings%PtfmVeloc(1)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(2)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(3)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(4)*R2D))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(5)*R2D))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(6)*R2D))//']')
      CALL WrScr(' PtfmAccel:           '//FLAG(DvrFlags%PtfmAccel)//         '     ['//TRIM(Num2LStr(DvrSettings%PtfmAccel(1)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(2)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(3)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(4)*R2D))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(5)*R2D))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(6)*R2D))//']')
   ELSE
      CALL WrScr(' PtfmCoord:           '//FLAG(DvrFlags%PtfmCoord)//         '     ['//TRIM(Num2LStr(DvrSettings%PtfmCoord(1)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(2)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(3)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(4)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(5)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmCoord(6)))//']')
      CALL WrScr(' PtfmVeloc:           '//FLAG(DvrFlags%PtfmVeloc)//         '     ['//TRIM(Num2LStr(DvrSettings%PtfmVeloc(1)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(2)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(3)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(4)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(5)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmVeloc(6)))//']')
      CALL WrScr(' PtfmAccel:           '//FLAG(DvrFlags%PtfmAccel)//         '     ['//TRIM(Num2LStr(DvrSettings%PtfmAccel(1)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(2)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(3)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(4)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(5)))//', '&
                                                                                      //TRIM(Num2LStr(DvrSettings%PtfmAccel(6)))//']')
   ENDIF
   CALL WrScr(' Degrees:             '//FLAG(DvrFlags%Degrees)//           '     PtfmCoord, PtfmVeloc, and PtfmAccel angles in degrees')
   CALL WrScr(' PointsDegrees:       '//FLAG(DvrFlags%PointsDegrees)//     '     PointsFile angles in degrees')
   CALL WrScr(' PointsOutputInit:    '//FLAG(DvrFlags%PointsOutputInit)//  '      Unit #:  '//TRIM(Num2LStr(DvrSettings%PointsOutputUnit)))
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
END MODULE OrcaDriver_Subs
