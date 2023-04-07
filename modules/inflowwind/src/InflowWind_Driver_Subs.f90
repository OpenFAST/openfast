!**********************************************************************************************************************************
!
!  MODULE: InflowWind_Driver_Subs  - This module contains subroutines used by the InflowWind Driver program
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2017 Envision Energy USA
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
MODULE InflowWind_Driver_Subs

   USE NWTC_Library
   USE InflowWind_Driver_Types
   USE InflowWind_IO
   USE IfW_FlowField
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
   CALL WrScr("     options:     "//SWChar//"ifw           -- treat <filename> as name of InflowWind input file")
   CALL WrScr("                                    (no driver input file)")
   CALL WrScr("")
   CALL WrScr("              The following options will overwrite values in the driver input file:")
   CALL WrScr("                  "//SwChar//"DT[#]          -- timestep                                        ")
   CALL WrScr("                  "//SwChar//"TStart[#]      -- start time                                      ")
   CALL WrScr("                  "//SwChar//"TSteps[#]      -- number of timesteps                             ")
   CALL WrScr("                  "//SwChar//"xrange[#:#]    -- range of x (#'s are reals)                      ")
   CALL WrScr("                  "//SwChar//"yrange[#:#]    -- range of y                                      ")
   CALL WrScr("                  "//SwChar//"zrange[#:#]    -- range in z (ground = 0.0)                       ")
   CALL WrScr("                  "//SwChar//"Dx[#]          -- spacing in x                                    ")
   CALL WrScr("                  "//SwChar//"Dy[#]          -- spacing in y                                    ")
   CALL WrScr("                  "//SwChar//"Dz[#]          -- spacing in z                                    ")
!   CALL WrScr("                  "//SwChar//"sum            -- summarize wind file info                        [N/A]")
!   CALL WrScr("                  "//SwChar//"FFT[X,Y,Z]     -- an fft over all t using specified DT at X,Y,Z   [N/A]")
   CALL WrScr("                  "//SwChar//"points[FILE]   -- calculates at x,y,z coordinates specified in a  ")
   CALL WrScr("                                    white space delimited FILE")
   CALL WrScr("                  "//SwChar//"v              -- verbose output ")
   CALL WrScr("                  "//SwChar//"vv             -- very verbose output ")
   CALL WrScr("                  "//SwChar//"HAWC           -- convert contents of <filename> to HAWC format ")
   CALL WrScr("                  "//SwChar//"Bladed         -- convert contents of <filename> to Bladed format ")
   CALL WrScr("                  "//SwChar//"vtk            -- convert contents of <filename> to vtk format ")
   CALL WrScr("                  "//SwChar//"accel          -- calculate wind acceleration in addition to velocity")
   CALL WrScr("                  "//SwChar//"BoxExceedAllow -- set flag to allow FF points outside wind box")
   CALL WrScr("                  "//SwChar//"help           -- print this help menu and exit")
   CALL WrScr("")
   CALL WrScr("   Notes:")
   CALL WrScr("   -- Unspecified ranges and resolutions default to what is in the file.")
   CALL WrScr("   -- If no XRange is specified, assumed to be only at X=0")
!   CALL WrScr("   -- Features marked [N/A] have not been implimented in this version.")
   CALL WrScr("   -- Options are not case sensitive.")
   CALL WrScr("")


END SUBROUTINE DispHelpText


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine retrieves the command line arguments and passes them to the
!! inflowwind_driver_subs::parsearg routine for processing.
SUBROUTINE RetrieveArgs( CLSettings, CLFlags, ErrStat, ErrMsg )

   USE NWTC_Library
   USE InflowWind_Driver_Types

   IMPLICIT NONE

      ! Storing the arguments
   TYPE( IfWDriver_Flags ),            INTENT(  OUT)  :: CLFlags      !< Flags indicating which command line arguments were specified
   TYPE( IfWDriver_Settings ),         INTENT(  OUT)  :: CLSettings           !< Command line arguments passed in

      ! Error Handling
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

      ! Local variable     IF ( CLFlags%Summarys
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
   CLFlags%DvrIptFile =  .FALSE.
   ErrStat                    =  ErrID_None
   ErrStatTmp                 =  ErrID_None
   ErrMsg                     =  ''
   ErrMsgTmp                  =  ''
   ifwFlag                    =  .FALSE.
   FileNameGiven              =  .FALSE.
   FileName                   =  ''


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
         CALL ParseArg( CLSettings, CLFlags, ArgUC(2:), Arg(2:), ifwFlag, ErrStatTmp, ErrMsgTmp )
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
      CLFlags%IfWIptFile =  .TRUE.
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
      USE InflowWind_Driver_Types
      USE InflowWind_Types

      IMPLICIT NONE

         ! Storing the arguments
      TYPE( IfWDriver_Flags ),            INTENT(INOUT)  :: CLFlags              ! Flags indicating which arguments were specified
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
      ErrMsg      = ''

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
         IF       ( TRIM(ThisArgUC) == "IFW" )   THEN
            ifwFlagSet              = .TRUE.             ! More logic in the routine that calls this one to set things.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "SUM" )   THEN
            CLFlags%Summary         = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "VV"  )   THEN
            CLFlags%VVerbose        = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "V"   )   THEN
            CLFlags%Verbose         = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "HAWC"   )   THEN
            CLFlags%WrHAWC       = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "BLADED"   )   THEN
            CLFlags%WrBladed     = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "VTK"   )   THEN
            CLFlags%WrVTK        = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "UNIFORM"   )   THEN
            CLFlags%WrUniform    = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "BOXEXCEEDALLOW"   )   THEN
            CLFlags%BoxExceedAllowF = .TRUE.
            RETURN
         ELSEIF   ( TRIM(ThisArgUC) == "ACCEL"   )   THEN
            CLFlags%OutputAccel    = .TRUE.
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
            CLFlags%FFTcalc            = .TRUE.
            CLSettings%FFTcoord(1)     = TempReal
         ELSE
            CLFlags%FFTcalc            = .FALSE.
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
            CLFlags%FFTcalc            = .TRUE.
            CLSettings%FFTcoord(2)     = TempReal
         ELSE
            CLFlags%FFTcalc            = .FALSE.
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
            CLFlags%FFTcalc            = .TRUE.
            CLSettings%FFTcoord(3)     = TempReal
         ELSE
            CLFlags%FFTcalc            = .FALSE.
            IF ( ErrStatTmp == ErrID_Warn ) THEN
               CALL SetErrStat( ErrStatTmp," Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ELSE
               CALL SetErrStat( ErrID_Fatal," Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'.", &
                  ErrStat, ErrMsg, 'ParseArg')
            ENDIF
            RETURN
         ENDIF


         IF ( CLSettings%FFTcoord(3) <= 0.0_ReKi ) THEN
            CALL SetErrStat( ErrID_Warn,' FFT coordinate ['//TRIM(Num2LStr(CLSettings%FFTcoord(1)))//','//  &
               TRIM(Num2LStr(CLSettings%FFTcoord(1)))//','//TRIM(Num2LStr(CLSettings%FFTcoord(1)))//        &
               '] is at or below ground level where there is no wind.  Ingoring.',   &
               ErrStat,ErrMsg,'UpdateSettingsWithCL' )
            CLFlags%FFTcalc =  .FALSE.
         ENDIF

          ! "XRANGE[#:#]"
      ELSEIF ( ThisArgUC(1:Delim1) == "XRANGE["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF     ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%XRange = .TRUE.
            CLSettings%XRange(1)    = TempReal
         ELSE
            CLFlags%XRange          = .FALSE.
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
            CLFlags%XRange          = .TRUE.
            CLSettings%XRange(2)    = TempReal
         ELSE
            CLFlags%XRange          = .FALSE.
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
            CLSettings%XRange(1)    = 0.0
            CLSettings%XRange(2)    = 0.0
            CLFlags%XRange          = .FALSE.
            CALL SetErrStat(ErrID_Warn," Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring.",   &
               ErrStat, ErrMsg, 'ParseArgs')
         ENDIF



         ! "YRANGE[#:#]"
      ELSEIF ( ThisArgUC(1:Delim1) == "YRANGE["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF     ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%YRange          = .TRUE.
            CLSettings%YRange(1)    = TempReal
         ELSE
            CLFlags%YRange          = .FALSE.
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
            CLFlags%YRange          = .TRUE.
            CLSettings%YRange(2)    = TempReal
         ELSE
            CLFlags%YRange          = .FALSE.
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
            CLSettings%YRange(1)    = 0.0
            CLSettings%YRange(2)    = 0.0
            CLFlags%YRange          = .FALSE.
            CALL SetErrStat(ErrID_Warn," Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring.",   &
               ErrStat, ErrMsg, 'ParseArgs')
         ENDIF



         ! "ZRANGE[#:#]"
      ELSEIF ( ThisArgUC(1:Delim1) == "ZRANGE["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArgUC(Delim1+1:DelimSep-1), ErrStatTmp )
         IF     ( ErrStatTmp == ErrID_None ) THEN
            CLFlags%ZRange          = .TRUE.
            CLSettings%ZRange(1)    = TempReal
         ELSE
            CLFlags%ZRange          = .FALSE.
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
            CLFlags%ZRange          = .TRUE.
            CLSettings%ZRange(2)    = TempReal
         ELSE
            CLFlags%ZRange          = .FALSE.
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
            CLSettings%ZRange(1)    = 0.0
            CLSettings%ZRange(2)    = 0.0
            CLFlags%ZRange          = .FALSE.
            CALL SetErrStat(ErrID_Warn," Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring.",   &
               ErrStat, ErrMsg, 'ParseArgs')
         ENDIF


         ! "DX[#]"
      ELSEIF( ThisArgUC(1:Delim1) == "DX["      ) THEN
         TempReal = StringToReal( ThisArgUC(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_None ) THEN
            CLFlags%Dx           = .TRUE.
            CLSettings%GridDelta(1) = abs(TempReal)
         ELSE
            CLFlags%Dx           = .FALSE.
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
            CLFlags%Dy           = .TRUE.
            CLSettings%GridDelta(2) = abs(TempReal)
         ELSE
            CLFlags%Dy           = .FALSE.
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
            CLFlags%Dz           = .TRUE.
            CLSettings%GridDelta(3) = abs(TempReal)
         ELSE
            CLFlags%Dz           = .FALSE.
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
            CLFlags%Dt           = .TRUE.
            CLSettings%DT        = abs(TempReal)
         ELSE
            CLFlags%Dt           = .FALSE.
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
            CLFlags%NumTimeSteps  = .TRUE.
            CLSettings%NumTimeSteps       = nint(abs(TempReal))
         ELSE
            CLFlags%NumTimeSteps  = .FALSE.
            CLSettings%NumTimeSteps = 1_IntKi
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
            CLFlags%TStart          = .TRUE.
            CLSettings%TStart       = abs(TempReal)
         ELSE
            CLFlags%TStart          = .FALSE.
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
   INTEGER(IntKi)                                     :: TmpIntAr3(3)         ! Temporary array for reading in a pair of integer values from the input file
   REAL(ReKi)                                         :: TmpRealAr3(3)        ! Temporary array for reading in a pair of real values from the input file
   REAL(ReKi)                                         :: GridCtrCoord(3)      ! Center coordinate of the grid read in


      ! Local error handling
   INTEGER(IntKi)                                     :: ios                  !< I/O status
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error messages for calls
   CHARACTER(*), PARAMETER                            :: RoutineName = 'ReadDvrIptFile'
   CHARACTER(1024)                                    :: PriPath                                   ! Path name of the primary file

      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1

   FileName = TRIM(DvrFileName)

   CALL GetNewUnit( UnIn )
   CALL OpenFInpFile( UnIn, FileName, ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,' Failed to open InflowWind Driver input file: '//FileName,   &
         ErrStat,ErrMsg,RoutineName)
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL WrScr( 'Opening InflowWind Driver input file:  '//FileName )

   CALL GetPath( FileName, PriPath )    ! Input files will be relative to the path where the primary input file is located.

   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, FileName,' InflowWind Driver input file header line 1', ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CLOSE( UnIn )
      RETURN
   ENDIF


   CALL ReadCom( UnIn, FileName, 'InflowWind Driver input file header line 2', ErrStatTmp, ErrMsgTmp )
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
      CALL ReadCom( UnIn, FileName,' InflowWind Driver input file header line 1', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF


      CALL ReadCom( UnIn, FileName, 'InflowWind Driver input file header line 2', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
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
   !  Driver setup section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' Driver setup section, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


      ! Name of InflowWind input file
   CALL ReadVar( UnIn, FileName,DvrSettings%IfWIptFileName,'IfWIptFileName',' InflowWind input filename',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%IfWIptFile  =  .TRUE.
   ENDIF
   
   IF ( PathIsRelative( DvrSettings%IfWIptFileName ) ) DvrSettings%IfWIptFileName = TRIM(PriPath)//TRIM(DvrSettings%IfWIptFileName)

   !-------------------------------------------------------------------------------------------------
   !  File Conversion Options section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,'File Conversion Options Section Header', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)

      ! WrHAWC
   CALL ReadVar( UnIn, FileName, DvrFlags%WrHAWC, 'WrHAWC', 'Convert wind data to HAWC2 format?', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      
      ! WrBladed
   CALL ReadVar( UnIn, FileName, DvrFlags%WrBladed, 'WrBladed', 'Convert wind data to Bladed format?', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)

      ! WrVTK
   CALL ReadVar( UnIn, FileName, DvrFlags%WrVTK, 'WrVTK', 'Convert wind data to VTK format?', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)

      ! WrUniform
   CALL ReadVar( UnIn, FileName, DvrFlags%WrUniform, 'WrUniform', 'Convert wind data to Uniform Wind format?', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


   !-------------------------------------------------------------------------------------------------
   !  Tests of Interpolation Options section
   !-------------------------------------------------------------------------------------------------
   CALL ReadCom( UnIn, FileName,'Tests of Interpolation Options Section Header', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      CALL SetErrStat(ErrStatTmp, ErrMsgTmp,ErrStat,ErrMsg,RoutineName)

   
      ! Number of timesteps
   CALL ReadVar( UnIn, FileName,NumTimeStepsChr,'NumTimeStepsChr',' Character string for number of timesteps to read.',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
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
      READ (NumTimeStepsChr,*,IOSTAT=IOS)   DvrSettings%NumTimeSteps
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'NumTimeSteps',NumType, ErrStatTmp, ErrMsgTmp )
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
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ELSE
      DvrFlags%TStart   =  .TRUE.
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
      READ (DTChr,*,IOSTAT=IOS)   DvrSettings%DT
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'DT',NumType, ErrStatTmp, ErrMsgTmp )
         RETURN
      ELSE     ! Was ok, so set the flags
         DvrFlags%DT         =  .TRUE.
         DvrFlags%DTDefault  =  .FALSE.
      ENDIF
   ENDIF
   
   
      ! Summarize the extents in the windfile
   CALL ReadVar( UnIn, FileName,DvrFlags%Summary,'Summary',' Summarize data extents in the windfile', &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
!   ELSE
!      DvrFlags%Summary  =  .TRUE.
   ENDIF


      ! Summarize everything in a summary file/
   CALL ReadVar( UnIn, FileName,DvrFlags%SummaryFile,'SummaryFile',' Summarize the results in a .sum file', &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
!   ELSE
!      DvrFlags%SummaryFile =  .TRUE.
   ENDIF

      ! Flag to allow sampling outside grid
   CALL ReadVar( UnIn, FileName,DvrFlags%BoxExceedAllowF,'BoxExceedAllow',' Allow point sampling outside grid', &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF

#ifdef UNUSED_INPUTFILE_LINES
   !-------------------------------------------------------------------------------------------------
   !  FFT calculations
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' FFT calculations, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


       ! FFTcalc    -- FFTcalc of the windfield needed.
   CALL ReadVar( UnIn, FileName,DvrFlags%FFTcalc,'FFTcalc',' Perform an FFT?',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
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
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
   ELSE
      CALL ReadCom( UnIn, FileName,' Skipping the FFT coordinate since not doint an FFT.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
        RETURN
      ENDIF
   ENDIF

#endif

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


      ! Points input file (unused if .not. DvrFlags%PointsFile)
   CALL ReadVar( UnIn, FileName,DvrSettings%PointsFileName,'PointsFileName',' Points file input filename',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF

   IF ( PathIsRelative( DvrSettings%PointsFileName ) ) DvrSettings%PointsFileName = TRIM(PriPath)//TRIM(DvrSettings%PointsFileName)

      ! CalcAccel - calculate wind acceleration (unused if .not. DvrFlags%PointsFile)
   CALL ReadVar( UnIn, FileName,DvrFlags%OutputAccel, 'CalcAccel', ' Calc and output wind acceleration',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF

   !-------------------------------------------------------------------------------------------------
   !  gridded data output
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName,' Gridded data output, comment line', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   ENDIF


       ! WindGrid    -- Gridded data output
   CALL ReadVar( UnIn, FileName,DvrFlags%WindGrid,'WindGrid',' Output a grid of data?',   &
      ErrStatTmp,ErrMsgTmp, UnEchoLocal )
   IF ( ErrStatTmp /= ErrID_None ) THEN
      CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
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
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF

         ! Read the DY and DZ stepsize
      CALL ReadAry ( UnIn, FileName, TmpRealAr3, 3, 'GridDX, GridDY, GridDZ', &
         'GridDX, GridDY, GridDZ', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF

         ! Save the DY and DZ values
      DvrSettings%GridDelta(1)   =  abs(TmpRealAr3(1))   ! X direction
      DvrSettings%GridDelta(2)   =  abs(TmpRealAr3(2))   ! Y direction
      DvrSettings%GridDelta(3)   =  abs(TmpRealAr3(3))   ! Z direction
      DvrFlags%Dx                =  .TRUE.               ! read in value for the X direction gridding
      DvrFlags%Dy                =  .TRUE.               ! read in value for the Y direction gridding
      DvrFlags%Dz                =  .TRUE.               ! read in value for the Z direction gridding


         ! Read the number of points in the Y and Z directions
      CALL ReadAry ( UnIn, FileName, TmpIntAr3, 3, 'GridNx, GridNY, GridNZ', &
         'GridNx, GridNY, GridNZ', ErrStatTmp, ErrMsgTmp, UnEchoLocal)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF

         ! Save the GridNY and GridNZ values
      DvrSettings%GridN(1)   =  TmpIntAr3(1)          ! X direction
      DvrSettings%GridN(2)   =  TmpIntAr3(2)          ! Y direction
      DvrSettings%GridN(3)   =  TmpIntAr3(3)          ! Z direction
      DvrFlags%XRange            =  .TRUE.               ! read in value for the X direction gridding
      DvrFlags%YRange            =  .TRUE.               ! read in value for the Y direction gridding
      DvrFlags%ZRange            =  .TRUE.               ! read in value for the Z direction gridding

      
      
         ! Check that valid values of Dx, Dy, and Dz were read in.
         ! Check GridDx
      IF ( EqualRealNos(DvrSettings%GridDelta(1), 0.0_ReKi) ) THEN
         DvrFlags%Dx                =  .FALSE.
         DvrFlags%XRange            =  .FALSE.
         DvrSettings%GridDelta(1)   = 0.0_ReKi
         CALL SetErrStat(ErrID_Warn,' Grid spacing in X direction is 0.  Ignoring.',ErrStat,ErrMsg,RoutineName)
      ENDIF

         ! Check GridDy
      IF ( EqualRealNos(DvrSettings%GridDelta(2), 0.0_ReKi) ) THEN
         DvrFlags%Dy                =  .FALSE.
         DvrFlags%YRange            =  .FALSE.
         DvrSettings%GridDelta(2)   = 0.0_ReKi
         CALL SetErrStat(ErrID_Warn,' Grid spacing in Y direction is 0.  Ignoring.',ErrStat,ErrMsg,RoutineName)
      ENDIF

         ! Check GridDz
      IF ( EqualRealNos(DvrSettings%GridDelta(3), 0.0_ReKi) ) THEN
         DvrFlags%Dz                =  .FALSE.
         DvrFlags%ZRange            =  .FALSE.
         DvrSettings%GridDelta(3)   = 0.0_ReKi
         CALL SetErrStat(ErrID_Warn,' Grid spacing in Z direction is 0.  Ignoring.',ErrStat,ErrMsg,RoutineName)
      ENDIF
      
      
         ! Now need to set the XRange, YRange, and ZRange values based on what we read in
         ! For XRange, check that we have an actual value for the number of points
      IF ( (DvrSettings%GridN(1) <= 0) .OR. (.NOT. DvrFlags%XRange) ) THEN
         DvrSettings%XRange   =  GridCtrCoord(1)
         DvrFlags%Dx          =  .FALSE.
         DvrFlags%XRange      =  .FALSE.

         IF ( DvrSettings%GridN(1) < 0 )  THEN
            CALL SetErrStat(ErrID_Warn,' Negative number for number of grid points along X direction.  Ignoring.',ErrStat,ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat(ErrID_Warn,' No points along X direction.  Ignoring.',ErrStat,ErrMsg, RoutineName)
         ENDIF

         DvrSettings%GridN(1) =  1_IntKi              ! Set to 1 for easier indexing.

      ELSE
            ! Set the XRange values
         DvrSettings%XRange(1)   =  GridCtrCoord(1) - (REAL(DvrSettings%GridN(1) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(1)
         DvrSettings%XRange(2)   =  GridCtrCoord(1) + (REAL(DvrSettings%GridN(1) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(1)
         DvrFlags%XRange         =  .TRUE.
      ENDIF


         ! For YRange, check that we have an actual value for the number of points
      IF ( (DvrSettings%GridN(2) <= 0) .OR. (.NOT. DvrFlags%YRange) ) THEN
         DvrSettings%YRange   =  GridCtrCoord(2)
         DvrFlags%Dy          =  .FALSE.
         DvrFlags%YRange      =  .FALSE.

         IF ( DvrSettings%GridN(2) < 0 )  THEN
            CALL SetErrStat(ErrID_Warn,' Negative number for number of grid points along Y direction.  Ignoring.',ErrStat,ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat(ErrID_Warn,' No points along Y direction.  Ignoring.',ErrStat,ErrMsg, RoutineName)
         ENDIF

         DvrSettings%GridN(2) =  1_IntKi              ! Set to 1 for easier indexing.

      ELSE
            ! Set the YRange values
         DvrSettings%YRange(1)   =  GridCtrCoord(2) - (REAL(DvrSettings%GridN(2) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(2)
         DvrSettings%YRange(2)   =  GridCtrCoord(2) + (REAL(DvrSettings%GridN(2) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(2)
         DvrFlags%YRange         =  .TRUE.
      ENDIF

         ! For ZRange, check that we have an actual value for the number of points, set to ctr point if negative or zero.
      IF ( (DvrSettings%GridN(3) <= 0) .OR. (.NOT. DvrFlags%ZRange) ) THEN
         DvrSettings%ZRange   =  abs(GridCtrCoord(3))       ! shouldn't have a negative value anyhow
         DvrFlags%Dz          =  .FALSE.
         DvrFlags%ZRange      =  .FALSE.

         IF ( DvrSettings%GridN(3) < 0 )  THEN
            CALL SetErrStat(ErrID_Warn,' Negative number for number of grid points along Z direction.  Ignoring.',ErrStat,ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat(ErrID_Warn,' No points along Z direction.  Ignoring.',ErrStat,ErrMsg, 'ReadDvrIptFile')
         ENDIF

         DvrSettings%GridN(3) =  1_IntKi              ! Set to 1 for easier indexing.

      ELSE
            ! Set the ZRange values
         DvrSettings%ZRange(1)   =  GridCtrCoord(3) - (REAL(DvrSettings%GridN(3) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(3)
         DvrSettings%ZRange(2)   =  GridCtrCoord(3) + (REAL(DvrSettings%GridN(3) - 1_IntKi ) / 2.0_ReKi ) * DvrSettings%GridDelta(3)
         DvrFlags%ZRange         =  .TRUE.
      ENDIF
   
      
   ELSE ! read these lines as comments (actually, we don't need to read them)
      
      
      DvrSettings%GridDelta = 0.0_ReKi
      DvrFlags%Dx                =  .FALSE.
      DvrFlags%Dy                =  .FALSE.
      DvrFlags%Dz                =  .FALSE.
      
      DvrSettings%GridN = 0.0_ReKi
      DvrFlags%XRange            =  .FALSE.
      DvrFlags%YRange            =  .FALSE.
      DvrFlags%ZRange            =  .FALSE.
      
         ! Skip the next three entries of the gridded data section.
      CALL ReadCom( UnIn, FileName,' Skipping the gridded data section since not calculating it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the gridded data section since not calculating it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         CALL CleanupEchoFile( EchoFileContents, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      ENDIF
      CALL ReadCom( UnIn, FileName,' Skipping the gridded data section since not calculating it.', ErrStatTmp, ErrMsgTmp, UnEchoLocal )
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
   ErrMsg      =  ''
   ErrStatTmp  =  ErrID_None
   ErrMsgTmp   =  ''


   DvrFlags%WrHAWC    = DvrFlags%WrHAWC    .or. CLFlags%WrHAWC             ! create file if specified in either place
   DvrFlags%WrBladed  = DvrFlags%WrBladed  .or. CLFlags%WrBladed           ! create file if specified in either place
   DvrFlags%WrVTK     = DvrFlags%WrVTK     .or. CLFlags%WrVTK              ! create file if specified in either place
   DvrFlags%WrUniform = DvrFlags%WrUniform .or. CLFlags%WrUniform          ! create file if specified in either place
   DvrFlags%BoxExceedAllowF   = DvrFlags%BoxExceedAllowF .or. CLFlags%BoxExceedAllowF  ! flag to allow points beyond box for FF
   DvrFlags%OutputAccel = DvrFlags%OutputAccel .or. CLFlags%OutputAccel          ! calculate acceleration if specified in either place

!      ! Due to the complexity, we are handling overwriting driver input file settings with
!      ! command line settings and the instance where no driver input file is read separately.
!   IF ( DVRIPT .AND. DvrFlags%WindGrid ) THEN

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
         DvrFlags%DT       =  .TRUE.
      ENDIF
      DvrSettings%DT       =  CLSettings%DT
      DvrFlags%DTDefault   =  .FALSE.
   ENDIF

      ! Check NumTimeSteps
   IF ( CLFlags%NumTimeSteps ) THEN
      IF ( DvrFlags%NumTimeSteps .AND. ( DvrSettings%NumTimeSteps /= CLSettings%NumTimeSteps ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for NumTimeSteps with '//  &
            TRIM(Num2LStr(CLSettings%NumTimeSteps))//'.',&
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ELSE
         DvrFlags%NumTimeSteps      =  .TRUE.
      ENDIF
      DvrSettings%NumTimeSteps      =  CLSettings%NumTimeSteps
      DvrFlags%NumTimeStepsDefault  =  .FALSE.
   ENDIF

      ! Make sure there is at least one timestep
   DvrSettings%NumTimeSteps   =  MAX(DvrSettings%NumTimeSteps,1_IntKi)


      !--------------------------------------------
      ! Did we change the ranges for the WindGrid?
      !--------------------------------------------

      ! XRange -- There are three bits of information required for the full specification of the Y
      !           data: Dx, XRange, and Ny (the number of steps).  In the driver input file, Dx and Ny are
      !           given, and XRange is calculated using the specified center point.  From the command line,
      !           Dx, and XRange are specified.
      !           In order for the command line values to be used, we must know which were specified or not
      !           as it will change the values we use for calculating.
      !
      !  First, check if WindGrid is set.  If not, set it.
      !
      !     XRange      Dx       What to change
      !        y        n        Assume Dx good.  Calculate new Ny.  Increase XRange(2) if necessary and warn.
      !        n        y        Assume XRange good (if Ny == 0, issue warning about no points, set to none)
      !        y        y        Calculate Nx from XRange and Dx, reset XRange(2) after calculating (extend
      !                          a bit if necessary).  Also set WindGrid to true if it is false.
      !        n        n        Nothing.
      !

   IF       ( CLFlags%XRange .AND. ( .NOT. CLFlags%Dx ) )   THEN

      DvrFlags%Dx       =  .TRUE.
      DvrFlags%XRange   =  .TRUE.
      WindGridModify    =  .TRUE.

         ! Calculate new Nx.
      DvrSettings%XRange(1)      =  minval(CLSettings%XRange)     ! just in case the order is funky
      DvrSettings%XRange(2)      =  maxval(CLSettings%XRange)

         ! Set number of points in X direction
      DvrSettings%GridN(1) =  CEILING( abs( ( DvrSettings%XRange(2) - DvrSettings%XRange(1) ) / DvrSettings%GridDelta(1) ) + 1_IntKi )

         ! Now adjust the upper limit of the range slightly for the integer number of steps
      DvrSettings%XRange(2)=  DvrSettings%XRange(1)   + ( DvrSettings%GridN(1) - 1_IntKi ) * DvrSettings%GridDelta(2)

         ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
      IF ( .NOT. EqualRealNos(DvrSettings%XRange(2),maxval(CLSettings%XRange)) )  THEN
         CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of XRange to '//TRIM(Num2LStr(DvrSettings%XRange(2)))//   &
            ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(1)))//' points with a requested spacing '//       &
            TRIM(Num2LStr(DvrSettings%GridDelta(1)))//'.',                                                           &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF


!         ! Only the XRange was specfied, so throw warning and don't change anything.
!      CALL SetErrStat( ErrID_Warn, ' XRange argument given, but no Dx.  Ignoring since no points to calculate.',  &
!         ErrStat,ErrMsg,'UpdateSettingsWithCL')
!      DvrFlags%XRange   =  .FALSE.     ! Should already be false, but just in case
!      DvrFlags%Dx       =  .FALSE.     ! Should already be false, but just in case

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

      DvrFlags%Dy       =  .TRUE.
      DvrFlags%YRange   =  .TRUE.
      WindGridModify    =  .TRUE.

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

         DvrFlags%Dy       =  .TRUE.
         DvrFlags%YRange   =  .TRUE.
         WindGridModify    =  .TRUE.

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
      DvrFlags%Dy       =  .TRUE.
      DvrFlags%YRange   =  .TRUE.
      WindGridModify    =  .TRUE.

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

   IF       ( CLFlags%ZRange .AND. ( .NOT. CLFlags%Dz ) )   THEN

      DvrFlags%Dz       =  .TRUE.
      DvrFlags%ZRange   =  .TRUE.
      WindGridModify    =  .TRUE.

         ! Calculate new Ny.
      DvrSettings%ZRange(1)      =  minval(CLSettings%ZRange)     ! just in case the order is funky
      DvrSettings%ZRange(2)      =  maxval(CLSettings%ZRange)

         ! Set number of points in Y direction
      DvrSettings%GridN(3) =  CEILING( abs( ( DvrSettings%ZRange(2) - DvrSettings%ZRange(1) ) / DvrSettings%GridDelta(3) ) + 1_IntKi )

         ! Now adjust the upper limit of the range slightly for the integer number of steps
      DvrSettings%ZRange(2)=  DvrSettings%ZRange(1)   + ( DvrSettings%GridN(3) - 1_IntKi ) * DvrSettings%GridDelta(3)

         ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
      IF ( .NOT. EqualRealNos(DvrSettings%ZRange(2),maxval(CLSettings%ZRange)) )  THEN
         CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of ZRange to '//TRIM(Num2LStr(DvrSettings%ZRange(2)))//   &
            ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(3)))//' points with a requested spacing '//       &
            TRIM(Num2LStr(DvrSettings%GridDelta(3)))//'.',                                                           &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF

   ELSEIF   ( ( .NOT. CLFlags%ZRange ) .AND. CLFlags%Dz )   THEN

         ! Make sure we have points to calculate
      IF ( DvrSettings%GridN(3)  == 0 )   THEN
            ! We only issue this warning.
         CALL SetErrStat( ErrID_Warn,' No points specified in driver input file for Z direction.  Ignoring Dz input.',  &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
         DvrFlags%Dz       =  .FALSE.
         DvrFlags%ZRange   =  .FALSE.
      ELSE

         DvrFlags%Dz       =  .TRUE.
         DvrFlags%ZRange   =  .TRUE.
         WindGridModify    =  .TRUE.

            ! Save the Dz value
         DvrSettings%GridDelta(3)   =  abs(CLSettings%GridDelta(3))  ! just in case there was a negative number

            ! Calculate the new value of Ny
         DvrSettings%GridN(3) =  CEILING( abs( ( DvrSettings%ZRange(2) - DvrSettings%ZRange(1) ) / DvrSettings%GridDelta(3) ) + 1_IntKi )

            ! Now adjust the upper limit of the range slightly for the integer number of steps
         DvrSettings%ZRange(2)=  DvrSettings%ZRange(1)   + ( DvrSettings%GridN(3) - 1_IntKi ) * DvrSettings%GridDelta(3)

            ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
         IF ( .NOT. EqualRealNos(DvrSettings%ZRange(2),maxval(CLSettings%ZRange)) )  THEN
            CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of ZRange to '//TRIM(Num2LStr(DvrSettings%ZRange(2)))//   &
               ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(3)))//' points with a requested spacing '//       &
               TRIM(Num2LStr(DvrSettings%GridDelta(3)))//'.',                                                           &
               ErrStat,ErrMsg,'UpdateSettingsWithCL')
         ENDIF

      ENDIF


   ELSEIF   ( CLFlags%ZRange .AND. CLFlags%Dz ) THEN

         ! Set flags if not already set
      DvrFlags%Dz       =  .TRUE.
      DvrFlags%ZRange   =  .TRUE.
      WindGridModify    =  .TRUE.

         ! Copy over the range and stepsize
      DvrSettings%ZRange(1)      =  minval(CLSettings%ZRange)     ! just in case the order is funky
      DvrSettings%ZRange(2)      =  maxval(CLSettings%ZRange)
      DvrSettings%GridDelta(3)   =  abs(CLSettings%GridDelta(3))  ! just in case there was a negative number

         ! Set number of points in Y direction
      DvrSettings%GridN(3) =  CEILING( abs( ( DvrSettings%ZRange(2) - DvrSettings%ZRange(1) ) / DvrSettings%GridDelta(3) ) + 1_IntKi )

         ! Now adjust the upper limit of the range slightly for the integer number of steps
      DvrSettings%ZRange(2)=  DvrSettings%ZRange(1)   + ( DvrSettings%GridN(3) - 1_IntKi ) * DvrSettings%GridDelta(3)

         ! Now issue a warning if we had to change the upper limit, otherwise silently adopt it.
      IF ( .NOT. EqualRealNos(DvrSettings%ZRange(2),maxval(CLSettings%ZRange)) )  THEN
         CALL SetErrStat(ErrID_Warn,' Adjusted upper limit of ZRange to '//TRIM(Num2LStr(DvrSettings%ZRange(2)))//   &
            ' so that there are '//TRIM(Num2LStr(DvrSettings%GridN(3)))//' points with a requested spacing '//       &
            TRIM(Num2LStr(DvrSettings%GridDelta(3)))//'.',                                                           &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF

   ELSE
      ! Nothing was specified, so we don't do anything
   ENDIF



      !--------------------------------------------
      ! Did we change the FFT info?
      !--------------------------------------------

   IF ( CLFlags%FFTcalc ) THEN
         !
      IF ( CLSettings%FFTcoord(3) <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Warn,' FFT coordinate ['//TRIM(Num2LStr(CLSettings%FFTcoord(1)))//','//  &
            TRIM(Num2LStr(CLSettings%FFTcoord(1)))//','//TRIM(Num2LStr(CLSettings%FFTcoord(1)))//        &
            '] is at or below ground level where there is no wind.  Ingoring.',   &
            ErrStat,ErrMsg,'UpdateSettingsWithCL' )
         DvrFlags%FFTcalc     =  .FALSE.
      ELSE
            ! If we are overriding driver input file settings, tell user
         IF ( DvrFlags%FFTcalc ) THEN
            CALL SetErrStat( ErrID_Warn,' Overriding driver input file settings for FFT calculation.',   &
               ErrStat,ErrMsg,'UpdateSettingsWithCL' )
         ENDIF
         DvrSettings%FFTcoord =  CLSettings%FFTcoord
         DvrFlags%FFTcalc     =  .TRUE.
      ENDIF

   ENDIF

      !--------------------------------------------
      ! Did we request a summary file?
      !--------------------------------------------

   IF ( CLFlags%Summary ) THEN
      DvrFlags%Summary  =  .TRUE.
   ENDIF

      !--------------------------------------------
      ! Did we request a different Points file?
      !--------------------------------------------

   IF ( CLFlags%PointsFile ) THEN
         ! If a name was given in the driver input file, then warn the user.
      IF ( DvrFlags%PointsFile ) THEN
         CALL SetErrStat( ErrID_Warn,' Overriding driver input file settings for Points file.',  &
            ErrStat,ErrMsg,'UpdateSettingsWithCL' )
      ENDIF
      DvrFlags%PointsFile        =  .TRUE.
      DvrSettings%PointsFileName =  CLSettings%PointsFileName
   ENDIF



      !--------------------------------------------
      ! If there was no driver input file, we need to set a few things.
      !--------------------------------------------

   IF ( .NOT. DVRIPT ) THEN

         ! Do we need to set the NumTimeStepsDefault flag?
      IF ( .NOT. DvrFlags%NumTimeSteps )  THEN
         DvrFlags%NumTimeStepsDefault  =  .TRUE.
         CALL SetErrStat( ErrID_Info,' The number of timesteps is not specified.  Defaulting to what is in the wind file.', &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF

         ! Do we need to set the DTDefault flag?
      IF ( .NOT. DvrFlags%DT )  THEN
         DvrFlags%DTDefault  =  .TRUE.
         CALL SetErrStat( ErrID_Info,' The timestep size is not specified.  Defaulting to what is in the wind file.', &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF


!FIXME:  Anything else?
      IF ( .NOT. DvrFlags%ZRange ) THEN
         DvrSettings%ZRange   =  50.0
         CALL SetErrStat( ErrID_Info,' ZRange not set.  Using value of 50.0.',   &
            ErrStat,ErrMsg,'UpdateSettingsWithCL')
      ENDIF
   ENDIF


      ! If no DT value has been set (DEFAULT requested), we need to set a default to pass into IfW
   IF ( .NOT. DvrFlags%DT ) THEN
      DvrSettings%DT =  0.025_DbKi     ! This value gets passed into the IfW_Init routine, so something must be set.
   ENDIF



END SUBROUTINE UpdateSettingsWithCL


SUBROUTINE ReadPointsFile( PointsFileName, CoordList, ErrStat, ErrMsg )

   CHARACTER(1024),                    INTENT(IN   )  :: PointsFileName       !< Name of the points file to read
   REAL(ReKi), ALLOCATABLE,            INTENT(  OUT)  :: CoordList(:,:)       !< The coordinates we read in
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< The error status
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< The message for the status

      ! Local variables
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error message for calls
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   INTEGER(IntKi)                                     :: FiUnitPoints         !< Unit number for points file to open

   INTEGER(IntKi)                                     :: NumDataColumns       !< Number of data columns
   INTEGER(IntKi)                                     :: NumDataPoints        !< Number of lines of data (one point per line)
   INTEGER(IntKi)                                     :: NumHeaderLines       !< Number of header lines to ignore

   INTEGER(IntKi)                                     :: I                    !< Generic counter

      ! Initialization of subroutine
   ErrMsg      =  ''
   ErrMsgTmp   =  ''
   ErrStat     =  ErrID_None
   ErrStatTmp  =  ErrID_None


      ! Now open file
   CALL GetNewUnit(    FiUnitPoints )
   CALL OpenFInpFile(   FiUnitPoints,  TRIM(PointsFileName), ErrStatTmp, ErrMsgTmp )   ! Unformatted input file
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'ReadPointsFile')
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF

      ! Find out how long the file is
   CALL GetFileLength( FiUnitPoints, PointsFileName, NumDataColumns, NumDataPoints, NumHeaderLines, ErrMsgTmp, ErrStatTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'ReadPointsFile')
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF
   IF ( NumDataColumns /= 3 ) THEN
      CALL SetErrStat( ErrID_Fatal,' Expecting three columns in '//TRIM(PointsFileName)//' corresponding to '//   &
         'X, Y, and Z coordinates.  Instead found '//TRIM(Num2LStr(NumDataColumns))//'.', &
         ErrStat, ErrMsg, 'ReadPointsFile')
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Allocate the storage for the data
   CALL AllocAry( CoordList, 3, NumDataPoints, "Array of Points data", ErrStatTmp, ErrMsgTmp )
   IF ( ErrStatTmp >= AbortErrLev ) THEN
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'ReadPointsFile')
      CLOSE( FiUnitPoints )
      RETURN
   ENDIF


      ! Read in the headers and throw them away
   DO I=1,NumHeaderLines
      CALL ReadCom( FiUnitPoints, PointsFileName,' Points file header line', ErrStatTmp, ErrMsgTmp )
      IF ( ErrStatTmp /= ErrID_None ) THEN
         CALL SetErrStat(ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadPointsFile')
         CLOSE( FiUnitPoints )
         RETURN
      ENDIF
   ENDDO

      ! Read in the datapoints
   DO I=1,NumDataPoints
      CALL ReadAry ( FiUnitPoints, PointsFileName, CoordList(:,I), 3, 'CoordList', &
         'Coordinate point from Points file', ErrStatTmp, ErrMsgTmp)
      IF ( ErrStat /= ErrID_None ) THEN
         CALL SetErrStat( ErrID_Fatal,ErrMsgTmp,ErrStat,ErrMsg,'ReadPointsFile')
         CLOSE( FiUnitPoints )
         RETURN
      ENDIF
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
                           TRIM(DataFileName)//'.', ErrStat, ErrMsg, 'GetFileLength')
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
                           ErrStat, ErrMsg, 'GetFileLength')
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



SUBROUTINE WindGridMessage( Settings, ToFile, Msg, MsgLen )

   TYPE(IfWDriver_Settings),           INTENT(IN   )  :: Settings
   LOGICAL,                            INTENT(IN   )  :: ToFile         !< Prepend comment character
   CHARACTER(2048),                    INTENT(  OUT)  :: Msg
   INTEGER(IntKi),                     INTENT(  OUT)  :: MsgLen

      ! Local Variables
   CHARACTER(11)                                      :: TmpNumString
   INTEGER(IntKi)                                     :: ErrStatTmp

   Msg   =  ''

   IF ( ToFile ) THEN
      Msg='#  '
   ELSE
      Msg="Requested wind grid data will be written to "//TRIM(Settings%WindGridOutput%Name)//'.'
   ENDIF
   Msg   =  TRIM(Msg)//"  Requested data:"//NewLine

   ! Header info:
   Msg   =  TRIM(Msg)
   IF ( ToFile ) Msg=TRIM(Msg)//'#'
   Msg   =           TRIM(Msg)//"    Dimension             Range             Stepsize       Num. points"//NewLine
   IF ( ToFile ) Msg=TRIM(Msg)//'#'
   Msg   =           TRIM(Msg)//"   -------------------------------------------------------------------"//NewLine
   ! X direction
   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%XRange(1)
   IF ( ToFile ) THEN
      Msg=  TRIM(Msg)//"#       X"
   ELSE
      Msg=  TRIM(Msg)//"        X"
   ENDIF
   MsgLen=  LEN_TRIM(Msg)
   Msg   =  Msg(1:MsgLen)//"       "//TmpNumString(1:10)
   MsgLen=  MsgLen+7+10

   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%XRange(2)
   Msg   =  Msg(1:MsgLen)//" -> "//TmpNumString(1:10)
   MsgLen=  MsgLen+4+10

   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%GridDelta(1)
   Msg   =  Msg(1:MsgLen)//"    "//TmpNumString(1:10)
   MsgLen=  MsgLen+4+10

   WRITE(TmpNumString,'(i6)',IOSTAT=ErrStatTmp)             Settings%GridN(1)
   Msg   =  Msg(1:MsgLen)//"      "//TmpNumString(1:6)//NewLine
   MsgLen=  MsgLen+6+6


   ! Y direction
   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%YRange(1)
   IF ( ToFile ) THEN
      Msg=  TRIM(Msg)//"#       Y"
   ELSE
      Msg=  TRIM(Msg)//"        Y"
   ENDIF
   MsgLen=  LEN_TRIM(Msg)
   Msg   =  Msg(1:MsgLen)//"       "//TmpNumString(1:10)
   MsgLen=  MsgLen+7+10

   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%YRange(2)
   Msg   =  Msg(1:MsgLen)//" -> "//TmpNumString(1:10)
   MsgLen=  MsgLen+4+10

   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%GridDelta(2)
   Msg   =  Msg(1:MsgLen)//"    "//TmpNumString(1:10)
   MsgLen=  MsgLen+4+10

   WRITE(TmpNumString,'(i6)',IOSTAT=ErrStatTmp)             Settings%GridN(2)
   Msg   =  Msg(1:MsgLen)//"      "//TmpNumString(1:6)//NewLine
   MsgLen=  MsgLen+6+6


   ! Z direction
   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%ZRange(1)
   IF ( ToFile ) THEN
      Msg=  TRIM(Msg)//"#       Z"
   ELSE
      Msg=  TRIM(Msg)//"        Z"
   ENDIF
   MsgLen=  LEN_TRIM(Msg)
   Msg   =  Msg(1:MsgLen)//"       "//TmpNumString(1:10)
   MsgLen=  MsgLen+7+10

   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%ZRange(2)
   Msg   =  Msg(1:MsgLen)//" -> "//TmpNumString(1:10)
   MsgLen=  MsgLen+4+10

   WRITE(TmpNumString,'(f7.2)',IOSTAT=ErrStatTmp)           Settings%GridDelta(3)
   Msg   =  Msg(1:MsgLen)//"    "//TmpNumString(1:10)
   MsgLen=  MsgLen+4+10

   WRITE(TmpNumString,'(i6)',IOSTAT=ErrStatTmp)             Settings%GridN(3)
   Msg   =  Msg(1:MsgLen)//"      "//TmpNumString(1:6)//NewLine
   MsgLen=  MsgLen+6+6


   ! T direction
   WRITE(TmpNumString,'(f10.4)',IOSTAT=ErrStatTmp)           Settings%TStart
   IF ( ToFile ) THEN
      Msg=  TRIM(Msg)//"#       T"
   ELSE
      Msg=  TRIM(Msg)//"        T"
   ENDIF
   MsgLen=  LEN_TRIM(Msg)
   Msg   =  Msg(1:MsgLen)//"      "//TmpNumString(1:11)
   MsgLen=  MsgLen+7+11

   WRITE(TmpNumString,'(f10.4)',IOSTAT=ErrStatTmp)           Settings%TStart+Settings%DT*Settings%NumTimeSteps
   Msg   =  Msg(1:MsgLen)//"->"//TmpNumString(1:11)
   MsgLen=  MsgLen+4+11

   WRITE(TmpNumString,'(f10.4)',IOSTAT=ErrStatTmp)           Settings%DT
   Msg   =  Msg(1:MsgLen)//" "//TmpNumString(1:11)
   MsgLen=  MsgLen+4+11

   WRITE(TmpNumString,'(i8)',IOSTAT=ErrStatTmp)             Settings%NumTimeSteps
   Msg   =  Msg(1:MsgLen)//" "//TmpNumString(1:8)  !//NewLine
   MsgLen=  MsgLen+6+8

END SUBROUTINE


!> This subroutine outputs the results of the WindGrid calculations information at each timestep.
SUBROUTINE WindGridVel_OutputWrite (OutFile, Settings, GridXYZ, GridVel, TIME, ErrStat, ErrMsg)
   TYPE(OutputFile),                   INTENT(INOUT)  :: OutFile
   TYPE(IfWDriver_Settings),           INTENT(IN   )  :: Settings             !< Settings for IfW driver
   REAL(ReKi),          ALLOCATABLE,   INTENT(IN   )  :: GridXYZ(:,:)         !< The position grid passed in
   REAL(ReKi),          ALLOCATABLE,   INTENT(IN   )  :: GridVel(:,:)         !< The velocity grid passed in
   REAL(DbKi),                         INTENT(IN   )  :: TIME                 !< The current time
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None

         ! Temporary local variables
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary variable for the status of error message
   CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary variable for the error message
   INTEGER(IntKi)                                     :: LenErrMsgTmp         !< Length of ErrMsgTmp (for getting WindGrid info)

   CHARACTER(52)                                      :: WindVelFmt           !< Format specifier for the output file for wave elevation series
   INTEGER(IntKi)                                     :: I                    !< generic counter


   WindVelFmt = "(3(F14.7,3x),3(F14.7,3x))"

   ErrMsg      = ''
   ErrStat     = ErrID_None
   ErrMsgTmp   = ''
   ErrStatTmp  = ErrID_None


      ! If it hasn't been initially written to, do this then exit. Otherwise set a few things and continue.
   IF ( .NOT. OutFile%Initialized ) THEN

      CALL GetNewUnit( OutFile%Unit )
      CALL OpenFOutFile( OutFile%Unit, TRIM(OutFile%Name), ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'WindGridVel_OutputWrite' )
      IF ( ErrStat >= AbortErrLev ) RETURN

      OutFile%Initialized =  .TRUE.

         ! Write header section
      WRITE( OutFile%Unit,'(A)', IOSTAT=ErrStatTmp )   '## This file was generated by '//TRIM(GetNVD(Settings%ProgInfo))//  &
            ' on '//CurDate()//' at '//CurTime()//'.'
      WRITE( OutFile%Unit,'(A)', IOSTAT=ErrStatTmp )   '## This file contains the wind velocity at a grid of points at each '// &
                                                   'requested timestep'
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '## It is arranged as blocks of X,Y,Z,U,V,W at each timestep'
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '## Each block is separated by two blank lines for use in gnuplot'
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      CALL WindGridMessage( Settings, .TRUE., ErrMsgTmp, LenErrMsgTmp )
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  ErrMsgTmp(1:LenErrMsgTmp)
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '#           X              Y              Z  '//  &
                                                   '            U              V              W'
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '#          (m)            (m)            (m) '//  &
                                                   '          (m/s)          (m/s)          (m/s)'
   ELSE

      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp ) NewLine//NewLine
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp ) '# Time: '//TRIM(Num2LStr(TIME))

      DO I = 1,SIZE(GridXYZ,DIM=2)

         WRITE (OutFile%Unit,WindVelFmt, IOSTAT=ErrStatTmp )    GridXYZ(1,I),GridXYZ(2,I),GridXYZ(3,I),GridVel(1,I),GridVel(2,I),GridVel(3,I)

      ENDDO

   ENDIF

END SUBROUTINE WindGridVel_OutputWrite


SUBROUTINE PointData_OutputWrite (OutFile, Settings, GridXYZ, GridVel, GridAcc, TIME, IsVel, ErrStat, ErrMsg)

   TYPE(OutputFile),                   INTENT(INOUT)  :: OutFile
   TYPE(IfWDriver_Settings),           INTENT(IN   )  :: Settings             !< Settings for IfW driver
   REAL(ReKi),                         INTENT(IN   )  :: GridXYZ(:,:)         !< The position grid passed in
   REAL(ReKi),                         INTENT(IN   )  :: GridVel(:,:)         !< The velocity grid passed in
   REAL(ReKi), allocatable,            INTENT(IN   )  :: GridAcc(:,:)         !< The acceleration grid passed in
   REAL(DbKi),                         INTENT(IN   )  :: TIME                 !< The current time
   LOGICAL,                            INTENT(IN   )  :: IsVel                !< Is this velocity output
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary variable for the status of error message
   CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary variable for the error message
   INTEGER(IntKi)                                     :: I, J                 !< Generic counter

   CHARACTER(61)                                      :: NameUnitFmt          !< Format specifier for the output file for channel names and units
   CHARACTER(61)                                      :: PointsVelFmt         !< Format specifier for the output file for wind point location and velocity

   NameUnitFmt = "( *(A16,3X) )"
   PointsVelFmt = "( *(F16.8,3X) )"

   ErrMsg      = ''
   ErrStat     = ErrID_None
   ErrMsgTmp   = ''
   ErrStatTmp  = ErrID_None


      ! If it hasn't been initially written to, do this then exit. Otherwise set a few things and continue.
   IF ( .NOT. OutFile%Initialized ) THEN

      CALL GetNewUnit( OutFile%Unit )
      CALL OpenFOutFile( OutFile%Unit, TRIM(OutFile%Name), ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, 'PointData_OutputWrite' )
      IF ( ErrStat >= AbortErrLev ) RETURN

      OutFile%Initialized =  .TRUE.

         ! Write header section
      WRITE( OutFile%Unit,'(A)', IOSTAT=ErrStatTmp )   '## This file was generated by '//TRIM(GetNVD(Settings%ProgInfo))//  &
            ' on '//CurDate()//' at '//CurTime()//'.'
      if (allocated(GridAcc)) then
         WRITE( OutFile%Unit,'(A)', IOSTAT=ErrStatTmp )   '## This file contains the wind velocity and acceleration at the '//   &
                                                      TRIM(Num2LStr(SIZE(GridXYZ,DIM=2)))//' points specified in the '// &
                                                      'file '//TRIM(Settings%PointsFileName)//'.'
      else
         WRITE( OutFile%Unit,'(A)', IOSTAT=ErrStatTmp )   '## This file contains the wind velocity at the '//   &
                                                      TRIM(Num2LStr(SIZE(GridXYZ,DIM=2)))//' points specified in the '// &
                                                      'file '//TRIM(Settings%PointsFileName)//'.'
      end if
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      WRITE (OutFile%Unit,'(A)', IOSTAT=ErrStatTmp  )  '# '
      if (allocated(GridAcc)) then
         WRITE (OutFile%Unit, NameUnitFmt, IOSTAT=ErrStatTmp  )  'T', 'X', 'Y', 'Z', 'U', 'V', 'W', 'UA', 'VA', 'WA'
         WRITE (OutFile%Unit, NameUnitFmt, IOSTAT=ErrStatTmp  )  '(s)', '(m)', '(m)', '(m)', '(m/s)', '(m/s)', '(m/s)', '(m/s/s)', '(m/s/s)', '(m/s/s)'
      else
         WRITE (OutFile%Unit, NameUnitFmt, IOSTAT=ErrStatTmp  )  'T', 'X', 'Y', 'Z', 'U', 'V', 'W'
         WRITE (OutFile%Unit, NameUnitFmt, IOSTAT=ErrStatTmp  )  '(s)', '(m)', '(m)', '(m)', '(m/s)', '(m/s)', '(m/s)'
      end if
   ELSE

      if (allocated(GridAcc)) then
         DO I = 1,SIZE(GridXYZ,DIM=2)
            WRITE (OutFile%Unit, PointsVelFmt, IOSTAT=ErrStatTmp) TIME, (GridXYZ(J,I),j=1,3), (GridVel(J,I),j=1,3), (GridAcc(J,I),j=1,3)
         ENDDO
      else
         DO I = 1,SIZE(GridXYZ,DIM=2)
            WRITE (OutFile%Unit, PointsVelFmt, IOSTAT=ErrStatTmp) TIME, (GridXYZ(J,I),j=1,3), (GridVel(J,I),j=1,3)
         ENDDO
      end if

   ENDIF

END SUBROUTINE PointData_OutputWrite


subroutine IfW_WriteUniform(FF, FileRootName, ErrStat, ErrMsg)

   type(FlowFieldType), intent(in)  :: FF             !< Parameters
   character(*), intent(in)         :: FileRootName   !< RootName for output files
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = "IfW_WriteUniform"
   type(UniformFieldType)           :: UF
   integer(IntKi)                   :: unit
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get new unit for writing file
   call GetNewUnit(unit, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Switch based on field type
   select case (FF%FieldType)

   case (Uniform_FieldType)

      call Uniform_WriteHH(FF%Uniform, FileRootName, unit, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   case (Grid3D_FieldType)

      call Grid3D_to_Uniform(FF%Grid3D, UF, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat < AbortErrLev) then
         call Uniform_WriteHH(UF, FileRootName, unit, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if

   case default

      ErrStat = ErrID_Warn
      ErrMsg = RoutineName//': Field type '//TRIM(Num2LStr(FF%FieldType))// &
               ' cannot be converted to UniformWind format.'
   end select

end subroutine IfW_WriteUniform

subroutine IfW_WriteHAWC(FF, FileRootName, ErrStat, ErrMsg)

   type(FlowFieldType), intent(in)  :: FF             !< Parameters
   character(*), intent(in)         :: FileRootName   !< RootName for output files
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = "IfW_Convert2HAWC"
   type(Grid3DFieldType)            :: G3D
   integer(IntKi)                   :: unit
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2

   ! Get new unit for writing file
   call GetNewUnit(unit, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Switch based on field type
   select case (FF%FieldType)

   case (Uniform_FieldType)

      call Uniform_to_Grid3D(FF%Uniform, FF%VelInterpCubic, G3D, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat < AbortErrLev) then
         call Grid3D_WriteHAWC(G3D, FileRootName, unit, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if

   case (Grid3D_FieldType)

      call Grid3D_WriteHAWC(FF%Grid3D, FileRootName, unit, ErrStat, ErrMsg)

   case default

      ErrStat = ErrID_Warn
      ErrMsg = RoutineName//': Field Type '//TRIM(Num2LStr(FF%FieldType))// &
               ' cannot be converted to HAWC format.'

   end select

end subroutine

subroutine IfW_WriteBladed(FF, FileRootName, ErrStat, ErrMsg)

   type(FlowFieldType), intent(in)  :: FF             !< Parameters
   character(*), intent(in)         :: FileRootName   !< RootName for output files
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = "IfW_WriteBladed"
   type(Grid3DFieldType)            :: G3D
   integer(IntKi)                   :: unit
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get new unit for writing file
   call GetNewUnit(unit, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Switch based on field type
   select case (FF%FieldType)

   case (Uniform_FieldType)

      call Uniform_to_Grid3D(FF%Uniform, FF%VelInterpCubic, G3D, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat < AbortErrLev) then
         call Grid3D_WriteBladed(G3D, FileRootName, unit, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if

   case (Grid3D_FieldType)

      call Grid3D_WriteBladed(FF%Grid3D, FileRootName, unit, ErrStat, ErrMsg)

   case default

      ErrStat = ErrID_Warn
      ErrMsg = RoutineName//': Field type '//TRIM(Num2LStr(FF%FieldType))// &
               ' cannot be converted to Bladed format.'

   end select

end subroutine


subroutine IfW_WriteVTK(FF, FileRootName, ErrStat, ErrMsg)

   type(FlowFieldType), intent(in)  :: FF             !< Parameters
   character(*), intent(in)         :: FileRootName   !< RootName for output files
   integer(IntKi), intent(out)      :: ErrStat        !< Error status of the operation
   character(*), intent(out)        :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter          :: RoutineName = "IfW_WriteVTK"
   type(Grid3DFieldType)            :: G3D
   integer(IntKi)                   :: unit
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Get new unit for writing file
   call GetNewUnit(unit, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Switch based on field type
   select case (FF%FieldType)

   case (Uniform_FieldType)

      call Uniform_to_Grid3D(FF%Uniform, FF%VelInterpCubic, G3D, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat < AbortErrLev) then
         call Grid3D_WriteVTK(G3D, FileRootName, unit, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if

   case (Grid3D_FieldType)

      call Grid3D_WriteVTK(FF%Grid3D, FileRootName, unit, ErrStat, ErrMsg)

   case default

      ErrStat = ErrID_Warn
      ErrMsg = RoutineName//': Field type '//TRIM(Num2LStr(FF%FieldType))// &
               ' cannot be converted to VTK format.'

   end select

end subroutine IfW_WriteVTK


!> This routine exists only to support the development of the module.  It will not be needed after the module is complete.
SUBROUTINE  printSettings( DvrFlags, DvrSettings )
      ! The arguments
   TYPE( IfWDriver_Flags ),            INTENT(IN   )  :: DvrFlags           !< Flags indicating which settings were set
   TYPE( IfWDriver_Settings ),         INTENT(IN   )  :: DvrSettings        !< Stored settings
   integer(IntKi)                                     :: i

   CALL WrsCr(TRIM(GetNVD(DvrSettings%ProgInfo)))
   CALL WrScr(' DvrIptFile:          '//FLAG(DvrFlags%DvrIptFile)//        '      '//TRIM(DvrSettings%DvrIptFileName))
   CALL WrScr(' IfWIptFile:          '//FLAG(DvrFlags%IfWIptFile)//        '      '//TRIM(DvrSettings%IfWIptFileName))
   CALL WrScr(' PointsFile:          '//FLAG(DvrFlags%PointsFile)//        '      '//TRIM(DvrSettings%PointsFileName))
   CALL WrScr(' FFTOutputName:       '//FLAG(DvrFlags%FFTcalc)//           '      '//TRIM(DvrSettings%FFTOutput%Name))
   CALL WrScr(' WindGridOutName:     '//FLAG(DvrFlags%WindGrid)//          '      '//TRIM(DvrSettings%WindGridOutput%Name))
   CALL WrScr(' Summary:             '//FLAG(DvrFlags%Summary))
   CALL WrScr(' SummaryFile:         '//FLAG(DvrFlags%SummaryFile)//       '      '//TRIM(DvrSettings%SummaryFileName))
   CALL WrScr(' TStart:              '//FLAG(DvrFlags%TStart)//            '      '//TRIM(Num2LStr(DvrSettings%TStart)))
   IF ( DvrFlags%DTDefault) THEN
      CALL WrScr(' DT:                  '//FLAG(DvrFlags%DT)//                '      DEFAULT')
   ELSE
      CALL WrScr(' DT:                  '//FLAG(DvrFlags%DT)//                '      '//TRIM(Num2LStr(DvrSettings%DT)))
   ENDIF
   IF ( DvrFlags%NumTimeStepsDefault) THEN
      CALL WrScr(' NumTimeSteps:        '//FLAG(DvrFlags%NumTimeSteps)//      '      DEFAULT')
   ELSE
      CALL WrScr(' NumTimeSteps:        '//FLAG(DvrFlags%NumTimeSteps)//      '      '//TRIM(Num2LStr(DvrSettings%NumTimeSteps)))
   ENDIF
   CALL WrScr(' FFTcalc:             '//FLAG(DvrFlags%FFTcalc)//           '     ['//TRIM(Num2LStr(DvrSettings%FFTcoord(1)))//', '&
                                                                                   //TRIM(Num2LStr(DvrSettings%FFTcoord(2)))//', '&
                                                                                   //TRIM(Num2LStr(DvrSettings%FFTcoord(3)))//']')
   CALL WrScr(' WindGrid:            '//FLAG(DvrFlags%WindGrid))
if (DvrFlags%WindGrid) then
   CALL WrScr(' GridN:                      '                                      //TRIM(Num2LStr(DvrSettings%GridN(1)))//' x '   &
                                                                                   //TRIM(Num2LStr(DvrSettings%GridN(2)))//' x '   &
                                                                                   //TRIM(Num2LStr(DvrSettings%GridN(3))))
   CALL WrScr(' XRange:              '//FLAG(DvrFlags%XRange)//            '      '//TRIM(Num2LStr(DvrSettings%XRange(1)))//' -- ' &
                                                                                   //TRIM(Num2LStr(DvrSettings%XRange(2))))
   CALL WrScr(' YRange:              '//FLAG(DvrFlags%YRange)//            '      '//TRIM(Num2LStr(DvrSettings%YRange(1)))//' -- ' &
                                                                                   //TRIM(Num2LStr(DvrSettings%YRange(2))))
   CALL WrScr(' ZRange:              '//FLAG(DvrFlags%ZRange)//            '      '//TRIM(Num2LStr(DvrSettings%ZRange(1)))//' -- ' &
                                                                                   //TRIM(Num2LStr(DvrSettings%ZRange(2))))
   CALL WrScr(' Dx:                  '//FLAG(DvrFlags%Dx)//                '      '//TRIM(Num2LStr(DvrSettings%GridDelta(1))))
   CALL WrScr(' Dy:                  '//FLAG(DvrFlags%Dy)//                '      '//TRIM(Num2LStr(DvrSettings%GridDelta(2))))
   CALL WrScr(' Dz:                  '//FLAG(DvrFlags%Dz)//                '      '//TRIM(Num2LStr(DvrSettings%GridDelta(3))))
   CALL WrScr(' WindGridOutputInit:  '//FLAG(DvrSettings%WindGridOutput%Initialized)//   '      Unit #:  '//TRIM(Num2LStr(DvrSettings%WindGridOutput%Unit)))
end if
   CALL WrScr(' FFTOutputInit:       '//FLAG(DvrSettings%FFTOutput%Initialized)//        '      Unit #:  '//TRIM(Num2LStr(DvrSettings%FFTOutput%Unit)))
   CALL WrScr(' PointsVelOutputInit: '//FLAG(DvrSettings%PointsVelOutput%Initialized)//  '      Unit #:  '//TRIM(Num2LStr(DvrSettings%PointsVelOutput%Unit)))
   CALL WrScr(' PointsAccOutputInit: '//FLAG(DvrSettings%PointsVelOutput%Initialized)//  '      Unit #:  '//TRIM(Num2LStr(DvrSettings%PointsVelOutput%Unit)))
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
