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
   CHARACTER(1024),                    INTENT(  OUT)  :: ErrMsg

   ErrStat  =   ErrID_None
   ErrMsg   =   ''


      !  Statement about usage
   CALL WrScr("")
   CALL WrScr("  Syntax:  InlowWind_Driver <filename> [options]")
   CALL WrScr("")
   CALL WrScr("       where:     <filename>     -- Name of wind file")
   CALL WrScr("     options:     "//SwChar//"type[<type>]  -- type of the file, where <type> is               [N/A]")
   CALL WrScr("                                 Uniform  -- Uniform                                     ")
   CALL WrScr("                                    FF    -- Full Field                                  ")
!   CALL WrScr("                                    CTP   -- Coherent turbulence wind field on top of another")
   CALL WrScr("                                    HAWC  -- HAWC formatted file                         ")
   CALL WrScr("                                    User  -- User Defined                                ")
   CALL WrScr("                  "//SwChar//"height[#]     -- height of the hub                         ")
   CALL WrScr("                  "//SwChar//"width[#]      -- width of the windfield                    ")
   CALL WrScr("                  "//SwChar//"x[#:#]        -- range of x (#'s are reals)                      [N/A]")
   CALL WrScr("                  "//SwChar//"y[#:#]        -- range of y                                      [N/A]")
   CALL WrScr("                  "//SwChar//"z[#:#]        -- range in z (ground = 0.0)                       [N/A]")
   CALL WrScr("                  "//SwChar//"t[#:#]        -- range of time                                   [N/A]")
   CALL WrScr("                  "//SwChar//"xres[#]       -- resolution in x                                 [N/A]")
   CALL WrScr("                  "//SwChar//"yres[#]       -- resolution in y                                 [N/A]")
   CALL WrScr("                  "//SwChar//"zres[#]       -- resolution in z                                 [N/A]")
   CALL WrScr("                  "//SwChar//"tres[#]       -- resolution in time                              [N/A]")
   CALL WrScr("                  "//SwChar//"paraprint     -- make an output file for ParaView                [N/A]")
   CALL WrScr("                  "//SwChar//"summary       -- summarize in  .sum file                         [N/A]")
   CALL WrScr("                  "//SwChar//"fft[X,Y,Z]    -- an fft over all t at X,Y,Z (outputs .fft file)  [N/A]")
   CALL WrScr("                  "//SwChar//"points[FILE]  -- calculates at x,y,z coordinates specified in a  [N/A]")
   CALL WrScr("                                    white space delimited FILE")
   CALL WrScr("                  "//SwChar//"help          -- print this help menu and exits")
   CALL WrScr("")
   CALL WrScr("   If the type is not specified, attempts are made to figure out what it is.")
   CALL WrScr("   Unspecified ranges and resolutions default to what is in the file.")
   CALL WrScr("   Features marked [N/A] have not been implimented in this version.")
!FIXME: Does the CTP get used with another? If so, specify in comment at end.
   CALL WrScr("")


END SUBROUTINE DispHelpText


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
SUBROUTINE RetrieveArgs( Settings, SettingsFlags, ErrStat, ErrMsg )
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
      ! Iterate through the input arguments !
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

   USE NWTC_Library
   USE InflowWind_Driver_Types

   IMPLICIT NONE

      ! Storing the arguments
   TYPE( InflowWind_Driver_ArgFlags ),        INTENT(  OUT)  :: SettingsFlags        ! Flags indicating which arguments were specified
   TYPE( InflowWind_Driver_Args ),            INTENT(  OUT)  :: Settings             ! Arguments passed in

      ! Error Handling
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(1024),                    INTENT(  OUT)  :: ErrMsg

      ! Local variables
   INTEGER(IntKi)                                     :: i
   CHARACTER(1024)                                    :: arg
   INTEGER(IntKi)                                     :: InputArgs
   LOGICAL                                            :: FileNameGiven


      ! initialize some things
   FileNameGiven = .FALSE.
   ErrStat = ErrID_None


      ! Check how many arguments are passed in
   InputArgs = COMMAND_ARGUMENT_COUNT()

      ! exit if we don't have enough
   IF (InputArgs == 0) THEN
      ErrMsg   = "Insufficient Arguments.  Use option "//SwChar//"help for help menu."
      ErrStat  = ErrID_Fatal     ! Cannot continue
      RETURN
   ENDIF


      ! Loop through all the arguments, and store them
   DO i=1,InputArgs
      CALL get_command_argument(i, arg)

         ! Check to see if it is a control parameter or the filename
      IF ( INDEX( SwChar, arg(1:1) ) > 0 ) THEN

            ! check to see if we asked for help
         IF ( arg(2:5) == "help" ) THEN
            CALL DispHelpText( ErrStat, ErrMsg )
            CALL ProgExit(0)
         ENDIF

            ! Check the argument and put it where it belongs
            ! chop the SwChar off before passing the argument
         CALL ParseArg( Settings, SettingsFlags, arg(2:), ErrStat, ErrMsg )

         IF ( ErrStat == ErrID_Warn ) CALL ProgWarn( ErrMsg )
         IF ( ErrStat == ErrID_Fatal ) RETURN


      ELSE

            ! since there is no switch character, assume it is the filename, unless we already set one
         IF ( FileNameGiven ) THEN
            ErrMsg   = "Multiple input filenames given: "//TRIM(Settings%InputFileName)//", "//TRIM(arg)
            ErrStat  = ErrID_Fatal
            RETURN
         ELSE
            Settings%InputFileName = TRIM(arg)
            FileNameGiven = .TRUE.
         ENDIF

      ENDIF
   END DO


      ! Check the arguments passed in:
   IF ( .NOT. FileNameGiven ) THEN
      ErrMsg   = "No filename given for file to open."
      ErrStat  = ErrID_Fatal
      RETURN
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
      INTEGER(IntKi)                                     :: TempIO         ! Temporary variable to hold the error status

         read( StringIn, *, iostat=TempIO) StringToReal

            ! If that isn't a number, only warn since we can continue by skipping this value
         IF ( TempIO .ne. 0 ) ErrSTat  = ErrID_Warn

   END FUNCTION StringToReal



   !-------------------------------------------------------------------------------
   SUBROUTINE ParseArg( Settings, SettingsFlags, ThisArg, ErrStat, ErrMsg )
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
         ! Parse and store the input argument  !
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

      USE NWTC_Library
      USE InflowWind_Driver_Types
      USE InflowWind_Types

      IMPLICIT NONE

         ! Storing the arguments
      TYPE( InflowWind_Driver_ArgFlags ),        INTENT(INOUT)  :: SettingsFlags        ! Flags indicating which arguments were specified
      TYPE( InflowWind_Driver_Args ),            INTENT(INOUT)  :: Settings             ! Arguments passed in

      CHARACTER(*),                       INTENT(IN   )  :: ThisArg              ! The current argument

         ! Error Handling
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(1024),                    INTENT(  OUT)  :: ErrMsg


         ! local variables
      INTEGER(IntKi)                                     :: Delim1               ! where the [ is
      INTEGER(IntKi)                                     :: Delim2               ! where the ] is
      INTEGER(IntKi)                                     :: DelimSep             ! where the : is
      INTEGER(IntKi)                                     :: DelimSep2            ! where the : is
      INTEGER(IntKi)                                     :: DelimSep3            ! where the : is
      REAL(ReKi)                                         :: TempReal             ! temp variable to hold a real
      INTEGER(IntKi)                                     :: TempIO               ! temp variable to store the IO error



         ! Initialize some things
      ErrStat  = ErrID_None
      ErrMsg   = ""
      TempIO   = 0

         ! Get the delimiters -- returns 0 if there isn't one
      Delim1   = INDEX(ThisArg,'[')
      Delim2   = INDEX(ThisArg,']')
      DelimSep = INDEX(ThisArg,':')


         ! check that if there is an opening bracket, then there is a closing one
      IF ( (Delim1 > 0 ) .and. (Delim2 < Delim1) ) THEN
         ErrMsg   = TRIM(ErrMsg)//NewLine//"Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
         ErrStat  = ErrID_Warn
         RETURN
      ENDIF

         ! check that if there is a colon, then there are brackets
      IF ( (DelimSep > 0) .and. (Delim1 == 0) ) THEN
         ErrMsg   = TRIM(ErrMsg)//NewLine//"Syntax error in option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
         ErrStat  = ErrID_Warn
         RETURN
      ENDIF


         ! Now go through the full list of possible options. Store as appropriate

         ! "type[#]"
      IF    ( ThisArg(1:Delim1) == "type["      ) THEN
         SELECT CASE (ThisArg(Delim1+1:Delim2-1))
            CASE ('Steady')    ! hub height
               SettingsFlags%WindFileType = .TRUE.
               Settings%WindFileType      = Steady_WindNumber

            CASE ('Uniform')    ! hub height
               SettingsFlags%WindFileType = .TRUE.
               Settings%WindFileType      = Uniform_WindNumber

            CASE ('FF')    ! full field
               SettingsFlags%WindFileType = .TRUE.
               Settings%WindFileType      = FF_WindNumber

!            CASE ('CTP')   ! coherent turbulence
!               SettingsFlags%WindFileType = .TRUE.
!               Settings%WindFileType      = CTP_WindNumber

            CASE ('HAWC')  ! HAWC compatible
               SettingsFlags%WindFileType = .TRUE.
               Settings%WindFileType      = HAWC_WindNumber

            CASE ('User')    ! User Defined
               SettingsFlags%WindFileType = .TRUE.
               Settings%WindFileType      = User_WindNumber

            CASE DEFAULT
               ErrMsg   = "Invalid wind type. Ignoring option: '"//SwChar//TRIM(ThisArg)//"'."
               ErrStat  = ErrID_Warn

         END SELECT



         ! "height[#]"
      ELSEIF( ThisArg(1:Delim1) == "height["    ) THEN
         TempReal = StringToReal( ThisArg(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%Height = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%Height = .TRUE.
            Settings%Height      = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%Height = .FALSE.
         ENDIF


         ! "width[#]"
      ELSEIF( ThisArg(1:Delim1) == "width["     ) THEN
         TempReal = StringToReal( ThisArg(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%Width  = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%Width  = .TRUE.
            Settings%Width       = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%Width  = .FALSE.
         ENDIF



         ! "x[#:#]"
      ELSEIF( ThisArg(1:Delim1) == "x["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArg(Delim1+1:DelimSep-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%XRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%XRange = .TRUE.
            Settings%XRange(1)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%XRange = .FALSE.
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArg(DelimSep+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%XRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%XRange = .TRUE.
            Settings%XRange(2)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%XRange = .FALSE.
         ENDIF

            ! Check the order of values
         IF ( Settings%XRange(1) > Settings%Xrange(2) ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring."
            ErrStat  = ErrID_Warn
            Settings%XRange(1)   = 0.0
            Settings%XRange(2)   = 0.0
            SettingsFlags%XRange = .FALSE.
         ENDIF


         ! "y[#:#]"
      ELSEIF( ThisArg(1:Delim1) == "y["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArg(Delim1+1:DelimSep-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%YRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%YRange = .TRUE.
            Settings%YRange(1)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%YRange = .FALSE.
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArg(DelimSep+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%YRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%YRange = .TRUE.
            Settings%YRange(2)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%YRange = .FALSE.
         ENDIF

            ! Check the order of values
         IF ( Settings%YRange(1) > Settings%Yrange(2) ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring."
            ErrStat  = 1
            Settings%YRange(1)   = 0.0
            Settings%YRange(2)   = 0.0
            SettingsFlags%YRange = .FALSE.
         ENDIF


         ! "z[#:#]"
      ELSEIF( ThisArg(1:Delim1) == "z["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArg(Delim1+1:DelimSep-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%ZRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%ZRange = .TRUE.
            Settings%ZRange(1)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%ZRange = .FALSE.
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArg(DelimSep+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%ZRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%ZRange = .TRUE.
            Settings%ZRange(2)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%ZRange = .FALSE.
         ENDIF

            ! Check the order of values
         IF ( Settings%ZRange(1) > Settings%Zrange(2) ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring."
            ErrStat  = 1
            Settings%ZRange(1)   = 0.0
            Settings%ZRange(2)   = 0.0
            SettingsFlags%ZRange = .FALSE.
         ENDIF


         ! "t[#:#]"
      ELSEIF( ThisArg(1:Delim1) == "t["         ) THEN

            ! First Value
         TempReal = StringToReal( ThisArg(Delim1+1:DelimSep-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%TRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%TRange = .TRUE.
            Settings%TRange(1)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%TRange = .FALSE.
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArg(DelimSep+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%TRange = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%TRange = .TRUE.
            Settings%TRange(2)   = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%TRange = .FALSE.
         ENDIF

            ! Check the order of values
         IF ( Settings%TRange(1) > Settings%Trange(2) ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Unexpected order of values in option '"//SwChar//TRIM(ThisArg)//"'. Ingoring."
            Settings%TRange(1)   = 0.0
            Settings%TRange(2)   = 0.0
            ErrStat  = 1
            SettingsFlags%TRange = .FALSE.
         ENDIF


         ! "xres[#]"
      ELSEIF( ThisArg(1:Delim1) == "xres["      ) THEN
         TempReal = StringToReal( ThisArg(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%XRes   = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%XRes   = .TRUE.
            Settings%XRes        = abs(TempReal)
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%XRes   = .FALSE.
         ENDIF

         ! "yres[#]"
      ELSEIF( ThisArg(1:Delim1) == "yres["      ) THEN
         TempReal = StringToReal( ThisArg(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%YRes   = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%YRes   = .TRUE.
            Settings%YRes        = abs(TempReal)
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%YRes   = .FALSE.
         ENDIF

         ! "zres[#]"
      ELSEIF( ThisArg(1:Delim1) == "zres["      ) THEN
         TempReal = StringToReal( ThisArg(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%ZRes   = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%ZRes   = .TRUE.
            Settings%ZRes        = abs(TempReal)
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%ZRes   = .FALSE.
         ENDIF

         ! "tres[#]"
      ELSEIF( ThisArg(1:Delim1) == "tres["      ) THEN
         TempReal = StringToReal( ThisArg(Delim1+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%TRes   = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%TRes   = .TRUE.
            Settings%TRes        = abs(TempReal)
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%TRes   = .FALSE.
         ENDIF

         ! "paraprint"
      ELSEIF( ThisArg(1:9)      == "paraprint"  ) THEN
         SettingsFlags%ParaPrint = .TRUE.
         ErrMsg   = TRIM(ErrMsg)//NewLine//"Feature not implimented. Ignoring option '"//SwChar//TRIM(ThisArg)//"'."
         ErrStat  = ErrID_Warn

         ! "summary"
      ELSEIF( ThisArg(1:8)      == "summary"    ) THEN
         SettingsFlags%Summary   = .TRUE.
         ErrMsg   = TRIM(ErrMsg)//NewLine//"Feature not implimented. Ignoring option '"//SwChar//TRIM(ThisArg)//"'."
         ErrStat  = ErrID_Warn

         ! "fft[X,Y,Z]"
      ELSEIF( ThisArg(1:Delim1) == "fft["       ) THEN
         DelimSep = INDEX(ThisArg,',')
         DelimSep2= INDEX(ThisArg(DelimSep+1:),',') + DelimSep

            ! First Value
         TempReal = StringToReal( ThisArg(Delim1+1:DelimSep-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%fft    = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%fft    = .TRUE.
            Settings%fft(1)      = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%fft    = .FALSE.
         ENDIF

            ! Second Value
         TempReal = StringToReal( ThisArg(DelimSep+1:DelimSep2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%fft    = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%fft    = .TRUE.
            Settings%fft(2)      = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%fft    = .FALSE.
         ENDIF

            ! Third Value
         TempReal = StringToReal( ThisArg(DelimSep2+1:Delim2-1), ErrStat )
         IF ( ErrStat == ErrID_Warn ) THEN
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Invalid number in option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            SettingsFlags%fft    = .FALSE.
            RETURN
         ELSEIF ( ErrStat == ErrID_None ) THEN
            SettingsFlags%fft    = .TRUE.
            Settings%fft(3)      = TempReal
         ELSE
            ErrMsg   = TRIM(ErrMsg)//NewLine//"Something failed in parsing option '"//SwChar//TRIM(ThisArg)//"'. Ignoring."
            ErrStat  = ErrID_Fatal
            SettingsFlags%fft    = .FALSE.
         ENDIF

         ! "points[FILE]"
      ELSEIF( ThisArg(1:Delim1) == "points["    ) THEN
         SettingsFlags%PointsFile= .TRUE.
         Settings%PointsFileName = ThisArg(Delim1+1:Delim2-1)

      ELSE
         ErrMsg  = TRIM(ErrMsg)//NewLine//"Unrecognized option: '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options."
         ErrStat = ErrID_Warn
      ENDIF

   END SUBROUTINE ParseArg
   !-------------------------------------------------------------------------------



END SUBROUTINE RetrieveArgs

!!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
END MODULE InflowWind_Driver_Subs
