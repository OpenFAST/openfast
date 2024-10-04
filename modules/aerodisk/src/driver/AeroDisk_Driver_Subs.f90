!**********************************************************************************************************************************
!
!  MODULE: AeroDisk_Driver_Subs  - This module contains subroutines used by the AeroDisk Driver program
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
!
!    This file is part of AeroDisk.
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
MODULE AeroDisk_Driver_Subs

   USE NWTC_Library
   USE AeroDisk_Driver_Types
   IMPLICIT NONE

!  NOTE: This is loosely based on the InflowWind driver code.

CONTAINS
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> Print out help information
SUBROUTINE DispHelpText()
      !  Statement about usage
   CALL WrScr("")
   CALL WrScr("  Syntax:  AeroDisk_Driver <filename> [options]")
   CALL WrScr("")
   CALL WrScr("       where:     <filename>     -- Name of driver input file to use")
   CALL WrScr("     options:     "//SWChar//"adsk          -- treat <filename> as name of AeroDisk input file")
   CALL WrScr("                                    (no driver input file)")
   CALL WrScr("")
   CALL WrScr("              The following options will overwrite values in the driver input file:")
   CALL WrScr("                  "//SwChar//"DT[#]         -- timestep                                        ")
   CALL WrScr("                  "//SwChar//"TStart[#]     -- start time                                      ")
   CALL WrScr("                  "//SwChar//"TSteps[#]     -- number of timesteps                             ")
   CALL WrScr("                  "//SwChar//"v             -- verbose output ")
   CALL WrScr("                  "//SwChar//"vv            -- very verbose output ")
   CALL WrScr("                  "//SwChar//"NonLinear     -- only return non-linear portion of reaction force")
   CALL WrScr("                  "//SwChar//"help          -- print this help menu and exit")
   CALL WrScr("")
   CALL WrScr("   Notes:")
   CALL WrScr("   -- Options are not case sensitive.")
   CALL WrScr("")
!FIXME: update this
END SUBROUTINE DispHelpText


subroutine InitSettingsFlags( ProgInfo, CLSettings, CLFlags )
   implicit none
      ! Storing the arguments
   type( ProgDesc ),                   intent(in   )  :: ProgInfo
   type( ADskDriver_Settings ),        intent(  out)  :: CLSettings        !< Command line arguments passed in
   type( ADskDriver_Flags ),           intent(  out)  :: CLFlags           !< Flags indicating which command line arguments were specified

      ! Set some CLSettings to null/default values
   CLSettings%DvrIptFileName   =  ""             ! No input name name until set
   CLSettings%ADskIptFileName  =  ""             ! No ADsk input file name until set
   CLSettings%NumTimeSteps     =  0_IntKi
   CLSettings%DT               =  0.0_DbKi
   CLSettings%TStart           =  0.0_ReKi
   CLSettings%ProgInfo         =  ProgInfo       ! Driver info

      ! Set some CLFlags to null/default values
   CLFlags%DvrIptFile          =  .FALSE.        ! Driver     input filename given as command line argument
   CLFlags%ADskIptFile         =  .FALSE.        ! AeroDisk input filename given as command line argument
   CLFlags%TStart              =  .FALSE.        ! specified time to start at
   CLFlags%NumTimeSteps        =  .FALSE.        ! specified a number of timesteps
   CLFlags%NumTimeStepsDefault =  .FALSE.        ! specified 'DEFAULT' for number of timesteps
   CLFlags%DT                  =  .FALSE.        ! specified a resolution in time
   CLFlags%DTDefault           =  .FALSE.        ! specified 'DEFAULT' for resolution in time
   CLFlags%Verbose             =  .FALSE.        ! Turn on verbose error reporting?
   CLFlags%VVerbose            =  .FALSE.        ! Turn on very verbose error reporting?

end subroutine InitSettingsFlags

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!> This subroutine retrieves the command line arguments and passes them to the
!! AeroDisk_driver_subs::parsearg routine for processing.
SUBROUTINE RetrieveArgs( CLSettings, CLFlags, ErrStat, ErrMsg )
      ! Storing the arguments
   type( ADskDriver_Flags ),           intent(  out)  :: CLFlags      !< Flags indicating which command line arguments were specified
   type( ADskDriver_Settings ),        intent(  out)  :: CLSettings           !< Command line arguments passed in
   integer(IntKi),                     intent(  out)  :: ErrStat
   CHARACTER(*),                       intent(  out)  :: ErrMsg

      ! Local variable
   integer(IntKi)                                     :: i                    !< Generic counter
   character(1024)                                    :: Arg                  !< argument given
   character(1024)                                    :: ArgUC                !< Upper case argument to check
   integer(IntKi)                                     :: NumInputArgs         !< Number of argements passed in from command line
   logical                                            :: adskFlag             !< The -adsk flag was set
   character(1024)                                    :: FileName             !< Filename from the command line.
   logical                                            :: FileNameGiven        !< Flag indicating if a filename was given.
   integer(IntKi)                                     :: ErrStatTmp           !< Temporary error status (for calls)
   character(1024)                                    :: ErrMsgTmp            !< Temporary error message (for calls)

      ! initialize some things
   CLFlags%DvrIptFile =  .FALSE.
   ErrStat            =  ErrID_None
   ErrStatTmp         =  ErrID_None
   ErrMsg             =  ''
   ErrMsgTmp          =  ''
   adskFlag           =  .FALSE.
   FileNameGiven      =  .FALSE.
   FileName           =  ''

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
            CALL DispHelpText()
            CALL ProgExit(0)
         ENDIF


            ! Check the argument and put it where it belongs
            ! chop the SwChar off before passing the argument
         CALL ParseArg( CLSettings, CLFlags, ArgUC(2:), Arg(2:), adskFlag, ErrStatTmp, ErrMsgTmp )
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

      ! Was the -adsk flag set?  If so, the filename is the AeroDisk input file.  Otherwise
      ! it is the driver input file.
   IF ( adskFlag ) THEN
      CLSettings%ADskIptFileName  =  TRIM(FileName)
      CLFlags%ADskIptFile =  .TRUE.
      call GetRoot( CLSettings%ADskIptFileName, CLSettings%OutRootName )
      CLFlags%OutRootName =  .TRUE.
   ELSE
      CLSettings%DvrIptFileName  =  TRIM(FileName)
      CLFlags%DvrIptFile =  .TRUE.
   ENDIF



   !-------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
      CONTAINS


   !-------------------------------------------------------------------------------
   !> Convert a string to a real number
   FUNCTION StringToReal( StringIn, ErrStat )
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(IN   )  :: StringIn

      REAL(ReKi)                                         :: StringToReal
      INTEGER(IntKi)                                     :: ErrStatTmp         ! Temporary variable to hold the error status

         read( StringIn, *, iostat=ErrStatTmp) StringToReal

            ! If that isn't a number, only warn since we can continue by skipping this value
         IF ( ErrStatTmp .ne. 0 ) ErrStat  = ErrID_Warn

   END FUNCTION StringToReal



   !-------------------------------------------------------------------------------
   SUBROUTINE ParseArg( CLSettings, CLFlags, ThisArgUC, ThisArg, adskFlagSet, ErrStat, ErrMsg )
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
         ! Parse and store the input argument  !
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

      USE NWTC_Library
      USE AeroDisk_Driver_Types
      USE AeroDisk_Types

      IMPLICIT NONE

         ! Storing the arguments
      TYPE( ADskDriver_Flags ),            INTENT(INOUT)  :: CLFlags             ! Flags indicating which arguments were specified
      TYPE( ADskDriver_Settings ),         INTENT(INOUT)  :: CLSettings          ! Arguments passed in

      CHARACTER(*),                       INTENT(IN   )  :: ThisArgUC            ! The current argument (upper case for testing)
      CHARACTER(*),                       INTENT(IN   )  :: ThisArg              ! The current argument (as passed in for error messages)
      LOGICAL,                            INTENT(INOUT)  :: adskFlagSet          ! Was the -adsk flag given?

         ! Error Handling
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


         ! local variables
      INTEGER(IntKi)                                     :: Delim1               ! where the [ is
      INTEGER(IntKi)                                     :: Delim2               ! where the ] is
      INTEGER(IntKi)                                     :: DelimSep             ! where the : is
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
            ! check to see if the filename is the name of the ADsk input file
         IF       ( ThisArgUC(1:4) == "ADSK" )   THEN
            adskFlagSet              = .TRUE.             ! More logic in the routine that calls this one to set things.
            RETURN
         ELSEIF   ( ThisArgUC(1:2) == "VV"  )   THEN
            CLFlags%VVerbose        = .TRUE.
            RETURN
         ELSEIF   ( ThisArgUC(1:1) == "V"   )   THEN
            CLFlags%Verbose         = .TRUE.
            RETURN
         ELSE
            CALL SetErrStat( ErrID_Warn," Unrecognized option '"//SwChar//TRIM(ThisArg)//"'. Ignoring. Use option "//SwChar//"help for list of options.",  &
               ErrStat,ErrMsg,'ParseArg')
         ENDIF

      ENDIF


         ! "DT[#]"
      IF( ThisArgUC(1:Delim1) == "DT["      ) THEN
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
!FIXME: add in the other inputs here.

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
SUBROUTINE ParseDvrIptFile( DvrFileName, DvrFileInfo, DvrFlags, DvrSettings, ProgInfo, CaseTime, CaseData, ErrStat, ErrMsg )

   CHARACTER(1024),                    INTENT(IN   )  :: DvrFileName
   type(FileInfoType),                 INTENT(IN   )  :: DvrFileInfo          ! Input file stored in FileInfoType structure
   TYPE(ADskDriver_Flags),             INTENT(INOUT)  :: DvrFlags
   TYPE(ADskDriver_Settings),          INTENT(INOUT)  :: DvrSettings
   TYPE(ProgDesc),                     INTENT(IN   )  :: ProgInfo
   real(DbKi),          allocatable,   intent(  out)  :: CaseTime(:)
   real(ReKi),          allocatable,   intent(  out)  :: CaseData(:,:)
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              ! returns a non-zero value when an error occurs
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)                                     :: CurLine              ! Current line in parsing
   INTEGER(IntKi)                                     :: TabLines             ! Number of lines in the table
   integer(IntKi)                                     :: i                    !< generic loop counter
   real(DbKi)                                         :: TmpDb7(7)            !< temporary real array
   CHARACTER(1024)                                    :: RootName             ! Root name of AeroDisk driver input file

      ! Input file echoing
   LOGICAL                                            :: EchoFileContents     ! Do we echo the driver file out or not?
   INTEGER(IntKi)                                     :: UnEc          ! The local unit number for this module's echo file
   CHARACTER(1024)                                    :: EchoFileName         ! Name of AeroDisk driver echo file

      ! Time steps
   CHARACTER(1024)                                    :: InputChr             ! Character string for timesteps and input file names (to handle DEFAULT or NONE value)

      ! Local error handling
   INTEGER(IntKi)                                     :: ios                  !< I/O status
   INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp            !< Temporary error messages for calls


      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEc = -1

   call GetRoot( DvrFileName, RootName )

   !======  General  ====================================================================================
   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( DvrFileInfo, CurLine, 'Echo', EchoFileContents, ErrStatTmp, ErrMsgTmp )
      if (Failed()) return;

   if ( EchoFileContents ) then
      CALL OpenEcho ( UnEc, TRIM(RootName)//'.ech', ErrStatTmp, ErrMsgTmp )
         if (Failed()) return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDisk driver input file: '//trim(DvrFileName)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') DvrFileInfo%Lines(1)
      WRITE(UnEc, '(A)') DvrFileInfo%Lines(2)
      WRITE(UnEc, '(A)') DvrFileInfo%Lines(3)

      CurLine = 4
      call ParseVar( DvrFileInfo, CurLine, 'Echo', EchoFileContents, ErrStatTmp, ErrMsgTmp, UnEc )
         if (Failed()) return
   endif

   !======  Primary file and rootname  ==================================================================
   if ( EchoFileContents )   WRITE(UnEc, '(A)') DvrFileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! AeroDisk input file
   call ParseVar( DvrFileInfo, CurLine, "ADskIptFile", DvrSettings%ADskIptFileName, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%ADskIptFile  =  .TRUE.

      ! AeroDisk output root name
   call ParseVar( DvrFileInfo, CurLine, "OutRootName", DvrSettings%OutRootName, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%OutRootName  =  .TRUE.
 

   !======  Geometry and Environment  ====================================================================
   if ( EchoFileContents )   WRITE(UnEc, '(A)') DvrFileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! AirDens
   call ParseVar( DvrFileInfo, CurLine, "AirDens", DvrSettings%AirDens, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%AirDens   =  .TRUE.

      ! RotorRad
   call ParseVar( DvrFileInfo, CurLine, "RotorRad", DvrSettings%RotorRad, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%RotorRad   =  .TRUE.

      ! RotorHeight
   call ParseVar( DvrFileInfo, CurLine, "RotorHeight", DvrSettings%RotorHeight, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%RotorHeight   =  .TRUE.

      ! ShftTilt
   call ParseVar( DvrFileInfo, CurLine, "ShftTilt", DvrSettings%ShftTilt, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%ShftTilt   =  .TRUE.


   !======  Case analysis  ========== ====================================================================
   if ( EchoFileContents )   WRITE(UnEc, '(A)') DvrFileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! TStart    -- start time
   call ParseVar( DvrFileInfo, CurLine, "TStart", DvrSettings%TStart, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 
   DvrFlags%TStart   =  .TRUE.


      ! DT    -- Timestep size for the driver to take (or DEFAULT for what the file contains)
   call ParseVar( DvrFileInfo, CurLine, "DT", InputChr, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return 

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( InputChr )
   IF ( TRIM(InputChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%DT         =  .TRUE.
      DvrFlags%DTDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.
      READ (InputChr,*,IOSTAT=IOS)   DvrSettings%DT
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'DT',NumType, ErrStatTmp, ErrMsgTmp )
         if (Failed()) return
      ELSE     ! Was ok, so set the flags
         DvrFlags%DT         =  .TRUE.
         DvrFlags%DTDefault  =  .FALSE.
      ENDIF
   ENDIF


      ! Number of timesteps
   call ParseVar( DvrFileInfo, CurLine, "NumTimeSteps", InputChr, ErrStatTmp, ErrMsgTmp, UnEc )
   if (Failed()) return

      ! Check if we asked for the DEFAULT (use what is in the file)
   CALL Conv2UC( InputChr )
   IF ( TRIM(InputChr) == 'DEFAULT' ) THEN     ! we asked for the default value
      DvrFlags%NumTimeSteps         =  .FALSE.
      DvrFlags%NumTimeStepsDefault  =  .TRUE.         ! This flag tells us to use the inflow wind file values
   ELSE
         !  We probably have a number if it isn't 'DEFAULT', so do an internal read and check to
         !  make sure that it was appropriately interpretted.
      READ (InputChr,*,IOSTAT=IOS)   DvrSettings%NumTimeSteps
      IF ( IOS /= 0 )  THEN  ! problem in the read, so parse the error.
         CALL CheckIOS ( IOS, '', 'NumTimeSteps',NumType, ErrStatTmp, ErrMsgTmp )
         if (Failed()) return
      ELSE     ! Was ok, so set the flags
         DvrFlags%NumTimeSteps         =  .TRUE.
         DvrFlags%NumTimeStepsDefault  =  .FALSE.
      ENDIF
   ENDIF


      ! Column headers
   if ( EchoFileContents )   WRITE(UnEc, '(A)') DvrFileInfo%Lines(CurLine)
   CurLine = CurLine + 1
   if ( EchoFileContents )   WRITE(UnEc, '(A)') DvrFileInfo%Lines(CurLine)
   CurLine = CurLine + 1


      ! Last line of table is assumed to be last line in file (or included file)
   TabLines = DvrFileInfo%NumLines - CurLine + 1
   call AllocAry( CaseTime,    TabLines, 'CaseTime', ErrStatTmp, ErrMsgTmp ); if (Failed()) return;
   call AllocAry( CaseData, 6, TabLines, 'CaseData', ErrStatTmp, ErrMsgTmp ); if (Failed()) return;
   do i=1,Tablines
      call ParseAry ( DvrFileInfo, CurLine, 'Coordinates', TmpDb7, 7, ErrStatTmp, ErrMsgTmp, UnEc )
         if (Failed())  return;
      ! Set Time_(s)
      CaseTime(i)     = TmpDb7(1)
      ! Set time, wind_x, wind_y, wind_z
      CaseData(1:3,i) = real(TmpDb7(2:4),ReKi)
      ! Set RotSpeed    (rpm -> rad/s)
      CaseData(4,i)   = real(TmpDb7(5),ReKi) * Pi / 30.0_ReKi
      ! Set Pitch       (deg -> rad)
      CaseData(5,i)   = real(TmpDb7(6),ReKi) * D2R
      ! Set Yaw         (deg -> rad)
      CaseData(6,i)   = real(TmpDb7(7),ReKi) * D2R
   enddo


      ! Close the echo and input file
   CALL CleanupEchoFile( EchoFileContents, UnEc )


CONTAINS

   !> Set error status, close stuff, and return
   logical function Failed()
      CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,'ParseDvrIptFile')
      if (ErrStat >= AbortErrLev) then
         CALL CleanupEchoFile( EchoFileContents, UnEc )
      endif
      Failed =    ErrStat >= AbortErrLev
   end function Failed

   !> Clean up the module echo file
   subroutine CleanupEchoFile( EchoFlag, UnEcho)
      logical,                      intent(in   )  :: EchoFlag          ! local version of echo flag
      integer(IntKi),               intent(in   )  :: UnEcho            !  echo unit number

         ! Close this module's echo file
      if ( EchoFlag ) then
         close(UnEcho)
      endif
   END SUBROUTINE CleanupEchoFile

END SUBROUTINE ParseDvrIptFile


!> This subroutine copies an command line (CL) settings over to the program settings.  Warnings are
!! issued if anything is changed from what the driver input file requested.
SUBROUTINE UpdateSettingsWithCL( DvrFlags, DvrSettings, CLFlags, CLSettings, DVRIPT, ErrStat, ErrMsg )

   TYPE(ADskDriver_Flags),             INTENT(INOUT)  :: DvrFlags
   TYPE(ADskDriver_Settings),          INTENT(INOUT)  :: DvrSettings
   TYPE(ADskDriver_Flags),             INTENT(IN   )  :: CLFlags
   TYPE(ADskDriver_Settings),          INTENT(IN   )  :: CLSettings
   LOGICAL,                            INTENT(IN   )  :: DVRIPT
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg


      ! Local variables
   INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status for calls
   CHARACTER(1024)                                    :: ErrMsgTmp         !< Temporary error status for calls
   LOGICAL                                            :: WindGridModify    !< Did we modify any of the WindGrid related settings?
   character(*), parameter                            :: RoutineName = 'UpdateSettingsWithCL'

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

      !  Check TStart
   IF ( CLFlags%TStart ) THEN
      IF ( DvrFlags%TStart .AND. ( .NOT. EqualRealNos(DvrSettings%TStart, CLSettings%TStart) ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for TStart with '//TRIM(Num2LStr(CLSettings%TStart))//'.', &
            ErrStat,ErrMsg,RoutineName)
      ELSE
         DvrFlags%TStart   =  .TRUE.
      ENDIF
      DvrSettings%TStart   =  CLSettings%TStart
   ENDIF

      ! Check DT
   IF ( CLFlags%DT ) THEN
      IF ( DvrFlags%DT .AND. ( .NOT. EqualRealNos(DvrSettings%DT, CLSettings%DT) ) )   THEN
         CALL SetErrStat( ErrID_Warn, ' Overriding driver input value for DT with '//TRIM(Num2LStr(CLSettings%DT))//'.',  &
            ErrStat,ErrMsg,RoutineName)
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
            ErrStat,ErrMsg,RoutineName)
      ELSE
         DvrFlags%NumTimeSteps      =  .TRUE.
      ENDIF
      DvrSettings%NumTimeSteps      =  CLSettings%NumTimeSteps
      DvrFlags%NumTimeStepsDefault  =  .FALSE.
   ENDIF

      ! Make sure there is at least one timestep
   DvrSettings%NumTimeSteps   =  MAX(DvrSettings%NumTimeSteps,1_IntKi)


      !--------------------------------------------
      ! If there was no driver input file, we need to set a few things.
      !--------------------------------------------

   IF ( .NOT. DVRIPT ) THEN

         ! Do we need to set the NumTimeStepsDefault flag?
      IF ( .NOT. DvrFlags%NumTimeSteps )  THEN
         DvrFlags%NumTimeStepsDefault  =  .TRUE.
         CALL SetErrStat( ErrID_Info,' The number of timesteps is not specified.  Defaulting to what is in the input series file.', &
            ErrStat,ErrMsg,RoutineName)
      ENDIF
   ENDIF


!FIXME: remove this after parsing rest of input file.
      ! If no DT value has been set (DEFAULT requested), we need to set a default to pass into ADsk
   IF ( .NOT. DvrFlags%DT ) THEN
      DvrSettings%DT =  0.025_DbKi     ! This value gets passed into the ADsk_Init routine, so something must be set.
   ENDIF


END SUBROUTINE UpdateSettingsWithCL



!> This routine exists only to support the development of the module.  It will not be needed after the module is complete.
SUBROUTINE  printSettings( DvrFlags, DvrSettings )
      ! The arguments
   TYPE( ADskDriver_Flags ),           INTENT(IN   )  :: DvrFlags           !< Flags indicating which settings were set
   TYPE( ADskDriver_Settings ),        INTENT(IN   )  :: DvrSettings        !< Stored settings

   CALL WrsCr(TRIM(GetNVD(DvrSettings%ProgInfo)))
   CALL WrScr(' DvrIptFile:          '//FLAG(DvrFlags%DvrIptFile)//        '      '//TRIM(DvrSettings%DvrIptFileName))
   CALL WrScr(' ADskIptFile:         '//FLAG(DvrFlags%ADskIptFile)//       '      '//TRIM(DvrSettings%ADskIptFileName))
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
   CALL WrScr(' Verbose:             '//FLAG(DvrFlags%Verbose))
   CALL WrScr(' VVerbose:            '//FLAG(DvrFlags%VVerbose))
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


SUBROUTINE Dvr_InitializeOutputFile(OutUnit,IntOutput,RootName,ErrStat,ErrMsg)
   integer(IntKi),              intent(  out):: OutUnit
   type(ADsk_InitOutputType),   intent(in   ):: IntOutput     ! Output for initialization routine
   integer(IntKi),              intent(  out):: ErrStat     ! Error status of the operation
   character(*),                intent(  out):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   character(*),                intent(in   ):: RootName
   integer(IntKi)                            :: i
   integer(IntKi)                            :: numOuts
   integer(IntKi)                            :: ErrStat2                     ! Temporary Error status
   character(ErrMsgLen)                      :: ErrMsg2                      ! Temporary Error message
   character(*), parameter                   :: RoutineName = 'Dvr_InitializeOutputFile'

   ErrStat = ErrID_none
   ErrMsg  = ""

   CALL GetNewUnit(OutUnit,ErrStat2,ErrMsg2)
   CALL OpenFOutFile ( OutUnit, trim(RootName)//'.out', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   write (OutUnit,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(IntOutput%Ver))
   write (OutUnit,'()' )    !print a blank line

   numOuts = size(IntOutput%WriteOutputHdr)
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................

   write (OutUnit,'()')
   write (OutUnit,'()')
   write (OutUnit,'()')

   call WrFileNR ( OutUnit, 'Time' )

   do i=1,NumOuts
      call WrFileNR ( OutUnit, tab//IntOutput%WriteOutputHdr(i) )
   end do ! i

   write (OutUnit,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................

   call WrFileNR ( OutUnit, '(s)' )

   do i=1,NumOuts
      call WrFileNR ( Outunit, tab//trim(IntOutput%WriteOutputUnt(i)) )
   end do ! i

   write (OutUnit,'()')


END SUBROUTINE Dvr_InitializeOutputFile

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Dvr_WriteOutputLine(t,OutUnit, OutFmt, Output)
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   integer(IntKi)         ,  intent(in   )   :: OutUnit              ! Status of error message
   character(*)           ,  intent(in   )   :: OutFmt
   type(ADsk_OutputType),    intent(in   )   :: Output
   integer(IntKi)                            :: errStat              ! Status of error message (we're going to ignore errors in writing to the file)
   character(ErrMsgLen)                      :: errMsg               ! Error message if ErrStat /= ErrID_None
   character(200)                            :: frmt                 ! A string to hold a format specifier
   character(15)                             :: tmpStr               ! temporary string to print the time output as text

   frmt = '"'//tab//'"'//trim(OutFmt)      ! format for array elements from individual modules

      ! time
   write( tmpStr, '(F15.6)' ) t
   call WrFileNR( OutUnit, tmpStr )
   call WrNumAryFileNR ( OutUnit, Output%WriteOutput,  frmt, errStat, errMsg )

     ! write a new line (advance to the next line)
   write (OutUnit,'()')
end subroutine Dvr_WriteOutputLine


END MODULE AeroDisk_Driver_Subs
